# hcr_prober/blast_wrapper.py
import subprocess, os, pandas as pd, tempfile
from loguru import logger

def create_blast_db(ref_fasta, db_path):
    """Creates a BLAST database from a FASTA file if it doesn't exist or is outdated."""
    if not ref_fasta: return None
    db_name_prefix = os.path.splitext(os.path.basename(ref_fasta))[0]
    db_dir = db_path if db_path else tempfile.gettempdir()
    db_name = os.path.join(db_dir, db_name_prefix)
    # A simple check using one of the BLAST db files to see if it's outdated.
    db_check_file = f'{db_name}.nsq'
    if os.path.exists(db_check_file) and os.path.getmtime(ref_fasta) < os.path.getmtime(db_check_file):
        return db_name
    logger.info(f'Creating/updating BLAST database for {os.path.basename(ref_fasta)}...')
    os.makedirs(os.path.dirname(db_name), exist_ok=True)
    cmd = ['makeblastdb', '-in', ref_fasta, '-dbtype', 'nucl', '-out', db_name, '-title', db_name_prefix]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except Exception as e:
        logger.critical(f'FATAL: Failed to create BLAST DB. Error: {e.stderr if hasattr(e, 'stderr') else e}')
        exit(1)
    return db_name

def _run_blast(probes, db_name, temp_dir, extra_args):
    """Runs BLASTn with a given set of probes against a database."""
    if not probes: return None
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta', dir=temp_dir) as f:
        for p in probes:
            f.write(f'>{p['pair_id']}\n{p['probe_dn_target']}nn{p['probe_up_target']}\n')
        query_path = f.name
    blast_out_path = f'{query_path}.blast.tsv'
    cmd = ['blastn', '-query', query_path, '-db', db_name, '-out', blast_out_path, '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore', '-task', 'blastn-short']
    if extra_args: cmd.extend(extra_args)
    try: subprocess.run(cmd, check=True, capture_output=True, text=True)
    except Exception as e: logger.error(f'BLAST search failed. Error: {e}'); return None
    finally: os.remove(query_path)
    return blast_out_path

def filter_probes_by_blast(probes, args, temp_dir):
    """Filters probes based on BLAST results, implementing various strategies."""
    if not probes or not args.blast_ref: return probes, {}

    # Step 1: Broad discovery BLAST search.
    blast_path = _run_blast(probes, args.blast_db_positive, temp_dir, args.blast_extra_args)
    blast_report = {'strong_hits': pd.DataFrame()}
    if not blast_path: return [], {'positive': blast_report}
    try:
        results = pd.read_csv(blast_path, sep='\t', names=['pair_id', 'hit_id', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'], dtype=str)
        results[['pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']] = results[['pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']].apply(pd.to_numeric)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        if os.path.exists(blast_path): os.remove(blast_path)
        return [], {'positive': blast_report}
    finally:
        if os.path.exists(blast_path): os.remove(blast_path)

    # Find all plausible hits using user-defined thresholds.
    plausible_hits = results[(results['bitscore'] >= args.min_bitscore) & (results['evalue'] <= args.max_evalue)].copy()
    blast_report['strong_hits'] = plausible_hits # We will report all plausible hits for transparency.

    if plausible_hits.empty:
        logger.warning(f'({args.job_name}) No probes passed the initial BLAST filter (bitscore >= {args.min_bitscore}, evalue <= {args.max_evalue}).')
        return [], {'positive': blast_report}

    # Step 2: Apply the selected filtering strategy.
    strategy = args.positive_selection_strategy
    logger.info(f'Applying positive selection strategy: "{strategy}"')
    passed_probes_set = set()

    if strategy == 'any-strong-hit':
        passed_probes_set = set(plausible_hits['pair_id'])

    elif strategy == 'specific-id':
        target_id = args.target_transcript_id
        logger.info(f'Filtering for probes that uniquely hit target ID: "{target_id}"')
        hits_on_target = plausible_hits[plausible_hits['hit_id'] == target_id]
        probes_hitting_target = set(hits_on_target['pair_id'])
        hits_off_target = plausible_hits[plausible_hits['hit_id'] != target_id]
        probes_hitting_others = set(hits_off_target['pair_id'])
        passed_probes_set = probes_hitting_target - probes_hitting_others

    elif strategy == 'best-coverage':
        # This strategy identifies the best transcript based on coverage and quality, then selects unique probes for it.
        # A) Calculate scores for each transcript.
        transcript_scores = []
        for hit_id, group in plausible_hits.groupby('hit_id'):
            # Breadth of coverage: how many unique probes hit this transcript?
            breadth = group['pair_id'].nunique()
            # Average quality: what is the mean bitscore of these hits?
            avg_quality = group['bitscore'].mean()
            transcript_scores.append({'id': hit_id, 'breadth': breadth, 'quality': avg_quality})
        
        if not transcript_scores:
            logger.warning('No plausible hits found to determine the best-supported transcript.')
            passed_probes_set = set()
        else:
            # B) Rank transcripts: primary key is breadth, secondary is quality.
            ranked_transcripts = sorted(transcript_scores, key=lambda x: (x['breadth'], x['quality']), reverse=True)
            best_transcript_id = ranked_transcripts[0]['id']
            logger.info(f'Identified "{best_transcript_id}" as best-supported transcript (Breadth: {ranked_transcripts[0]['breadth']}, Avg. Quality: {ranked_transcripts[0]['quality']:.2f}).')

            # C) Final strict selection: find probes unique to the best transcript using a fixed high-quality threshold.
            HIGH_QUALITY_BITSCORE = 75.0  # Internal threshold for final selection, independent of user args.
            strong_hits = plausible_hits[plausible_hits['bitscore'] >= HIGH_QUALITY_BITSCORE]
            
            hits_on_best = strong_hits[strong_hits['hit_id'] == best_transcript_id]
            probes_on_best = set(hits_on_best['pair_id'])
            
            hits_off_best = strong_hits[strong_hits['hit_id'] != best_transcript_id]
            probes_off_best = set(hits_off_best['pair_id'])
            
            passed_probes_set = probes_on_best - probes_off_best
            
    # Step 3: Filter the original probe list.
    final_probes = [p for p in probes if p['pair_id'] in passed_probes_set]
    logger.info(f'({args.job_name}) {len(final_probes)} of {len(probes)} pairs passed POSITIVE screen using "{strategy}" strategy.')
    return final_probes, {'positive': blast_report}