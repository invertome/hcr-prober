# hcr_prober/blast_wrapper.py
import subprocess, os, pandas as pd, tempfile
from loguru import logger

def create_blast_db(ref_fasta, db_path):
    """Creates a BLAST database from a FASTA file if it doesn't exist or is outdated."""
    if not ref_fasta: return None
    if not db_path: db_path = os.path.dirname(ref_fasta)
    db_name = os.path.join(db_path, os.path.basename(ref_fasta))
    if os.path.exists(f'{db_name}.nsq'):
        try:
            if os.path.getmtime(ref_fasta) < os.path.getmtime(f'{db_name}.nsq'): return db_name
        except FileNotFoundError: pass
    logger.info(f'Creating/updating BLAST database for {ref_fasta}...')
    os.makedirs(db_path, exist_ok=True)
    cmd = ['makeblastdb', '-in', ref_fasta, '-dbtype', 'nucl', '-out', db_name, '-title', os.path.basename(ref_fasta)]
    try: subprocess.run(cmd, check=True, capture_output=True, text=True)
    except Exception as e: logger.critical(f'FATAL: Failed to create BLAST DB. Error: {e}'); exit(1)
    return db_name

def _run_blast(probes, db_name, temp_dir, extra_args):
    """Runs BLASTn with a given set of probes against a database."""
    if not probes: return None
    # BLAST the combined 52bp target sequence (DN_arm + 'nn' + UP_arm)
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
    """
    Filters probes based on BLAST results against a reference database.
    Implements three strategies for positive selection:
    1. any-strong-hit: Keeps probes with any hit meeting score thresholds.
    2. most-likely-transcript: Identifies the most-hit transcript and keeps probes unique to it.
    3. specific-transcript: Keeps probes unique to a user-specified transcript ID.
    """
    if not probes or not args.blast_ref: return probes, {}

    blast_path = _run_blast(probes, args.blast_db_positive, temp_dir, args.blast_extra_args)
    blast_report = {'strong_hits': pd.DataFrame()}

    if not blast_path: return [], {'positive': blast_report}
    try:
        results = pd.read_csv(blast_path, sep='\t', names=['pair_id', 'hit_id', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'], dtype=str)
        results[['pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']] = results[['pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']].apply(pd.to_numeric)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        if os.path.exists(blast_path): os.remove(blast_path)
        return [], {'positive': blast_report}
    if os.path.exists(blast_path): os.remove(blast_path)

    # Step 1: Find all high-quality hits based on bitscore and e-value.
    # This dataframe is used for all strategies and for the final report.
    strong_hits = results[(results['bitscore'] >= args.min_bitscore) & (results['evalue'] <= args.max_evalue)].copy()
    blast_report['strong_hits'] = strong_hits

    if strong_hits.empty:
        logger.warning(f'({args.job_name}) No probes passed the initial BLAST filter (bitscore >= {args.min_bitscore}, evalue <= {args.max_evalue}).')
        return [], {'positive': blast_report}

    # Step 2: Apply the selected filtering strategy.
    strategy = args.positive_selection_strategy
    logger.info(f'Applying positive selection strategy: "{strategy}"')
    passed_probes_set = set()

    if strategy == 'any-strong-hit':
        # Default behavior: keep any probe that has at least one strong hit.
        passed_probes_set = set(strong_hits['pair_id'])
    
    elif strategy == 'specific-transcript':
        # Keep probes that uniquely hit a user-defined transcript ID.
        target_id = args.target_transcript_id
        logger.info(f'Filtering for probes that uniquely hit target ID: "{target_id}"')
        # Find all probes that hit our specified target.
        hits_on_target = strong_hits[strong_hits['hit_id'] == target_id]
        probes_hitting_target = set(hits_on_target['pair_id'])
        # Find all probes that hit any other transcript.
        hits_off_target = strong_hits[strong_hits['hit_id'] != target_id]
        probes_hitting_others = set(hits_off_target['pair_id'])
        # The final set is probes that hit the target AND DO NOT hit any others.
        passed_probes_set = probes_hitting_target - probes_hitting_others

    elif strategy == 'most-likely-transcript':
        # Find the single transcript with the most hits, and keep probes unique to it.
        hit_counts = strong_hits['hit_id'].value_counts()
        if hit_counts.empty:
            logger.warning('No strong hits found to determine the most likely transcript.')
            passed_probes_set = set()
        else:
            most_likely_target_id = hit_counts.index[0]
            logger.info(f'Identified "{most_likely_target_id}" as the most likely transcript (with {hit_counts.iloc[0]} strong hits).')
            # Find all probes that hit this most-likely target.
            hits_on_target = strong_hits[strong_hits['hit_id'] == most_likely_target_id]
            probes_hitting_target = set(hits_on_target['pair_id'])
            # Find all probes that hit any other transcript.
            hits_off_target = strong_hits[strong_hits['hit_id'] != most_likely_target_id]
            probes_hitting_others = set(hits_off_target['pair_id'])
            # The final set is probes that hit the most likely target AND DO NOT hit any others.
            passed_probes_set = probes_hitting_target - probes_hitting_others
            
    # Step 3: Filter the original probe list based on the determined set of passing probe IDs.
    final_probes = [p for p in probes if p['pair_id'] in passed_probes_set]
    logger.info(f'({args.job_name}) {len(final_probes)} of {len(probes)} pairs passed POSITIVE screen using "{strategy}" strategy.')
    return final_probes, {'positive': blast_report}