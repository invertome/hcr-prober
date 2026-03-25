# hcr_prober/blast_wrapper.py
import subprocess, os, sys, pandas as pd, tempfile
from loguru import logger
def create_blast_db(ref_fasta, db_path):
    if not ref_fasta: return None
    db_name_prefix, db_dir = os.path.splitext(os.path.basename(ref_fasta))[0], db_path if db_path else tempfile.gettempdir()
    db_name = os.path.join(db_dir, db_name_prefix); db_check_file = f'{db_name}.nsq'
    if os.path.exists(db_check_file) and os.path.getmtime(ref_fasta) < os.path.getmtime(db_check_file): return db_name
    logger.info(f'Creating/updating BLAST database for {os.path.basename(ref_fasta)}...')
    os.makedirs(os.path.dirname(db_name), exist_ok=True)
    cmd = ['makeblastdb', '-in', ref_fasta, '-dbtype', 'nucl', '-out', db_name, '-title', db_name_prefix]
    try: subprocess.run(cmd, check=True, capture_output=True, text=True)
    except Exception as e: logger.critical(f'FATAL: Failed to create BLAST DB. Error: {e.stderr if hasattr(e, "stderr") else e}'); sys.exit(1)
    return db_name
def _run_blast(probes, db_name, temp_dir, extra_args):
    if not probes: return None
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta', dir=temp_dir) as f:
        for p in probes: f.write(f'>{p['pair_id']}\n{p['probe_dn_target']}nn{p['probe_up_target']}\n')
        query_path = f.name
    blast_out_path = f'{query_path}.blast.tsv'
    cmd = ['blastn', '-query', query_path, '-db', db_name, '-out', blast_out_path, '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore', '-task', 'blastn-short']
    if extra_args: cmd.extend(extra_args)
    try: subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logger.error(f'BLAST search failed with return code {e.returncode}. stderr: {e.stderr}')
        return None
    except Exception as e:
        logger.error(f'Unexpected error during BLAST: {e}')
        return None
    finally: os.remove(query_path)
    return blast_out_path
def filter_probes_by_blast(probes, args, temp_dir):
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
    finally: os.remove(blast_path)
    plausible_hits = results[(results['bitscore'] >= args.min_bitscore) & (results['evalue'] <= args.max_evalue)].copy()
    blast_report['strong_hits'] = plausible_hits
    if plausible_hits.empty: return [], {'positive': blast_report}
    strategy = args.positive_selection_strategy
    logger.info(f'Applying positive selection strategy: "{strategy}"')
    passed_probes_set = set()
    if strategy == 'any-strong-hit': passed_probes_set = set(plausible_hits['pair_id'])
    elif strategy == 'specific-id':
        target_id = args.target_transcript_id
        hits_on_target = plausible_hits[plausible_hits['hit_id'] == target_id]
        probes_hitting_target = set(hits_on_target['pair_id'])
        hits_off_target = plausible_hits[plausible_hits['hit_id'] != target_id]
        passed_probes_set = probes_hitting_target - set(hits_off_target['pair_id'])
    elif strategy == 'best-coverage':
        transcript_scores = []
        for hit_id, group in plausible_hits.groupby('hit_id'):
            transcript_scores.append({'id': hit_id, 'breadth': group['pair_id'].nunique(), 'quality': group['bitscore'].mean()})
        if not transcript_scores: passed_probes_set = set()
        else:
            ranked_transcripts = sorted(transcript_scores, key=lambda x: (x['breadth'], x['quality']), reverse=True)
            best_transcript_id = ranked_transcripts[0]['id']
            logger.info(f'Identified "{best_transcript_id}" as best-supported transcript (Breadth: {ranked_transcripts[0]['breadth']}, Avg. Quality: {ranked_transcripts[0]['quality']:.2f}).')
            strong_hits = plausible_hits[plausible_hits['bitscore'] >= args.min_bitscore]
            hits_on_best = strong_hits[strong_hits['hit_id'] == best_transcript_id]
            probes_on_best = set(hits_on_best['pair_id'])
            hits_off_best = strong_hits[strong_hits['hit_id'] != best_transcript_id]
            passed_probes_set = probes_on_best - set(hits_off_best['pair_id'])
    probe_dict = {p['pair_id']: p for p in probes}
    final_probes = [probe_dict[pid] for pid in passed_probes_set if pid in probe_dict]
    logger.info(f'({args.job_name}) {len(final_probes)} of {len(probes)} total candidates passed the POSITIVE screen.')
    return final_probes, {'positive': blast_report}
def run_negative_screen(probes, args, temp_dir, target_ids):
    """Screen probes against non-target transcripts. Reject probes with strong off-target hits."""
    if not probes:
        return probes, {}
    neg_bitscore = args.negative_bitscore if args.negative_bitscore is not None else args.min_bitscore
    neg_evalue = args.negative_evalue if args.negative_evalue is not None else args.max_evalue
    # Determine negative reference
    neg_ref = getattr(args, 'blast_negative_ref', None)
    if not neg_ref and args.blast_ref:
        # Auto-derive: extract non-target sequences from blast_ref
        from Bio import SeqIO
        neg_path = os.path.join(temp_dir, 'negative_ref.fasta')
        count = 0
        with open(neg_path, 'w') as f:
            for rec in SeqIO.parse(args.blast_ref, 'fasta'):
                if rec.id not in target_ids:
                    f.write(f'>{rec.id}\n{str(rec.seq)}\n')
                    count += 1
        if count == 0:
            logger.warning('No non-target sequences found for negative screen. Skipping.')
            return probes, {}
        neg_ref = neg_path
        logger.info(f'Auto-derived negative reference with {count} non-target sequences.')
    if not neg_ref:
        return probes, {}
    # Build negative DB and run BLAST
    neg_db = os.path.join(temp_dir, 'neg_blast_db')
    create_blast_db(neg_ref, neg_db)
    blast_output = _run_blast(probes, neg_db, temp_dir, getattr(args, 'blast_extra_args', []))
    if not blast_output:
        return probes, {}
    try:
        df = pd.read_csv(blast_output, sep='\t', names=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'])
    except Exception:
        return probes, {}
    # Filter for strong off-target hits
    strong_offtarget = df[(df['bitscore'] >= neg_bitscore) & (df['evalue'] <= neg_evalue)]
    report = {'off_target_hits': strong_offtarget}
    if strong_offtarget.empty:
        logger.info('Negative screen: no off-target hits found. All probes pass.')
        return probes, report
    # Reject probes with off-target hits
    reject_ids = set(strong_offtarget['qseqid'].unique())
    passed = [p for p in probes if p['pair_id'] not in reject_ids]
    rejected = len(probes) - len(passed)
    logger.info(f'Negative screen rejected {rejected} probes with off-target hits.')
    return passed, report