# hcr_prober/file_io.py
import os, json, pandas as pd, yaml, glob
import numpy as np
from . import visualization
from loguru import logger
from Bio import SeqIO

def read_fasta(file_path):
    if not file_path: return {}
    try: return {rec.id: str(rec.seq) for rec in SeqIO.parse(file_path, 'fasta')}
    except FileNotFoundError: logger.critical(f'FASTA not found: {file_path}'); exit(1)

def load_amplifiers(pkg_path):
    amplifiers, plugin_dir = {}, os.path.join(pkg_path, 'config', 'amplifiers')
    if not os.path.isdir(plugin_dir):
        logger.error(f'Amplifier dir not found: {plugin_dir}')
        return {}
    for pf in glob.glob(os.path.join(plugin_dir, '*.json')):
        try:
            # **FIX v1.9.5**: Corrected syntax for try/with block
            with open(pf, 'r') as f:
                amplifiers.update(json.load(f))
        except Exception as e:
            logger.warning(f'Could not load amplifier plugin \'{os.path.basename(pf)}\': {e}')
    if amplifiers:
        logger.info(f'Loaded {len(amplifiers)} amplifiers: {', '.join(amplifiers.keys())}')
    else:
        logger.error('No amplifiers were loaded.')
    return amplifiers

def load_config(config_path):
    if os.path.exists(config_path):
        with open(config_path, 'r') as f: return yaml.safe_load(f) or {}
    return {}

def write_outputs(probes, sequence, gene_name, amplifier, args, blast_reports, audit_trail):
    amp_dir = os.path.join(args.output_dir, gene_name, amplifier)
    os.makedirs(amp_dir, exist_ok=True)
    with open(os.path.join(amp_dir, f'{gene_name}_{amplifier}_summary.txt'), 'w') as f:
        f.write(f'HCR-prober v1.9.5 Summary for: {gene_name} | Amplifier: {amplifier}\n{'='*70}\n')
        if not probes: f.write('\n*** PIPELINE FAILED TO PRODUCE ANY FINAL PROBES ***\n')
        f.write('\n--- Run Parameters ---\n')
        f.write(f'  Target Sequence Length: {len(sequence)} nt\n')
        f.write(f'  Min Probe Distance: {args.min_probe_distance} nt\n')
        f.write(f'  GC Range: {args.min_gc}-{args.max_gc} %\n')
        f.write(f'  Tm Range: {args.min_tm}-{args.max_tm} C\n')
        if args.blast_ref:
            f.write(f'  BLAST Reference: {os.path.basename(args.blast_ref)}\n')
            f.write(f'  Discovery Bitscore Cutoff: {args.min_bitscore}\n')
            f.write(f'  Discovery E-value Cutoff: {args.max_evalue}\n')
            f.write(f'  Positive Selection Strategy: {args.positive_selection_strategy}\n')
            if args.positive_selection_strategy == 'specific-id': f.write(f'    Target Transcript ID: {args.target_transcript_id}\n')
        f.write('\n\n--- Filtering Funnel (Blueprint Architecture) ---\n')
        audit_key_map = {
            'initial_windows': 'Initial Windows', 'after_5prime_skip': 'After 5\' Skip', 'after_region_mask': 'After Region Masking', 'after_seq_mask': 'After Sequence Masking',
            'after_thermo_filter': 'Thermo-stable Candidates', 'after_gc_balance_filter': 'After GC Balance', 'after_tm_filter': 'After Tm',
            'after_blast_filter': 'Specific Candidates (Post-BLAST)', 'after_spacing_filter': 'Spatially-Diverse Blueprint Probes',
            'after_subsampling': 'After Final Subsampling'
        }
        ordered_keys = ['initial_windows', 'after_5prime_skip', 'after_region_mask', 'after_seq_mask', 'after_thermo_filter', 'after_gc_balance_filter', 'after_tm_filter', 'after_blast_filter', 'after_spacing_filter', 'after_subsampling']
        for key in ordered_keys:
            if key in audit_trail: f.write(f'  {audit_trail[key]:>7} {audit_key_map.get(key, key)}\n')
        f.write(f'  ---------------------------\n  {len(probes):>7} Final Probe Pairs\n')
        _write_blast_report_section(f, blast_reports, probes)
    if probes:
        logger.success(f'Writing {len(probes)} probe pairs for \'{gene_name} - {amplifier}\' to {amp_dir}/')
        pool_name = getattr(args, 'pool_name', None) or f'{amplifier}_{gene_name}_PP{len(probes)}'
        order_data = {'Pool name': [], 'Sequence': []}
        for p in sorted(probes, key=lambda x: x['pair_num']):
            order_data['Pool name'].extend([pool_name, pool_name]); order_data['Sequence'].extend([p['probe_dn_final'], p['probe_up_final']])
        pd.DataFrame(order_data).to_excel(os.path.join(amp_dir, f'{gene_name}_{amplifier}_order.xlsx'), index=False)
        with open(os.path.join(amp_dir, f'{gene_name}_{amplifier}_probes.fasta'), 'w') as f:
            for p in sorted(probes, key=lambda x: x['pair_num']): f.write(f'>{p['pair_id']}_A\n{p['probe_dn_final']}\n>{p['pair_id']}_B\n{p['probe_up_final']}\n')
        visualization.generate_svg_probe_map(probes, len(sequence), amplifier, gene_name, os.path.join(amp_dir, f'{gene_name}_{amplifier}_probe_map.svg'))
    else: logger.warning(f'No final probes for {gene_name} with amplifier {amplifier}. Report created.')

def _write_blast_report_section(f, blast_reports, final_probes):
    """Writes the detailed BLAST report, highlighting probes that were in the final set."""
    if not blast_reports: return

    # **FIX v1.9.6**: Use the preserved 'cand_id' to correctly identify final probes.
    final_probe_cand_ids = {p.get('cand_id') for p in final_probes}

    f.write('\n' + '='*70 + '\n--- DETAILED BLAST REPORT ---\n' + '='*70 + '\n')
    for screen_type, report in blast_reports.items():
        f.write(f'\n--- {screen_type.upper()} SCREEN ---\n')
        pd.set_option('display.max_rows', None); pd.set_option('display.width', 1000); pd.set_option('display.max_colwidth', None)

        if 'strong_hits' in report and not report['strong_hits'].empty:
            f.write('\n[+] Plausible BLAST Hits. Probes in the final set are marked with a *.\n\n')

            df_hits = report['strong_hits'].copy()

            # Use a blank space for non-selected probes for better alignment
            df_hits['SELECTED'] = np.where(df_hits['pair_id'].isin(final_probe_cand_ids), '*', ' ')

            cols = ['SELECTED'] + [col for col in df_hits.columns if col != 'SELECTED']
            df_hits = df_hits[cols]

            # **FIX v1.9.6**: Use to_string(index=False) to fix table formatting.
            f.write(df_hits.to_string(index=False) + '\n')
        else:
            f.write('\n[+] No BLAST hits passed discovery thresholds.\n')
