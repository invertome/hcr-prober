# hcr_prober/prober.py
import numpy as np
from loguru import logger
from .utils import sequence_utils as su, thermo_utils as tu
from . import file_io

def design_probes(sequence, args):
    rev_comp_seq = su.reverse_complement(sequence)
    audit = {}

    all_windows = [{'window_sequence': rev_comp_seq[i:i+args.window_size], 'start_pos_rev': i} for i in range(len(rev_comp_seq) - args.window_size + 1)]
    audit['initial_windows'] = len(all_windows)
    
    skip_rev_coord = len(rev_comp_seq) - args.skip_5prime
    all_windows = [w for w in all_windows if w['start_pos_rev'] < skip_rev_coord]
    audit['after_5prime_skip'] = len(all_windows)
    
    if getattr(args, 'mask_sequences', None):
        mask_sequences = list(file_io.read_fasta(args.mask_sequences).values())
        all_windows = [w for w in all_windows if not any(m.upper() in w['window_sequence'].upper() for m in mask_sequences)]
        audit['after_seq_mask'] = len(all_windows)

    thermo_passed = [w for w in all_windows if not su.has_homopolymer(w['window_sequence'], args.max_homopolymer) and args.min_gc <= tu.calculate_gc_content(w['window_sequence']) <= args.max_gc]
    audit['after_thermo_filter'] = len(thermo_passed)
    
    balanced_gc_passed = []
    for w in thermo_passed:
        p1 = w['window_sequence'][0:args.probe_len]; p2 = w['window_sequence'][args.probe_len + args.spacer_len:]
        if abs(tu.calculate_gc_content(p1) - tu.calculate_gc_content(p2)) <= args.max_gc_diff:
            w['probe_up_target'] = p2; w['probe_dn_target'] = p1; balanced_gc_passed.append(w)
    audit['after_gc_balance_filter'] = len(balanced_gc_passed)
    
    tm_passed = []
    for w in balanced_gc_passed:
        tm1, tm2 = tu.calculate_tm(w['probe_dn_target']), tu.calculate_tm(w['probe_up_target'])
        if args.min_tm <= tm1 <= args.max_tm and args.min_tm <= tm2 <= args.max_tm:
            w['tm_dn'] = tm1; w['tm_up'] = tm2; tm_passed.append(w)
    audit['after_tm_filter'] = len(tm_passed)

    probe_pairs_filtered, last_pos = [], -1000
    for cand in tm_passed:
        if cand['start_pos_rev'] >= last_pos + args.window_size:
            probe_pairs_filtered.append(cand); last_pos = cand['start_pos_rev']
    audit['after_overlap_filter'] = len(probe_pairs_filtered)
    
    return probe_pairs_filtered, audit

def format_probes(probes_to_format, gene_name, amplifier_name, amplifiers, sequence, args):
    if not probes_to_format: return []
    try: amp_data = amplifiers[amplifier_name]
    except KeyError: logger.error(f'FATAL: Amplifier \'{amplifier_name}\' not found.'); return []
    up_init, dn_init = amp_data['up'], amp_data['dn']
    up_spc, dn_spc = su.resolve_iupac_spacer(amp_data['upspc']), su.resolve_iupac_spacer(amp_data['dnspc'])

    final_probes_list = []
    for i, pair_data in enumerate(probes_to_format):
        final_probes_list.append({
            'pair_id': f'{gene_name}_pair_{i+1}', 'pair_num': i + 1,
            'probe_up_final': f'{up_init}{up_spc}{pair_data['probe_up_target']}',
            'probe_dn_final': f'{pair_data['probe_dn_target']}{dn_spc}{dn_init}',
            'probe_up_target': pair_data['probe_up_target'], 'probe_dn_target': pair_data['probe_dn_target'],
            'start_pos_on_sense': len(sequence) - pair_data['start_pos_rev'] - args.window_size
        })
    return final_probes_list

def subsample_probes(probes, num_to_keep):
    if num_to_keep >= len(probes): return probes
    logger.info(f'Subsampling from {len(probes)} pairs to a max of {num_to_keep} for even coverage.')
    indices = np.linspace(0, len(probes) - 1, num=num_to_keep, dtype=int)
    return [probes[i] for i in indices]