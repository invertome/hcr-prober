# hcr_prober/prober.py
import numpy as np
import copy
from loguru import logger
from .utils import sequence_utils as su, thermo_utils as tu
from . import file_io

def generate_thermo_candidates(sequence, args):
    """Generates all possible probe windows and filters them based on intrinsic thermodynamic properties."""
    rev_comp_seq = su.reverse_complement(sequence)
    audit = {}
    all_windows = [{'window_sequence': rev_comp_seq[i:i+args.window_size], 'start_pos_rev': i} for i in range(len(rev_comp_seq) - args.window_size + 1)]
    audit['initial_windows'] = len(all_windows)
    skip_rev_coord = len(rev_comp_seq) - args.skip_5prime
    all_windows = [w for w in all_windows if w['start_pos_rev'] < skip_rev_coord]
    audit['after_5prime_skip'] = len(all_windows)
    if getattr(args, 'mask_regions', None):
        mask_intervals = su.parse_mask_regions(args.mask_regions)
        if mask_intervals:
            original_len, unmasked_windows = len(sequence), []
            for w in all_windows:
                start_on_sense = original_len - w['start_pos_rev'] - args.window_size
                end_on_sense = start_on_sense + args.window_size
                is_masked = any(max(start_on_sense, mask_start) < min(end_on_sense, mask_end) for mask_start, mask_end in mask_intervals)
                if not is_masked: unmasked_windows.append(w)
            all_windows = unmasked_windows
            audit['after_region_mask'] = len(all_windows)
    if getattr(args, 'mask_sequences', None):
        mask_sequences = list(file_io.read_fasta(args.mask_sequences).values())
        all_windows = [w for w in all_windows if not any(m.upper() in w['window_sequence'].upper() for m in mask_sequences)]
        audit['after_seq_mask'] = len(all_windows)
    thermo_passed = [w for w in all_windows if not su.has_homopolymer(w['window_sequence'], args.max_homopolymer) and args.min_gc <= tu.calculate_gc_content(w['window_sequence']) <= args.max_gc]
    audit['after_thermo_filter'] = len(thermo_passed)
    balanced_gc_passed = []
    for w in thermo_passed:
        p1, p2 = w['window_sequence'][0:args.probe_len], w['window_sequence'][args.probe_len + args.spacer_len:]
        if abs(tu.calculate_gc_content(p1) - tu.calculate_gc_content(p2)) <= args.max_gc_diff:
            w['probe_up_target'], w['probe_dn_target'] = p2, p1; balanced_gc_passed.append(w)
    audit['after_gc_balance_filter'] = len(balanced_gc_passed)
    tm_passed = []
    for w in balanced_gc_passed:
        tm1, tm2 = tu.calculate_tm(w['probe_dn_target']), tu.calculate_tm(w['probe_up_target'])
        if args.min_tm <= tm1 <= args.max_tm and args.min_tm <= tm2 <= args.max_tm:
            w['tm_dn'], w['tm_up'] = tm1, tm2; tm_passed.append(w)
    audit['after_tm_filter'] = len(tm_passed)
    return tm_passed, audit

def select_spatially_diverse_probes(probes_to_filter, args):
    """Selects the maximum possible number of non-overlapping probes using a dynamic programming algorithm."""
    if not probes_to_filter: return []
    sorted_probes = sorted(probes_to_filter, key=lambda p: p['start_pos_rev'] + args.window_size)
    num_probes = len(sorted_probes)
    dp, backtrack = [1] * num_probes, [-1] * num_probes
    required_footprint = args.window_size + args.min_probe_distance
    for i in range(num_probes):
        for j in range(i):
            if sorted_probes[j]['start_pos_rev'] + required_footprint <= sorted_probes[i]['start_pos_rev']:
                if dp[j] + 1 > dp[i]:
                    dp[i], backtrack[i] = dp[j] + 1, j
    if not dp: return []
    max_len, best_end_index = 0, -1
    for i in range(num_probes):
        if dp[i] > max_len:
            max_len, best_end_index = dp[i], i
    selected_probes, curr_index = [], best_end_index
    while curr_index != -1:
        selected_probes.append(sorted_probes[curr_index])
        curr_index = backtrack[curr_index]
    selected_probes.reverse()
    logger.info(f"Optimal spacing filter selected {len(selected_probes)} probes from {len(probes_to_filter)} specific candidates.")
    return selected_probes

def format_probes_for_blast(thermo_candidates, gene_name, sequence, args):
    formatted_probes = []
    for i, p_data in enumerate(thermo_candidates):
        formatted_probes.append({
            'pair_id': f'{gene_name}_cand_{i+1}', 'probe_up_target': p_data['probe_up_target'],
            'probe_dn_target': p_data['probe_dn_target'],
            'start_pos_on_sense': len(sequence) - p_data['start_pos_rev'] - args.window_size,
            'start_pos_rev': p_data['start_pos_rev']
        })
    return formatted_probes

def finalize_probes(blueprint_probes, amplifier_name, amplifiers, args):
    final_probes = []
    amp_data = amplifiers.get(amplifier_name)
    if not amp_data: return []
    up_init, dn_init = amp_data['up'], amp_data['dn']
    up_spc, dn_spc = su.resolve_iupac_spacer(amp_data.get('upspc', '')), su.resolve_iupac_spacer(amp_data.get('dnspc', ''))
    gene_name = blueprint_probes[0]['pair_id'].split('_cand_')[0] if blueprint_probes else 'gene'
    sorted_blueprint = sorted(blueprint_probes, key=lambda p: p['start_pos_on_sense'])
    for i, probe in enumerate(sorted_blueprint):
        final_probe = copy.deepcopy(probe)
        # **FIX v1.9.6**: Preserve the original candidate ID for traceability in reports.
        final_probe['cand_id'] = probe['pair_id']
        final_probe['pair_num'] = i + 1
        final_probe['pair_id'] = f'{gene_name}_pair_{i+1}'
        final_probe['probe_up_final'] = f'{up_init}{up_spc}{probe["probe_up_target"]}'
        final_probe['probe_dn_final'] = f'{probe["probe_dn_target"]}{dn_spc}{dn_init}'
        final_probes.append(final_probe)
    return final_probes

def subsample_probes(probes, num_to_keep):
    if num_to_keep >= len(probes): return probes
    logger.info(f'Subsampling from {len(probes)} pairs to a max of {num_to_keep} for even coverage.')
    indices = np.linspace(0, len(probes) - 1, num=num_to_keep, dtype=int)
    return [probes[i] for i in indices]
