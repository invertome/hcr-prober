# hcr_prober/prober.py
import re
import numpy as np
import copy
import primer3
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
        if mask_sequences:
            mask_pattern = re.compile('|'.join(re.escape(m.upper()) for m in mask_sequences))
            all_windows = [w for w in all_windows if not mask_pattern.search(w['window_sequence'].upper())]
        audit['after_seq_mask'] = len(all_windows)
    thermo_passed = [w for w in all_windows if not su.has_homopolymer(w['window_sequence'], args.max_homopolymer) and args.min_gc <= tu.calculate_gc_content(w['window_sequence']) <= args.max_gc]
    audit['after_thermo_filter'] = len(thermo_passed)
    balanced_gc_passed = []
    for w in thermo_passed:
        p1, p2 = w['window_sequence'][0:args.probe_len], w['window_sequence'][args.probe_len + args.spacer_len:]
        gc1, gc2 = tu.calculate_gc_content(p1), tu.calculate_gc_content(p2)
        if abs(gc1 - gc2) <= args.max_gc_diff:
            w['probe_up_target'], w['probe_dn_target'] = p2, p1
            w['gc_dn'], w['gc_up'] = gc1, gc2
            balanced_gc_passed.append(w)
    audit['after_gc_balance_filter'] = len(balanced_gc_passed)
    tm_passed = []
    na = getattr(args, 'na_conc', 50)
    mg = getattr(args, 'mg_conc', 0)
    dntps = getattr(args, 'dntp_conc', 0)
    dnac = getattr(args, 'dna_conc', 25)
    for w in balanced_gc_passed:
        tm1 = tu.calculate_tm(w['probe_dn_target'], dnac1=dnac, dnac2=dnac, Na=na, Mg=mg, dNTPs=dntps)
        tm2 = tu.calculate_tm(w['probe_up_target'], dnac1=dnac, dnac2=dnac, Na=na, Mg=mg, dNTPs=dntps)
        if args.min_tm <= tm1 <= args.max_tm and args.min_tm <= tm2 <= args.max_tm:
            w['tm_dn'], w['tm_up'] = tm1, tm2
            tm_passed.append(w)
    audit['after_tm_filter'] = len(tm_passed)
    if hasattr(args, 'max_hairpin_dg'):
        tm_passed = filter_by_structure(tm_passed, args)
        audit['after_structure_filter'] = len(tm_passed)
    return tm_passed, audit

def select_spatially_diverse_probes(probes_to_filter, args):
    """Selects the maximum possible number of non-overlapping probes using dynamic programming with binary search."""
    if not probes_to_filter: return []
    from bisect import bisect_right
    sorted_probes = sorted(probes_to_filter, key=lambda p: p['start_pos_rev'] + args.window_size)
    num_probes = len(sorted_probes)
    # Build array of end positions for binary search
    end_positions = [p['start_pos_rev'] + args.window_size for p in sorted_probes]
    # dp[i] = max probes considering probes 0..i
    # choice[i] = 'take' or 'skip' for backtracking
    dp = [0] * num_probes
    backtrack = [-1] * num_probes
    choice = ['take'] * num_probes
    # Base case: taking probe 0 gives 1
    dp[0] = 1
    for i in range(1, num_probes):
        # Option 1: skip probe i — inherit dp[i-1]
        dp[i] = dp[i - 1]
        choice[i] = 'skip'
        backtrack[i] = i - 1
        # Option 2: take probe i — find last non-overlapping predecessor via bisect
        threshold = sorted_probes[i]['start_pos_rev'] - args.min_probe_distance
        j = bisect_right(end_positions, threshold) - 1
        take_val = (dp[j] if j >= 0 else 0) + 1
        if take_val > dp[i]:
            dp[i] = take_val
            choice[i] = 'take'
            backtrack[i] = j  # -1 if no predecessor
    # Backtrack from the end to reconstruct the solution
    selected_probes = []
    i = num_probes - 1
    while i >= 0:
        if choice[i] == 'take':
            selected_probes.append(sorted_probes[i])
            i = backtrack[i]
        else:
            i = i - 1
    selected_probes.reverse()
    logger.info(f"Optimal spacing filter selected {len(selected_probes)} probes from {len(probes_to_filter)} specific candidates.")
    return selected_probes

def filter_by_structure(candidates, args):
    """Filter probes by secondary structure: hairpins, homodimers, heterodimers."""
    max_hp = getattr(args, 'max_hairpin_dg', -3.0)
    max_homo = getattr(args, 'max_homodimer_dg', -5.0)
    max_hetero = getattr(args, 'max_heterodimer_dg', -5.0)
    passed = []
    for w in candidates:
        dn, up = w['probe_dn_target'], w['probe_up_target']
        hp_dn = primer3.calcHairpin(dn).dg / 1000.0
        hp_up = primer3.calcHairpin(up).dg / 1000.0
        if hp_dn < max_hp or hp_up < max_hp:
            continue
        homo_dn = primer3.calcHomodimer(dn).dg / 1000.0
        homo_up = primer3.calcHomodimer(up).dg / 1000.0
        if homo_dn < max_homo or homo_up < max_homo:
            continue
        hetero = primer3.calcHeterodimer(dn, up).dg / 1000.0
        if hetero < max_hetero:
            continue
        w['hairpin_dg_dn'] = round(hp_dn, 2)
        w['hairpin_dg_up'] = round(hp_up, 2)
        w['homodimer_dg_dn'] = round(homo_dn, 2)
        w['homodimer_dg_up'] = round(homo_up, 2)
        w['heterodimer_dg'] = round(hetero, 2)
        passed.append(w)
    return passed

def format_probes_for_blast(thermo_candidates, gene_name, sequence, args):
    formatted_probes = []
    for i, p_data in enumerate(thermo_candidates):
        probe = {
            'pair_id': f'{gene_name}_cand_{i+1}',
            'probe_up_target': p_data['probe_up_target'],
            'probe_dn_target': p_data['probe_dn_target'],
            'start_pos_on_sense': len(sequence) - p_data['start_pos_rev'] - args.window_size,
            'start_pos_rev': p_data['start_pos_rev'],
        }
        # Carry forward all cached thermo/structure values
        for key in ('gc_dn', 'gc_up', 'tm_dn', 'tm_up',
                     'hairpin_dg_dn', 'hairpin_dg_up', 'homodimer_dg_dn',
                     'homodimer_dg_up', 'heterodimer_dg'):
            if key in p_data:
                probe[key] = p_data[key]
        formatted_probes.append(probe)
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
