# hcr_prober/isoform_analyzer.py
import os
from loguru import logger
import numpy as np
from Bio import pairwise2

def group_sequences_by_prefix(sequences, delimiter='_'):
    groups = {}
    for seq_id, sequence in sequences.items():
        prefix = seq_id.split(delimiter)[0]
        if prefix not in groups: groups[prefix] = {}
        groups[prefix][seq_id] = sequence
    logger.info(f'Identified {len(groups)} gene group(s): {', '.join(groups.keys())}')
    return groups

def merge_intervals(intervals):
    if not intervals: return []
    intervals.sort(key=lambda x: x[0])
    merged = [intervals[0]]
    for current_start, current_end in intervals[1:]:
        last_start, last_end = merged[-1]
        if current_start <= last_end:
            merged[-1] = (last_start, max(last_end, current_end))
        else:
            merged.append((current_start, current_end))
    return merged

def invert_intervals(sequence_length, intervals_to_mask):
    if not intervals_to_mask: return [(0, sequence_length)]
    inverted, current_pos = [], 0
    sorted_intervals = merge_intervals(intervals_to_mask)
    for start, end in sorted_intervals:
        if start > current_pos: inverted.append((current_pos, start))
        current_pos = max(current_pos, end)
    if current_pos < sequence_length: inverted.append((current_pos, sequence_length))
    return inverted

def find_common_regions(isoform_group):
    if len(isoform_group) < 2:
        logger.warning(f'Gene group has only one sequence. Entire sequence treated as "common".')
        ref_id = list(isoform_group.keys())[0]
        return ref_id, [(0, len(isoform_group[ref_id]))]

    ref_id = max(isoform_group, key=lambda k: len(isoform_group[k]))
    ref_seq = isoform_group[ref_id]
    logger.info(f'Using "{ref_id}" (length: {len(ref_seq)}) as reference for finding common regions.')
    other_isoforms = {k: v for k, v in isoform_group.items() if k != ref_id}

    coverage_array = np.zeros(len(ref_seq), dtype=int)
    alignment_params = {'match': 2, 'mismatch': -1, 'open_gap': -5, 'extend_gap': -2}
    logger.info(f"Aligning {len(other_isoforms)} isoform(s) against reference with parameters: {alignment_params}")

    for query_id, query_seq in other_isoforms.items():
        alignments = pairwise2.align.localms(ref_seq, query_seq, alignment_params['match'], alignment_params['mismatch'], alignment_params['open_gap'], alignment_params['extend_gap'])
        if not alignments:
            logger.warning(f'No significant alignment found between {ref_id} and {query_id}. Cannot establish common regions.')
            return ref_id, []

        best_score = alignments[0].score
        high_scoring_alignments = [aln for aln in alignments if aln.score > best_score * 0.8]
        isoform_coverage = np.zeros(len(ref_seq), dtype=bool)
        for aln in high_scoring_alignments:
            ref_pos = aln.start
            for ref_char, query_char in zip(aln.seqA, aln.seqB):
                if ref_char != '-':
                    if query_char != '-': isoform_coverage[ref_pos] = True
                    ref_pos +=1
        coverage_array += isoform_coverage
    
    common_indices = np.where(coverage_array == len(other_isoforms))[0]
    if common_indices.size == 0:
        logger.info("No regions were found to be common across all isoforms.")
        return ref_id, []

    gaps = np.diff(common_indices) != 1
    split_points = np.where(gaps)[0] + 1
    interval_indices = np.split(common_indices, split_points)
    raw_intervals = [(indices[0], indices[-1] + 1) for indices in interval_indices if indices.size > 0]
    min_len = 52
    final_intervals = [iv for iv in raw_intervals if (iv[1] - iv[0]) >= min_len]
    merged_common = merge_intervals(final_intervals)
    total_common_len = sum(end - start for start, end in merged_common)
    logger.success(f'Found {len(merged_common)} common region(s) totaling {total_common_len} bp.')
    return ref_id, merged_common