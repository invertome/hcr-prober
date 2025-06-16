# hcr_prober/isoform_analyzer.py
import os, tempfile, subprocess
from loguru import logger

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
        if current_start <= last_end: merged[-1] = (last_start, max(last_end, current_end))
        else: merged.append((current_start, current_end))
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
        logger.warning(f'Gene group has only one sequence. Entire sequence treated as \'common\'.')
        ref_id = list(isoform_group.keys())[0]
        return ref_id, [(0, len(isoform_group[ref_id]))]

    ref_id = max(isoform_group, key=lambda k: len(isoform_group[k]))
    ref_seq = isoform_group[ref_id]
    logger.info(f'Using \'{ref_id}\' (length: {len(ref_seq)}) as reference for finding common regions.')
    other_isoforms = {k: v for k, v in isoform_group.items() if k != ref_id}
    
    with tempfile.TemporaryDirectory(prefix='isoform_blast_') as temp_dir:
        ref_fasta_path = os.path.join(temp_dir, 'ref.fasta')
        db_name = os.path.join(temp_dir, 'ref_db')
        with open(ref_fasta_path, 'w') as f: f.write(f'>{ref_id}\n{ref_seq}\n')
        subprocess.run(['makeblastdb', '-in', ref_fasta_path, '-dbtype', 'nucl', '-out', db_name], check=True, capture_output=True)
        
        all_common_intervals = []
        for query_id, query_seq in other_isoforms.items():
            query_fasta_path = os.path.join(temp_dir, 'query.fasta')
            blast_out_path = os.path.join(temp_dir, 'blast.tsv')
            with open(query_fasta_path, 'w') as f: f.write(f'>{query_id}\n{query_seq}\n')
            cmd = ['blastn', '-query', query_fasta_path, '-db', db_name, '-out', blast_out_path, '-outfmt', '6 sstart send', '-perc_identity', '98', '-qcov_hsp_perc', '90']
            subprocess.run(cmd, check=True, capture_output=True)
            
            isoform_intervals = []
            if os.path.exists(blast_out_path) and os.path.getsize(blast_out_path) > 0:
                with open(blast_out_path, 'r') as f_blast:
                    for line in f_blast:
                        start, end = sorted([int(p) for p in line.strip().split()])
                        isoform_intervals.append((start - 1, end))
            all_common_intervals.append(set(isoform_intervals))

    if not all_common_intervals: return ref_id, []
    intersected_intervals = list(set.intersection(*all_common_intervals))
    merged_common = merge_intervals(intersected_intervals)
    total_common_len = sum(end - start for start, end in merged_common)
    logger.info(f'Found {len(merged_common)} common region(s) totaling {total_common_len} bp.')
    return ref_id, merged_common