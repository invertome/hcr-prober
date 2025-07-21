# hcr_prober/main.py
import argparse, os, sys, shutil, tempfile, copy
from loguru import logger
from . import file_io, prober, blast_wrapper, isoform_analyzer, swapper

def setup_logging(): logger.remove(); logger.add(sys.stderr, format='<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <level>{message}</level>')

def check_dependencies():
    if not shutil.which('blastn') or not shutil.which('makeblastdb'):
        logger.critical('FATAL ERROR: NCBI BLAST+ is not installed or not in your system\'s PATH.'); sys.exit(1)
    logger.info('Dependency check passed: NCBI BLAST+ found.')

def add_shared_design_args(parser):
    proc_group = parser.add_argument_group('Processing & Performance Arguments')
    design_group = parser.add_argument_group('Core Design Parameters')
    thermo_group = parser.add_argument_group('Thermodynamic & Sequence Filters')
    blast_group = parser.add_argument_group('BLAST Specificity Filters')
    adv_group = parser.add_argument_group('Advanced Structural Parameters')
    proc_group.add_argument('--force', action='store_true', help='Force re-run and ignore cached results.')
    proc_group.add_argument('--db-path', help='Permanent directory to store/find BLAST databases.')
    design_group.add_argument('--amplifier', nargs='+', required=True, help='One or more HCR amplifier IDs.')
    design_group.add_argument('--max-probes', type=int, default=33, help='Maximum number of probe pairs to select.')
    design_group.add_argument('--skip-5prime', type=int, default=100, help='Nucleotides to exclude at the 5\' end.')
    design_group.add_argument('--mask-sequences', help='Path to a FASTA of sequences to exclude (e.g., repeats).')
    design_group.add_argument('--mask-regions', help='Comma-separated nt regions to exclude (e.g., \'1-100,500-650\').')
    design_group.add_argument('--min-probe-distance', type=int, default=2, help='Minimum distance (in nt) between probe footprints.')
    thermo_group.add_argument('--min-gc', type=float, default=40.0); thermo_group.add_argument('--max-gc', type=float, default=60.0)
    thermo_group.add_argument('--min-tm', type=float, default=40.0); thermo_group.add_argument('--max-tm', type=float, default=55.0)
    thermo_group.add_argument('--max-homopolymer', type=int, default=4); thermo_group.add_argument('--max-gc-diff', type=float, default=15.0)
    blast_group.add_argument('--blast-ref', help='Path to FASTA for positive BLAST screen.')
    blast_group.add_argument('--positive-selection-strategy', choices=['any-strong-hit', 'best-coverage', 'specific-id'], default='any-strong-hit', help='Global BLAST strategy for the design command.')
    blast_group.add_argument('--target-transcript-id', help='The exact transcript ID to target for the \'specific-id\' strategy.')
    blast_group.add_argument('--min-bitscore', type=float, default=75.0); blast_group.add_argument('--max-evalue', type=float, default=1e-10)
    blast_group.add_argument('--blast-extra-args', type=str, default='')
    adv_group.add_argument('--window-size', type=int, default=52); adv_group.add_argument('--probe-len', type=int, default=25); adv_group.add_argument('--spacer-len', type=int, default=2)

def create_probe_blueprint(gene_name, seq, temp_dir, args):
    """Performs all amplifier-independent steps (thermo, blast, spacing) once per gene."""
    thermo_candidates, audit_trail = prober.generate_thermo_candidates(seq, args)
    if not thermo_candidates: return None, None, audit_trail
    blast_formatted_probes = prober.format_probes_for_blast(thermo_candidates, gene_name, seq, args)
    specific_probes, blast_reports = blast_wrapper.filter_probes_by_blast(blast_formatted_probes, args, temp_dir)
    audit_trail['after_blast_filter'] = len(specific_probes)
    if not specific_probes: return None, blast_reports, audit_trail
    spaced_probes = prober.select_spatially_diverse_probes(specific_probes, args)
    audit_trail['after_spacing_filter'] = len(spaced_probes)
    return spaced_probes, blast_reports, audit_trail

def main():
    setup_logging()
    config = file_io.load_config('hcr-prober.yaml')
    parser = argparse.ArgumentParser(description='HCR-prober v1.9.5', formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest='command', required=True)
    p_design = subparsers.add_parser('design', help='Design probes for standard transcripts.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_design.add_argument('-i', '--input', required=True, help='Input FASTA file.')
    p_design.add_argument('-o', '--output-dir', default='hcr_prober_output')
    p_design.add_argument('--gene-name', help='Process only this transcript ID from input.')
    p_design.add_argument('--pool-name', help='Custom name for the probe pool (replaces default).')
    add_shared_design_args(p_design)

    p_iso = subparsers.add_parser('isoform-split', help='Design probes for common/unique regions of isoforms.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_iso.add_argument('-i', '--input', required=True, help='Input FASTA with ALL isoforms for ONE OR MORE genes.')
    p_iso.add_argument('-o', '--output-dir', default='hcr_isoform_output')
    p_iso.add_argument('--gene-prefix', nargs='+', required=True, help='Gene prefix(es) to process.')
    p_iso.add_argument('--delimiter', default='_', help='Delimiter for gene prefix.')
    iso_blast_group = p_iso.add_argument_group('Isoform-Split BLAST Strategy (v1.9.4+)')
    iso_blast_group.add_argument('--common-strategy', choices=['any-strong-hit', 'best-coverage', 'specific-id'], default='any-strong-hit', help='BLAST strategy for COMMON probes.')
    iso_blast_group.add_argument('--unique-strategy', choices=['any-strong-hit', 'best-coverage', 'specific-id'], default='best-coverage', help='BLAST strategy for UNIQUE probes.')
    add_shared_design_args(p_iso)

    p_swap = subparsers.add_parser('swap', help='Swap amplifiers on existing probe file(s).', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_swap.add_argument('--input-probes', required=True, help='Path to a single .xlsx file or a directory of them.')
    p_swap.add_argument('--output-dir', default='swapped_probes')
    p_swap.add_argument('--new-amplifier', required=True, help='ID of the new amplifier.')

    parser.set_defaults(**config)
    args = parser.parse_args()
    args.amplifiers = file_io.load_amplifiers(os.path.dirname(os.path.realpath(__file__)))
    if not args.amplifiers: sys.exit(1)

    if args.command in ['design', 'isoform-split']:
        check_dependencies()
        strat = getattr(args, 'positive_selection_strategy', None)
        common_strat = getattr(args, 'common_strategy', None)
        if (strat == 'specific-id' or common_strat == 'specific-id') and not args.target_transcript_id: logger.critical('FATAL: \'specific-id\' strategy requires --target-transcript-id.'); sys.exit(1)
        for amp in args.amplifier:
            if amp not in args.amplifiers: logger.critical(f'Amplifier \'{amp}\' not found.'); sys.exit(1)
        os.makedirs(args.output_dir, exist_ok=True)
        args.blast_db_positive = blast_wrapper.create_blast_db(args.blast_ref, args.db_path)
        args.blast_extra_args = args.blast_extra_args.split()

    if args.command == 'design':
        sequences = file_io.read_fasta(args.input)
        sequences_to_process = {k: v for k, v in sequences.items() if not args.gene_name or k == args.gene_name}
        if not sequences_to_process: logger.error(f'No sequences in \'{args.input}\' match --gene-name \'{args.gene_name}\'.'); sys.exit(1)
        logger.info(f'Starting design jobs for {len(sequences_to_process)} gene(s) and {len(args.amplifier)} amplifier(s).')
        with tempfile.TemporaryDirectory(prefix='hcr_prober_') as temp_dir:
            for gene_name, seq in sequences_to_process.items():
                logger.info(f'--- Creating probe blueprint for: {gene_name} ---')
                args.job_name = gene_name
                blueprint, blast_reports, audit_trail = create_probe_blueprint(gene_name, seq, temp_dir, args)
                if not blueprint:
                    logger.warning(f'Could not create a probe blueprint for {gene_name}. Writing empty report(s).')
                    for amp in args.amplifier: file_io.write_outputs([], seq, gene_name, amp, args, blast_reports or {}, audit_trail)
                    continue
                logger.success(f'Successfully created blueprint with {len(blueprint)} probes for {gene_name}.')
                for amp in args.amplifier:
                    finalized = prober.finalize_probes(blueprint, amp, args.amplifiers, args)
                    subsampled = prober.subsample_probes(finalized, args.max_probes)
                    amp_audit = {**audit_trail, 'after_subsampling': len(subsampled)}
                    file_io.write_outputs(subsampled, seq, gene_name, amp, args, blast_reports, amp_audit)

    elif args.command == 'isoform-split':
        sequences = file_io.read_fasta(args.input)
        grouped_isoforms = isoform_analyzer.group_sequences_by_prefix(sequences, args.delimiter)
        with tempfile.TemporaryDirectory(prefix='hcr_prober_') as temp_dir:
            for prefix in args.gene_prefix:
                if prefix not in grouped_isoforms: logger.warning(f'Gene prefix \'{prefix}\' not in input. Skipping.'); continue
                logger.info(f'========== Analyzing Gene Group: {prefix} ==========')
                iso_group = grouped_isoforms[prefix]
                ref_id, common_intervals = isoform_analyzer.find_common_regions(iso_group)
                ref_seq, unique_intervals = iso_group[ref_id], isoform_analyzer.invert_intervals(len(iso_group[ref_id]), common_intervals)
                common_args = copy.deepcopy(args)
                common_args.positive_selection_strategy = args.common_strategy
                logger.info(f'--- Creating COMMON probe blueprint for {prefix} (Strategy: {common_args.positive_selection_strategy}) ---')
                common_args.output_dir = os.path.join(args.output_dir, prefix, 'common_probes')
                common_args.mask_regions = ','.join([f'{s+1}-{e}' for s, e in unique_intervals]) if unique_intervals else None
                common_args.job_name = f'{prefix}_common'
                blueprint, blast_reports, audit_trail = create_probe_blueprint(f'{prefix}_common', ref_seq, temp_dir, common_args)
                if blueprint: logger.success(f'Created COMMON blueprint with {len(blueprint)} probes for {prefix}.')
                for amp in args.amplifier:
                    finalized = prober.finalize_probes(blueprint or [], amp, args.amplifiers, common_args)
                    subsampled = prober.subsample_probes(finalized, args.max_probes)
                    amp_audit = {**(audit_trail or {}), 'after_subsampling': len(subsampled)}
                    file_io.write_outputs(subsampled, ref_seq, f'{prefix}_common', amp, common_args, blast_reports or {}, amp_audit)
                unique_output_base, unique_mask = os.path.join(args.output_dir, prefix, 'isoform_specific_probes'), ','.join([f'{s+1}-{e}' for s, e in common_intervals]) if common_intervals else None
                for iso_id, iso_seq in iso_group.items():
                    unique_args = copy.deepcopy(args)
                    unique_args.positive_selection_strategy = args.unique_strategy
                    logger.info(f'--- Creating UNIQUE probe blueprint for {iso_id} (Strategy: {unique_args.positive_selection_strategy}) ---')
                    unique_args.output_dir, unique_args.mask_regions, unique_args.job_name = unique_output_base, unique_mask, iso_id
                    blueprint, blast_reports, audit_trail = create_probe_blueprint(iso_id, iso_seq, temp_dir, unique_args)
                    if blueprint: logger.success(f'Created UNIQUE blueprint with {len(blueprint)} probes for {iso_id}.')
                    for amp in args.amplifier:
                        finalized = prober.finalize_probes(blueprint or [], amp, args.amplifiers, unique_args)
                        subsampled = prober.subsample_probes(finalized, args.max_probes)
                        amp_audit = {**(audit_trail or {}), 'after_subsampling': len(subsampled)}
                        file_io.write_outputs(subsampled, iso_seq, iso_id, amp, unique_args, blast_reports or {}, amp_audit)

    elif args.command == 'swap': swapper.swap_amplifiers(args, args.amplifiers)
    logger.success('HCR-prober pipeline finished.')

if __name__ == '__main__': main()