# hcr_prober/main.py
import argparse, os, sys, json, shutil, itertools, tempfile
from loguru import logger
from . import file_io, prober, blast_wrapper, isoform_analyzer, swapper, visualization

def setup_logging():
    logger.remove(); logger.add(sys.stderr, format='<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <level>{message}</level>')

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
    thermo_group.add_argument('--min-gc', type=float, default=0.0, help='(Optional) Min %%GC. Default is off.')
    thermo_group.add_argument('--max-gc', type=float, default=100.0, help='(Optional) Max %%GC. Default is off.')
    thermo_group.add_argument('--min-tm', type=float, default=0.0, help='(Optional) Min Tm (Celsius). Default is off.')
    thermo_group.add_argument('--max-tm', type=float, default=100.0, help='(Optional) Max Tm (Celsius). Default is off.')
    thermo_group.add_argument('--max-homopolymer', type=int, default=4); thermo_group.add_argument('--max-gc-diff', type=float, default=15.0)
    blast_group.add_argument('--blast-ref', help='Path to FASTA for positive BLAST screen.')
    # NEW: Add argument for positive selection strategy
    blast_group.add_argument('--positive-selection-strategy', choices=['any-strong-hit', 'most-likely-transcript', 'specific-transcript'], default='any-strong-hit', help='Strategy for selecting probes based on BLAST hits.')
    # NEW: Add argument for the specific transcript ID, used with 'specific-transcript' strategy
    blast_group.add_argument('--target-transcript-id', help='The exact transcript ID to target when using the \'specific-transcript\' strategy.')
    blast_group.add_argument('--min-bitscore', type=float, default=75.0); blast_group.add_argument('--max-evalue', type=float, default=1e-10)
    blast_group.add_argument('--blast-extra-args', type=str, default='')
    adv_group.add_argument('--window-size', type=int, default=52); adv_group.add_argument('--probe-len', type=int, default=25); adv_group.add_argument('--spacer-len', type=int, default=2)

def run_design_pipeline(gene_name, seq, amplifier, output_dir, temp_dir, args):
    job_name = f'{gene_name}-{amplifier}'
    logger.info(f'--- Starting Job: {job_name} ---')
    args.job_name = job_name

    candidate_probes, audit_trail = prober.design_probes(seq, args)
    formatted_probes = prober.format_probes(candidate_probes, gene_name, amplifier, args.amplifiers, seq, args)
    audit_trail['after_formatting'] = len(formatted_probes)

    blast_passed_probes, blast_reports = blast_wrapper.filter_probes_by_blast(formatted_probes, args, temp_dir)
    audit_trail['after_blast_filter'] = len(blast_passed_probes)
    
    final_probes = prober.subsample_probes(blast_passed_probes, args.max_probes)
    audit_trail['after_subsampling'] = len(final_probes)
    
    file_io.write_outputs(final_probes, seq, gene_name, amplifier, args, blast_reports, audit_trail)

def main():
    setup_logging()
    config = file_io.load_config('hcr-prober.yaml')
    parser = argparse.ArgumentParser(description='HCR-prober v1.4.0', formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest='command', required=True)

    p_design = subparsers.add_parser('design', help='Design probes for standard transcripts.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_design.add_argument('-i', '--input', required=True, help='Input FASTA file of target transcript(s).')
    p_design.add_argument('-o', '--output-dir', default='hcr_prober_output')
    p_design.add_argument('--gene-name', help='Process only this transcript ID. If not set, all are processed.')
    p_design.add_argument('--pool-name', help='Custom name for the probe pool.')
    add_shared_design_args(p_design)

    p_iso = subparsers.add_parser('isoform-split', help='Design probes for common and unique regions of isoforms.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_iso.add_argument('-i', '--input', required=True, help='Input FASTA with ALL isoforms for ONE OR MORE genes.')
    p_iso.add_argument('-o', '--output-dir', default='hcr_isoform_output')
    p_iso.add_argument('--gene-prefix', nargs='+', required=True, help='Gene prefix(es) to process from the input file.')
    p_iso.add_argument('--delimiter', default='_', help='Delimiter separating gene prefix from isoform ID.')
    add_shared_design_args(p_iso)

    p_swap = subparsers.add_parser('swap', help='Swap amplifiers on existing probe file(s).', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_swap.add_argument('--input-probes', required=True, help='Path to a single .xlsx file or a directory containing them.')
    p_swap.add_argument('--output-dir', default='swapped_probes')
    p_swap.add_argument('--new-amplifier', required=True, help='ID of the new amplifier to swap to.')

    parser.set_defaults(**config)
    args = parser.parse_args()
    
    args.amplifiers = file_io.load_amplifiers(os.path.dirname(os.path.realpath(__file__)))
    if not args.amplifiers: sys.exit(1)

    if args.command in ['design', 'isoform-split']:
        check_dependencies()
        # Validate arguments
        if args.positive_selection_strategy == 'specific-transcript' and not args.target_transcript_id:
            logger.critical('FATAL: The \'specific-transcript\' strategy requires you to provide --target-transcript-id.')
            sys.exit(1)
        for amp in args.amplifier:
            if amp not in args.amplifiers: logger.critical(f'Amplifier \'{amp}\' not found.'); sys.exit(1)
        os.makedirs(args.output_dir, exist_ok=True)
        args.blast_db_positive = blast_wrapper.create_blast_db(args.blast_ref, args.db_path)
        args.blast_extra_args = args.blast_extra_args.split()
        
    if args.command == 'design':
        sequences = file_io.read_fasta(args.input)
        sequences_to_process = {k: v for k, v in sequences.items() if not args.gene_name or k == args.gene_name}
        if not sequences_to_process: logger.error(f'No sequences found in \'{args.input}\' that match --gene-name \'{args.gene_name}\'.'); sys.exit(1)
        job_list = list(itertools.product(sequences_to_process.items(), args.amplifier))
        logger.info(f'Starting {len(job_list)} total \'design\' job(s)...')
        with tempfile.TemporaryDirectory(prefix='hcr_prober_') as temp_dir:
            for (gene_name, seq), amplifier in job_list:
                run_design_pipeline(gene_name, seq, amplifier, args.output_dir, temp_dir, args)

    elif args.command == 'isoform-split':
        sequences = file_io.read_fasta(args.input)
        grouped_isoforms = isoform_analyzer.group_sequences_by_prefix(sequences, args.delimiter)
        with tempfile.TemporaryDirectory(prefix='hcr_prober_') as temp_dir:
            for prefix in args.gene_prefix:
                if prefix not in grouped_isoforms:
                    logger.warning(f'Gene prefix \'{prefix}\' not found. Skipping.'); continue
                logger.info(f'===== Analyzing Gene Group: {prefix} =====')
                iso_group = grouped_isoforms[prefix]
                ref_id, common_intervals = isoform_analyzer.find_common_regions(iso_group)
                ref_seq = iso_group[ref_id]

                unique_intervals = isoform_analyzer.invert_intervals(len(ref_seq), common_intervals)
                common_output_dir = os.path.join(args.output_dir, prefix, 'common_probes')
                args.mask_regions = ','.join([f'{s+1}-{e}' for s, e in unique_intervals]) if unique_intervals else None
                logger.info(f'--- Designing probes for COMMON regions of {prefix} ---')
                for amp in args.amplifier:
                    run_design_pipeline(f'{prefix}_common', ref_seq, amp, common_output_dir, temp_dir, args)
                
                for iso_id, iso_seq in iso_group.items():
                    unique_output_dir = os.path.join(args.output_dir, prefix, 'isoform_specific_probes')
                    args.mask_regions = ','.join([f'{s+1}-{e}' for s, e in common_intervals]) if common_intervals else None
                    logger.info(f'--- Designing probes for UNIQUE regions of {iso_id} ---')
                    for amp in args.amplifier:
                        run_design_pipeline(iso_id, iso_seq, amp, unique_output_dir, temp_dir, args)

    elif args.command == 'swap':
        swapper.swap_amplifiers(args, args.amplifiers)

    logger.success('HCR-prober pipeline finished.')

if __name__ == '__main__':
    main()