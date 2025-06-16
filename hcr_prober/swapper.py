# hcr_prober/swapper.py
import pandas as pd, os, glob
from loguru import logger
from .utils import sequence_utils as su

def _swap_amplifiers_on_file(input_path, output_path, args, amplifiers):
    try:
        df = pd.read_excel(input_path, engine='openpyxl')
        if 'Sequence' not in df.columns or 'Pool name' not in df.columns:
            logger.warning(f'Skipping \'{input_path}\': Missing required columns.')
            return False
    except Exception as e:
        logger.error(f'Failed to read or parse \'{input_path}\': {e}'); return False

    new_amp_data = amplifiers[args.new_amplifier]
    new_up, new_dn = new_amp_data['up'], new_amp_data['dn']
    new_upspc = su.resolve_iupac_spacer(new_amp_data['upspc'])
    new_dnspc = su.resolve_iupac_spacer(new_amp_data['dnspc'])
    
    initiator_map = {data['up']: ('up', data.get('upspc', '')) for data in amplifiers.values()}
    initiator_map.update({data['dn']: ('dn', data.get('dnspc', '')) for data in amplifiers.values()})
    
    new_sequences = []
    for seq in df['Sequence']:
        seq_upper, found = seq.upper(), False
        for old_initiator, (initiator_type, old_spacer_iupac) in initiator_map.items():
            spacer_len = len(old_spacer_iupac)
            if seq_upper.startswith(old_initiator.upper()):
                target_seq = seq[len(old_initiator) + spacer_len:]
                new_sequences.append(f'{new_up}{new_upspc}{target_seq}'); found = True; break
            elif seq_upper.endswith(old_initiator.upper()):
                target_seq = seq[:-(len(old_initiator) + spacer_len)]
                new_sequences.append(f'{target_seq}{new_dnspc}{new_dn}'); found = True; break
        if not found:
            logger.warning(f'Could not identify initiator for \'{seq[:10]}...\'. Keeping as-is.')
            new_sequences.append(seq)

    df['Sequence'] = new_sequences
    old_pool_name = df['Pool name'][0]
    try:
        old_amp = old_pool_name.split('_')[0]
        df['Pool name'] = old_pool_name.replace(old_amp, args.new_amplifier) if old_amp in amplifiers else f'{args.new_amplifier}_{old_pool_name}'
    except (IndexError, AttributeError):
        df['Pool name'] = f'{args.new_amplifier}_{old_pool_name}'

    df.to_excel(output_path, index=False)
    logger.success(f'Saved swapped file to: {output_path}')
    return True

def swap_amplifiers(args, amplifiers):
    if args.new_amplifier not in amplifiers:
        logger.error(f'New amplifier \'{args.new_amplifier}\' not found.'); return
    os.makedirs(args.output_dir, exist_ok=True)
    if os.path.isdir(args.input_probes):
        logger.info(f'Batch mode: Processing .xlsx files in \'{args.input_probes}\'')
        files_to_process = glob.glob(os.path.join(args.input_probes, '**', '*.xlsx'), recursive=True)
        if not files_to_process: logger.warning(f'No .xlsx files found in \'{args.input_probes}\'.'); return
        for file_path in files_to_process:
            logger.info(f'--- Processing: {os.path.basename(file_path)} ---')
            base_name = os.path.splitext(os.path.basename(file_path))[0]
            output_path = os.path.join(args.output_dir, f'{base_name}_swapped_to_{args.new_amplifier}.xlsx')
            _swap_amplifiers_on_file(file_path, output_path, args, amplifiers)
    elif os.path.isfile(args.input_probes):
        logger.info(f'Single file mode: Processing \'{args.input_probes}\'')
        base_name = os.path.splitext(os.path.basename(args.input_probes))[0]
        output_path = os.path.join(args.output_dir, f'{base_name}_swapped_to_{args.new_amplifier}.xlsx')
        _swap_amplifiers_on_file(args.input_probes, output_path, args, amplifiers)
    else:
        logger.error(f'Input path not found: \'{args.input_probes}\'')