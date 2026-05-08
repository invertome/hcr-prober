import subprocess, os, pytest, shutil

FIXTURE_DIR = os.path.join(os.path.dirname(__file__), 'fixtures')
FASTA = os.path.join(FIXTURE_DIR, 'test_sequence.fasta')
OUTPUT_DIR = '/tmp/hcr_prober_integration_test'

@pytest.fixture(autouse=True)
def cleanup():
    yield
    if os.path.exists(OUTPUT_DIR):
        shutil.rmtree(OUTPUT_DIR)

def run_design(*extra_args):
    """Helper to run hcr-prober design with common defaults."""
    cmd = [
        'hcr-prober', 'design', '-i', FASTA, '--amplifier', 'B1',
        '-o', OUTPUT_DIR, '--max-probes', '5',
        '--min-gc', '20', '--max-gc', '80', '--min-tm', '0', '--max-tm', '100',
        '--skip-5prime', '0',
    ] + list(extra_args)
    return subprocess.run(cmd, capture_output=True, text=True)

def test_design_basic():
    """Basic design command runs without errors."""
    result = run_design()
    assert result.returncode == 0, f"stderr: {result.stderr}"
    assert 'finished' in result.stderr.lower()

def test_design_produces_all_outputs():
    """Design should produce xlsx, fasta, svg, summary, and CSV."""
    result = run_design()
    assert result.returncode == 0
    found = {'xlsx': False, 'fasta': False, 'svg': False, 'txt': False, 'csv': False}
    for root, dirs, files in os.walk(OUTPUT_DIR):
        for f in files:
            if f.endswith('_order.xlsx'): found['xlsx'] = True
            if f.endswith('_probes.fasta'): found['fasta'] = True
            if f.endswith('_probe_map.svg'): found['svg'] = True
            if f.endswith('_summary.txt'): found['txt'] = True
            if f.endswith('_details.csv'): found['csv'] = True
    for fmt, present in found.items():
        assert present, f"Missing output: {fmt}"

def test_design_csv_has_thermo_columns():
    """Details CSV should have thermodynamic columns."""
    import pandas as pd
    result = run_design()
    assert result.returncode == 0
    csv_path = None
    for root, dirs, files in os.walk(OUTPUT_DIR):
        for f in files:
            if f.endswith('_details.csv'):
                csv_path = os.path.join(root, f)
    assert csv_path is not None
    df = pd.read_csv(csv_path)
    for col in ['gc_dn', 'gc_up', 'tm_dn', 'tm_up']:
        assert col in df.columns, f"Missing column: {col}"


def test_design_csv_has_start_and_end_position_columns():
    """Details CSV should include both start_pos_on_sense and
    end_pos_on_sense so each probe pair can be mapped back to its
    52-nt footprint on the input mRNA without external math.
    """
    import pandas as pd
    result = run_design()
    assert result.returncode == 0
    csv_path = None
    for root, dirs, files in os.walk(OUTPUT_DIR):
        for f in files:
            if f.endswith('_details.csv'):
                csv_path = os.path.join(root, f)
    assert csv_path is not None
    df = pd.read_csv(csv_path)
    assert 'start_pos_on_sense' in df.columns
    assert 'end_pos_on_sense' in df.columns
    if len(df):
        assert (df['end_pos_on_sense'] - df['start_pos_on_sense']).eq(52).all(), (
            f"end_pos_on_sense should always be start_pos_on_sense + 52 (window size); "
            f"diff was {(df['end_pos_on_sense'] - df['start_pos_on_sense']).unique()}"
        )

def test_design_with_structure_filter():
    """Design with structure filter flags should work."""
    result = run_design('--max-hairpin-dg', '-2.0', '--max-homodimer-dg', '-4.0')
    assert result.returncode == 0, f"stderr: {result.stderr}"

def test_design_with_custom_concentrations():
    """Design with custom Tm concentrations should work."""
    result = run_design('--na-conc', '100', '--mg-conc', '2.0')
    assert result.returncode == 0, f"stderr: {result.stderr}"

def test_design_summary_has_version():
    """Summary file should show the current package version."""
    from hcr_prober import __version__
    result = run_design()
    assert result.returncode == 0
    for root, dirs, files in os.walk(OUTPUT_DIR):
        for f in files:
            if f.endswith('_summary.txt'):
                with open(os.path.join(root, f)) as fh:
                    content = fh.read()
                    assert f'v{__version__}' in content
                    return
    pytest.fail("No summary file found")

def test_swap_still_works():
    """Swap command should still work."""
    result = run_design()
    assert result.returncode == 0
    xlsx = None
    for root, dirs, files in os.walk(OUTPUT_DIR):
        for f in files:
            if f.endswith('_order.xlsx'):
                xlsx = os.path.join(root, f)
    if xlsx is None:
        pytest.skip("No probes generated to swap")
    swap_dir = OUTPUT_DIR + '_swapped'
    result = subprocess.run(
        ['hcr-prober', 'swap', '--input-probes', xlsx, '--new-amplifier', 'B3',
         '--output-dir', swap_dir],
        capture_output=True, text=True
    )
    assert result.returncode == 0
    if os.path.exists(swap_dir):
        shutil.rmtree(swap_dir)

def test_version_in_help():
    """Help output should show the current package version."""
    from hcr_prober import __version__
    result = subprocess.run(['hcr-prober', '--help'], capture_output=True, text=True)
    assert __version__ in result.stdout


@pytest.mark.parametrize('subcommand', ['design', 'isoform-split', 'swap'])
def test_subcommand_help_formats_cleanly(subcommand):
    """Subcommand -h must print help and exit 0. argparse uses % for help-string
    interpolation; any literal % in a help string must be escaped as %% or
    formatting blows up only when the help formatter actually runs (which the
    rest of the test suite never triggers)."""
    result = subprocess.run(['hcr-prober', subcommand, '-h'], capture_output=True, text=True)
    assert result.returncode == 0, (
        f"{subcommand} -h exited {result.returncode}\nstderr: {result.stderr}"
    )
    assert 'ValueError' not in result.stderr
    assert 'usage:' in result.stdout
