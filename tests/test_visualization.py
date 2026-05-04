"""Tests for the SVG probe-map visualization."""
import re


def test_svg_has_unique_color_per_probe(tmp_path):
    """Phase 5.4: with N probes, the SVG should use N distinct fill colors
    (HSL hue stepping by 360/N) instead of recycling a 6-color palette."""
    from hcr_prober.visualization import generate_svg_probe_map

    n = 12
    probes = [{
        'pair_num': i + 1,
        'start_pos_on_sense': i * 60,
    } for i in range(n)]

    out = tmp_path / 'map.svg'
    generate_svg_probe_map(probes, seq_len=1000, amp='B1', gene='gene1',
                            out_path=str(out), window_size=52)

    content = out.read_text()
    fills = re.findall(r'fill="(hsl\([^"]+\))"', content)
    assert len(fills) == n, f'expected {n} probe rects, got {len(fills)}'
    assert len(set(fills)) == n, (
        f'expected {n} distinct colors, got {len(set(fills))}: {fills}'
    )
