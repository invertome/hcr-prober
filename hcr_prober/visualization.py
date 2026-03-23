# hcr_prober/visualization.py
def generate_svg_probe_map(probes, seq_len, amp, gene, out_path):
    """Generates an SVG image showing the locations of probes on the target transcript."""
    # **UPDATED v1.9.5**: Compact layout with probes on one side.
    width, height = 800, 120
    padding = 20
    track_height = 10
    # Position the track lower in the canvas.
    track_y = height - padding - 40
    colors = ["#3498db", "#e74c3c", "#2ecc71", "#f1c40f", "#9b59b6", "#1abc9c"]

    svg = [f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">']
    svg.append('<style>.title { font: bold 20px sans-serif; } .label { font: 12px sans-serif; }</style>')
    svg.append(f'<text x="{width/2}" y="{padding+5}" text-anchor="middle" class="title">Probe Map: {gene} ({amp})</text>')

    # Draw the main transcript track
    svg.append(f'<rect x="{padding}" y="{track_y}" width="{width - 2 * padding}" height="{track_height}" fill="#bdc3c7" />')
    svg.append(f'<text x="{padding}" y="{track_y + track_height + 15}" class="label">0</text>')
    svg.append(f'<text x="{width - padding}" y="{track_y + track_height + 15}" text-anchor="end" class="label">{seq_len}</text>')

    scale = (width - 2 * padding) / seq_len if seq_len > 0 else 0
    sorted_probes = sorted(probes, key=lambda p: p['start_pos_on_sense'])

    for i, probe in enumerate(sorted_probes):
        start_px = padding + probe['start_pos_on_sense'] * scale
        probe_width = (probe['start_pos_rev'] + 52 - probe['start_pos_rev']) * scale
        # **UPDATED v1.9.5**: Draw all probes above the track, not alternating.
        probe_y = track_y - track_height
        color = colors[i % len(colors)]
        svg.append(f'<rect x="{start_px}" y="{probe_y}" width="{max(1, probe_width)}" height="{track_height}" fill="{color}">')
        # Add a tooltip to each probe rectangle
        svg.append(f'  <title>Pair {probe["pair_num"]}\nStart: {probe["start_pos_on_sense"]}</title>')
        svg.append('</rect>')

    svg.append('</svg>')
    with open(out_path, 'w', encoding='utf-8') as f: f.write('\n'.join(svg))