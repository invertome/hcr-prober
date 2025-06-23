# hcr_prober/visualization.py
def generate_svg_probe_map(probes, seq_len, amp, gene, out_path):
    """Creates an SVG image showing probe locations on the target sequence."""
    # Basic parameters
    width, height = 800, 150
    padding = 20
    track_y = height // 2
    track_height = 10
    text_y = track_y - 20
    colors = ["#3498db", "#e74c3c", "#2ecc71", "#f1c40f", "#9b59b6", "#1abc9c"]

    # SVG header
    svg = [f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">']
    # Background and Title
    svg.append('<style>.title { font: bold 20px sans-serif; } .label { font: 12px sans-serif; }</style>')
    svg.append(f'<text x="{width/2}" y="{padding+5}" text-anchor="middle" class="title">Probe Map: {gene} ({amp})</text>')
    # Main sequence track
    svg.append(f'<rect x="{padding}" y="{track_y}" width="{width - 2 * padding}" height="{track_height}" fill="#bdc3c7" />')
    svg.append(f'<text x="{padding}" y="{track_y + track_height + 15}" class="label">0</text>')
    svg.append(f'<text x="{width - padding}" y="{track_y + track_height + 15}" text-anchor="end" class="label">{seq_len}</text>')
    
    # Draw probes
    scale = (width - 2 * padding) / seq_len if seq_len > 0 else 0
    for i, probe in enumerate(probes):
        start_px = padding + probe['start_pos_on_sense'] * scale
        probe_width = 52 * scale # Window size
        probe_y = track_y - track_height if i % 2 == 0 else track_y + track_height
        color = colors[i % len(colors)]
        svg.append(f'<rect x="{start_px}" y="{probe_y}" width="{max(1, probe_width)}" height="{track_height}" fill="{color}">')
        svg.append(f'  <title>Pair {probe["pair_num"]}\nStart: {probe["start_pos_on_sense"]}</title>')
        svg.append('</rect>')

    svg.append('</svg>')
    with open(out_path, 'w') as f: f.write('\n'.join(svg))