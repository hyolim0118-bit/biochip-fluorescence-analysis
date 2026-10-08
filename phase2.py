"""Automatic workflow and comparison; Phase 1 measurement functions are reused."""
import csv
import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

from auto_grid import DetectionError, detect_grids


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False), encoding='utf-8')


def block_metadata(block):
    return {key: value for key, value in block.items() if key not in ('grid', 'node_distance')}


def save_detection_qc(display, blocks, path):
    """Diagnostic only; successful measurement also saves Phase 1 exact-bounds QC."""
    image = Image.fromarray(display).convert('RGB')
    draw = ImageDraw.Draw(image)
    for block in blocks:
        color = (0, 255, 0) if block['status'] == 'OK' else (255, 180, 0)
        for center in block['grid'].reshape(-1, 2):
            x, y = np.floor(center+0.5).astype(int)
            draw.rectangle((int(x-5), int(y-5), int(x+5), int(y+5)), outline=color)
        x, y = block['grid'][0, 0]
        draw.text((int(x), max(0, int(y)-24)), f"Block {block['block']:02d} {block['status']}", fill=color)
    image.save(path)


def compare_grids(raw, automatic, manual, background_size, prefix):
    from main import analyze_image
    auto_rows, _ = analyze_image(raw, automatic, background_size)
    manual_rows, _ = analyze_image(raw, manual, background_size)
    displacement = automatic-manual
    rows = []
    for index, (auto, reference) in enumerate(zip(auto_rows, manual_rows)):
        row, col = divmod(index, 35)
        dx, dy = displacement[row, col]
        measured_dx = auto['Center_X']-reference['Center_X']
        measured_dy = auto['Center_Y']-reference['Center_Y']
        rows.append(dict(Spot=auto['Spot'], Auto_Center_X=float(automatic[row, col, 0]),
                         Auto_Center_Y=float(automatic[row, col, 1]),
                         Manual_Center_X=float(manual[row, col, 0]), Manual_Center_Y=float(manual[row, col, 1]),
                         Delta_X=float(dx), Delta_Y=float(dy), Grid_Displacement=float(np.hypot(dx, dy)),
                         ROI_Delta_X=measured_dx, ROI_Delta_Y=measured_dy,
                         Auto_Mean=auto['Mean'], Manual_Mean=reference['Mean'],
                         Mean_Difference=auto['Mean']-reference['Mean'],
                         Auto_Median=auto['Median'], Manual_Median=reference['Median'],
                         Median_Difference=auto['Median']-reference['Median']))
    with open(str(prefix)+'_comparison.csv', 'w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    distances = np.linalg.norm(displacement, axis=2)
    summary = dict(mean_displacement_pixels=float(distances.mean()),
                   median_displacement_pixels=float(np.median(distances)),
                   max_displacement_pixels=float(distances.max()),
                   mean_absolute_mean_difference=float(np.mean([abs(r['Mean_Difference']) for r in rows])),
                   mean_absolute_median_difference=float(np.mean([abs(r['Median_Difference']) for r in rows])))
    write_json(str(prefix)+'_comparison_summary.json', summary)
    return summary


def run_auto(raw, display, args):
    from main import (analyze_image, create_qc_image, generate_grid, prepare_rois,
                      save_results, select_four_corners)
    destination = args.output / 'auto'
    destination.mkdir(parents=True, exist_ok=True)
    overview_prefix = destination / args.input.stem
    try:
        blocks = detect_grids(raw)
    except DetectionError as error:
        write_json(str(overview_prefix)+'_detection.json', dict(status='GRID_MISMATCH', reason=str(error)))
        Image.fromarray(display).save(str(overview_prefix)+'_auto_grid_QC.png')
        if args.fallback == 'stop':
            raise DetectionError(f'{error} Retry with --grid manual.')
        print(f'GRID_MISMATCH: {error}\nOpening Phase 1 manual selection.')
        corners = select_four_corners(display)
        grid = generate_grid(corners)
        results, bounds = analyze_image(raw, grid, args.background_size)
        prefix = destination / (args.input.stem+'_manual_fallback')
        save_results(results, str(prefix)+'_results.csv')
        create_qc_image(display, bounds, str(prefix)+'_QC.png')
        write_json(str(prefix)+'_corners.json', dict(corners=dict(zip(['TL', 'TR', 'BL', 'BR'], corners.tolist())),
                                                    shape=list(raw.shape), source='manual_fallback'))
        return
    # Validate before publishing any automatic intensity results.
    for block in blocks:
        try:
            prepare_rois(block['grid'], raw.shape)
        except ValueError as error:
            block['status'] = 'OUT_OF_BOUNDS' if 'outside' in str(error) else 'GRID_MISMATCH'
            block['reason'] = str(error)
    save_detection_qc(display, blocks, str(overview_prefix)+'_auto_grid_QC.png')
    write_json(str(overview_prefix)+'_detection.json', dict(blocks=[block_metadata(b) for b in blocks]))
    selected = blocks
    if args.block is not None:
        selected = [b for b in blocks if b['block'] == args.block]
        if not selected:
            raise ValueError(f'--block must be in 1..{len(blocks)}.')
    manual_grid = None
    compare_block = None
    if args.compare_manual:
        data = json.loads(args.compare_manual.read_text(encoding='utf-8'))
        if not isinstance(data, dict) or not isinstance(data.get('corners'), dict) or not all(
                name in data['corners'] for name in ('TL', 'TR', 'BL', 'BR')):
            raise ValueError('Manual reference must contain TL/TR/BL/BR in a corners object.')
        if 'shape' in data and data['shape'] != list(raw.shape):
            raise ValueError('Manual reference image shape differs from input TIFF.')
        manual_grid = generate_grid([data['corners'][name] for name in ('TL', 'TR', 'BL', 'BR')])
        prepare_rois(manual_grid, raw.shape)
        compare_block = min(selected, key=lambda b: np.linalg.norm(np.asarray(b['center'])-manual_grid.mean(axis=(0, 1))))['block']
    skipped = []
    for block in selected:
        number, grid = block['block'], block['grid']
        directory = destination / f'block{number:02d}'
        directory.mkdir(parents=True, exist_ok=True)
        prefix = directory / args.input.stem
        print(f"Block {number:02d}: {block['status']}, matched={block['matched_fraction']:.1%}, RMS={block['rms_pixels']:.2f}px")
        save_detection_qc(display, [block], str(prefix)+'_auto_grid_QC.png')
        write_json(str(prefix)+'_grid.json', dict(**block_metadata(block), grid=grid.tolist(),
                   coordinate_origin='zero-based pixel center', rounding='floor(value + 0.5)',
                   roi_size=11, shape=list(raw.shape), dtype=str(raw.dtype)))
        # Keep QC separate from the unchanged measurement CSV.
        with open(str(prefix)+'_confidence.csv', 'w', newline='', encoding='utf-8-sig') as handle:
            writer = csv.writer(handle)
            writer.writerow(['Spot', 'Grid_Status', 'Spot_Status', 'Candidate_Distance_Pixels'])
            for index, distance in enumerate(block['node_distance'].ravel()):
                row, col = divmod(index, 35)
                status = 'OK' if distance < min(block['pitch_x'], block['pitch_y'])*0.2 else 'LOW_CONFIDENCE'
                writer.writerow([f'R{row+1:02d}C{col+1:02d}', block['status'], status, distance])
        if manual_grid is not None and number == compare_block and block['status'] not in ('OUT_OF_BOUNDS', 'GRID_MISMATCH'):
            summary = compare_grids(raw, grid, manual_grid, args.background_size, prefix)
            print(f"Manual comparison: mean displacement={summary['mean_displacement_pixels']:.3f}px")
        source = 'auto'
        if block['status'] != 'OK':
            if args.fallback == 'stop':
                skipped.append(number)
                continue
            print(f'Block {number:02d}: manual fallback. Select this block in the image.')
            corners = select_four_corners(display)
            grid = generate_grid(corners)
            source = 'manual_fallback'
            write_json(str(prefix)+'_manual_corners.json', dict(shape=list(raw.shape),
                       corners=dict(zip(['TL', 'TR', 'BL', 'BR'], corners.tolist()))))
        results, bounds = analyze_image(raw, grid, args.background_size)
        save_results(results, str(prefix)+'_results.csv')
        create_qc_image(display, bounds, str(prefix)+'_QC.png')
        write_json(str(prefix)+'_measurement.json', dict(source=source, grid=grid.tolist(),
                   background_size=args.background_size, roi_size=11, sd_ddof=1, input=str(args.input.resolve())))
    if skipped:
        raise DetectionError(f'Low-confidence blocks {skipped} saved as QC only. Use --grid manual or --fallback manual. Other OK blocks were measured.')
    print(f'Automatic workflow outputs: {destination.resolve()}')
