"""Fluorescence microarray measurement with manual or automatic 35×35 grids."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image, ImageDraw

ROWS = COLS = 35
ROI_SIZE = 11
FIELDS = ['Spot', 'Row', 'Column', 'Center_X', 'Center_Y', 'Mean', 'Median',
          'SD', 'CV', 'Min', 'Max', 'Total_Intensity', 'BG_Mean', 'BG_Median',
          'BG_SD', 'BG_CV', 'Pct_BG_1SD', 'Pct_BG_2SD']


def load_tiff(path):
    """Keep original dtype; reject ambiguous stacks and color images."""
    raw = tifffile.imread(path)
    if raw.ndim != 2 or raw.dtype.kind not in 'uif':
        raise ValueError('A single 2D grayscale TIFF is required (no stacks/RGB).')
    if not np.isfinite(raw).all():
        raise ValueError('TIFF contains NaN or infinite intensities.')
    raw.setflags(write=False)
    return raw


def make_display_image(raw):
    """Percentile stretch on a separate array, for visualization only."""
    low, high = np.percentile(raw, [1, 99.8])
    if high <= low:
        low, high = float(raw.min()), float(raw.max())
    if high <= low:
        return np.zeros(raw.shape, dtype=np.uint8)
    return (np.clip((raw.astype(np.float64) - low) / (high - low), 0, 1) * 255).astype(np.uint8)


def validate_corners(corners):
    points = np.asarray(corners, dtype=float)
    if points.shape != (4, 2) or not np.isfinite(points).all():
        raise ValueError('Exactly four finite (x, y) corners are required: TL, TR, BL, BR.')
    polygon = points[[0, 1, 3, 2]]
    edges = np.roll(polygon, -1, axis=0) - polygon
    following = np.roll(edges, -1, axis=0)
    cross = edges[:, 0] * following[:, 1] - edges[:, 1] * following[:, 0]
    if not np.all(cross > 1e-6):
        raise ValueError('Corners must form a convex TL → TR → BR → BL polygon in image coordinates.')
    return points


def generate_grid(corners):
    tl, tr, bl, br = validate_corners(corners)
    v, u = np.meshgrid(np.linspace(0, 1, ROWS), np.linspace(0, 1, COLS), indexing='ij')
    return ((1-u)[..., None]*(1-v)[..., None]*tl + u[..., None]*(1-v)[..., None]*tr
            + (1-u)[..., None]*v[..., None]*bl + u[..., None]*v[..., None]*br)


def roi_bounds(center, shape):
    """Half-up rounding; slices are exclusive at the upper bound."""
    x, y = np.floor(np.asarray(center) + 0.5).astype(int)
    x0, y0 = x - 5, y - 5
    x1, y1 = x + 6, y + 6
    if x0 < 0 or y0 < 0 or x1 > shape[1] or y1 > shape[0]:
        raise ValueError(f'11×11 ROI at ({x}, {y}) extends outside the image.')
    return x0, y0, x1, y1


def extract_roi(raw, bounds):
    x0, y0, x1, y1 = bounds
    return raw[y0:y1, x0:x1]


def calculate_foreground_statistics(roi):
    values = roi.astype(np.float64)
    mean, sd = float(values.mean()), float(values.std(ddof=1))
    return dict(Mean=mean, Median=float(np.median(values)), SD=sd,
                CV=sd / mean * 100 if mean != 0 else float('nan'),
                Min=float(values.min()), Max=float(values.max()),
                Total_Intensity=float(values.sum()))


def calculate_background_statistics(values):
    if values.size < 2:
        raise ValueError('Fewer than two background pixels remain; increase --background-size.')
    values = values.astype(np.float64)
    mean, sd = float(values.mean()), float(values.std(ddof=1))
    return dict(BG_Mean=mean, BG_Median=float(np.median(values)), BG_SD=sd,
                BG_CV=sd / mean * 100 if mean != 0 else float('nan'))


def calculate_threshold_percentages(roi, background):
    return {f'Pct_BG_{n}SD': float(np.count_nonzero(
        roi > background['BG_Mean'] + n * background['BG_SD']) / roi.size * 100)
        for n in (1, 2)}


def prepare_rois(grid, shape):
    bounds = [roi_bounds(center, shape) for center in grid.reshape(-1, 2)]
    occupied = np.zeros(shape, dtype=bool)
    for x0, y0, x1, y1 in bounds:
        if occupied[y0:y1, x0:x1].any():
            raise ValueError('Foreground ROIs overlap. Check corners and spot spacing.')
        occupied[y0:y1, x0:x1] = True
    return bounds, occupied


def analyze_image(raw, grid, background_size=21):
    if background_size <= ROI_SIZE or background_size % 2 != 1:
        raise ValueError('Background size must be odd and greater than 11.')
    bounds, occupied = prepare_rois(grid, raw.shape)
    results = []
    half = background_size // 2
    for index, box in enumerate(bounds):
        row, col = divmod(index, COLS)
        x0, y0, x1, y1 = box
        x, y = x0 + 5, y0 + 5
        bx0, by0 = max(0, x-half), max(0, y-half)
        bx1, by1 = min(raw.shape[1], x+half+1), min(raw.shape[0], y+half+1)
        bg = calculate_background_statistics(
            raw[by0:by1, bx0:bx1][~occupied[by0:by1, bx0:bx1]])
        roi = extract_roi(raw, box)
        result = dict(Spot=f'R{row+1:02d}C{col+1:02d}', Row=row+1, Column=col+1,
                      Center_X=x, Center_Y=y)
        result.update(calculate_foreground_statistics(roi))
        result.update(bg)
        result.update(calculate_threshold_percentages(roi, bg))
        results.append(result)
    return results, bounds


def select_four_corners(display):
    """Interactive selection with cursor-centered wheel zoom and right-drag pan."""
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Button
    from matplotlib.patches import Rectangle

    fig, ax = plt.subplots(figsize=(12, 9))
    fig.subplots_adjust(bottom=0.16)
    ax.imshow(display, cmap='gray', vmin=0, vmax=255, interpolation='nearest')
    full_xlim, full_ylim = ax.get_xlim(), ax.get_ylim()
    drag = [None]
    fig.text(0.5, 0.105, 'Wheel: zoom at cursor | Right drag: pan | Left click: select corner',
             ha='center', fontsize=10)
    points, overlays = [], []
    approved = [False]
    labels = ['TL', 'TR', 'BL', 'BR']

    def redraw(message=None):
        for artist in overlays:
            artist.remove()
        overlays.clear()
        for label, (x, y) in zip(labels, points):
            overlays.extend(ax.plot(x, y, 'r+', markersize=10))
            overlays.append(ax.text(x+6, y, label, color='yellow'))
        if len(points) == 4:
            try:
                grid = generate_grid(points)
                boxes, _ = prepare_rois(grid, display.shape)
                for x0, y0, x1, y1 in boxes:
                    patch = Rectangle((x0-0.5, y0-0.5), x1-x0, y1-y0,
                                      fill=False, edgecolor='lime', linewidth=0.45)
                    ax.add_patch(patch)
                    overlays.append(patch)
                message = 'Check all 1225 ROIs, then click Analyze. Undo/Reset to revise.'
            except ValueError as error:
                message = str(error)
        ax.set_title(message or f'Click {labels[len(points)]} spot center ({len(points)}/4). Wheel to zoom.')
        fig.canvas.draw_idle()

    def click(event):
        if event.inaxes != ax or event.button != 1 or len(points) >= 4:
            return
        if ax.get_navigate_mode() is not None:
            return
        if event.xdata is not None and event.ydata is not None:
            points.append((float(event.xdata), float(event.ydata)))
            redraw()

    def undo(event):
        if points:
            points.pop()
        redraw()

    def scroll(event):
        if event.inaxes != ax or event.xdata is None or event.ydata is None:
            return
        if ax.get_navigate_mode() is not None or drag[0] is not None:
            return
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        # Limit magnification so individual pixels remain easy to inspect.
        factor = 1.25 ** (-float(event.step))
        width = abs(xlim[1] - xlim[0])
        target = np.clip(width * factor, 5, display.shape[1] * 4)
        factor = target / width
        ax.set_xlim([event.xdata + (edge-event.xdata)*factor for edge in xlim])
        ax.set_ylim([event.ydata + (edge-event.ydata)*factor for edge in ylim])
        fig.canvas.draw_idle()

    def start_pan(event):
        if event.inaxes == ax and event.button == 3 and ax.get_navigate_mode() is None:
            drag[0] = (event.x, event.y, ax.get_xlim(), ax.get_ylim())

    def move_pan(event):
        if drag[0] is None or event.x is None or event.y is None:
            return
        start_x, start_y, xlim, ylim = drag[0]
        dx = (event.x-start_x) * (xlim[1]-xlim[0]) / ax.bbox.width
        dy = (event.y-start_y) * (ylim[1]-ylim[0]) / ax.bbox.height
        ax.set_xlim(xlim[0]-dx, xlim[1]-dx)
        ax.set_ylim(ylim[0]-dy, ylim[1]-dy)
        fig.canvas.draw_idle()

    def stop_pan(event):
        drag[0] = None

    def fit_image(event):
        drag[0] = None
        ax.set_xlim(full_xlim)
        ax.set_ylim(full_ylim)
        fig.canvas.draw_idle()

    def reset(event):
        points.clear()
        redraw()

    def approve(event):
        if len(points) != 4:
            redraw('Select all four corners before analysis.')
            return
        try:
            prepare_rois(generate_grid(points), display.shape)
        except ValueError as error:
            redraw(str(error))
            return
        approved[0] = True
        plt.close(fig)

    buttons = []
    for left, label, callback in [(0.12, 'Undo', undo), (0.32, 'Reset', reset),
                                  (0.52, 'Fit image', fit_image), (0.72, 'Analyze', approve)]:
        button = Button(fig.add_axes([left, 0.04, 0.16, 0.05]), label)
        button.on_clicked(callback)
        buttons.append(button)
    fig.canvas.mpl_connect('button_press_event', click)
    fig.canvas.mpl_connect('scroll_event', scroll)
    fig.canvas.mpl_connect('button_press_event', start_pan)
    fig.canvas.mpl_connect('motion_notify_event', move_pan)
    fig.canvas.mpl_connect('button_release_event', stop_pan)
    redraw()
    plt.show()
    if not approved[0]:
        raise RuntimeError('Selection cancelled; analysis was not approved.')
    return validate_corners(points)


def save_results(results, path):
    with open(path, 'w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(results)


def create_qc_image(display, bounds, path):
    """Full-resolution PNG: draw exactly the measured inclusive pixel extents."""
    image = Image.fromarray(display).convert('RGB')
    draw = ImageDraw.Draw(image)
    for x0, y0, x1, y1 in bounds:
        draw.rectangle((x0, y0, x1-1, y1-1), outline=(0, 255, 0), width=1)
    image.save(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path, default=Path('results'))
    parser.add_argument('--background-size', type=int, default=21)
    parser.add_argument('--grid', choices=['manual', 'auto'], default='manual')
    parser.add_argument('--fallback', choices=['manual', 'stop'], default='manual',
                        help='Automatic detection failure: open manual selection or save QC and stop.')
    parser.add_argument('--block', type=int, help='Analyze one auto-detected block (1-based reading order).')
    parser.add_argument('--compare-manual', type=Path, help='Phase 1 corners JSON to compare with the nearest auto block.')
    args = parser.parse_args()
    try:
        if args.background_size <= 11 or args.background_size % 2 != 1:
            raise ValueError('--background-size must be odd and greater than 11.')
        raw = load_tiff(args.input)
        print(f'Loaded {raw.shape[1]}×{raw.shape[0]}, dtype={raw.dtype}, raw range={raw.min()}..{raw.max()}')
        display = make_display_image(raw)
        if args.grid == 'auto':
            from phase2 import run_auto
            run_auto(raw, display, args)
            return
        if args.block is not None or args.compare_manual:
            raise ValueError('--block and --compare-manual require --grid auto.')
        corners = select_four_corners(display)
        grid = generate_grid(corners)
        results, bounds = analyze_image(raw, grid, args.background_size)
        args.output.mkdir(parents=True, exist_ok=True)
        prefix = args.output / args.input.stem
        save_results(results, str(prefix) + '_results.csv')
        create_qc_image(display, bounds, str(prefix) + '_QC.png')
        metadata = dict(input=str(args.input.resolve()), shape=list(raw.shape), dtype=str(raw.dtype),
                        corners=dict(zip(['TL', 'TR', 'BL', 'BR'], corners.tolist())),
                        grid=[ROWS, COLS], roi_size=ROI_SIZE, background_size=args.background_size,
                        sd_ddof=1, coordinate_origin='zero-based pixel center', rounding='floor(value + 0.5)')
        Path(str(prefix) + '_corners.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
        print(f'Saved {len(results)} spots to {args.output.resolve()}')
    except (ValueError, RuntimeError, OSError, tifffile.TiffFileError) as error:
        parser.exit(1, f'Error: {error}\n')


if __name__ == '__main__':
    main()
