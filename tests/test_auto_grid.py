import tempfile
import unittest
from pathlib import Path

import numpy as np

from auto_grid import DetectionError, detect_grids
from main import analyze_image
from phase2 import compare_grids


def synthetic_chip(seed=8, missing=0.2, distorted=True):
    rng = np.random.default_rng(seed)
    row, col = np.mgrid[:35, :35]
    grid = np.stack([80+20*col-1.4*row, 35+1.8*col+20.3*row], axis=-1)
    if distorted:
        grid[:, :, 0] += 1.4*np.sin(np.pi*col/34)*np.sin(np.pi*row/34)
        grid[:, :, 1] += 0.8*((col-17)/17)**2
    yy, xx = np.mgrid[:820, :820]
    image = 200+0.07*xx+0.1*yy+rng.normal(0, 2, (820, 820))
    for x, y in grid.reshape(-1, 2):
        if rng.random() < missing:
            continue
        amp = rng.uniform(400, 7000)
        x0, y0 = max(0, int(x)-8), max(0, int(y)-8)
        x1, y1 = min(820, int(x)+9), min(820, int(y)+9)
        image[y0:y1, x0:x1] += amp*np.exp(-((xx[y0:y1, x0:x1]-x)**2+(yy[y0:y1, x0:x1]-y)**2)/7)
    # Dust outside the lattice and a saturated artifact inside it.
    for x, y in [(770, 200), (750, 730), (367, 352)]:
        image += 65000*np.exp(-((xx-x)**2+(yy-y)**2)/30)
    return np.clip(image, 0, 65535).astype(np.uint16), grid


class AutomaticGridTests(unittest.TestCase):
    def test_missing_spots_skew_and_smooth_distortion(self):
        raw, expected = synthetic_chip()
        before = raw.copy()
        blocks = detect_grids(raw)
        self.assertEqual(len(blocks), 1)
        block = blocks[0]
        self.assertEqual(block['status'], 'OK')
        errors = np.linalg.norm(block['grid']-expected, axis=2)
        self.assertLess(float(np.median(errors)), 1.0)
        self.assertLess(float(np.max(errors)), 2.0)
        self.assertGreater(block['local_correction_max'], 0)
        np.testing.assert_array_equal(raw, before)
        rows, _ = analyze_image(raw, block['grid'])
        self.assertEqual(len(rows), 1225)

    def test_blank_is_not_a_grid(self):
        with self.assertRaises(DetectionError):
            detect_grids(np.full((800, 800), 500, dtype=np.uint16))

    def test_comparison_identical_grids(self):
        raw, grid = synthetic_chip(missing=0, distorted=False)
        with tempfile.TemporaryDirectory() as directory:
            summary = compare_grids(raw, grid, grid, 21, Path(directory)/'chip')
            self.assertEqual(summary['max_displacement_pixels'], 0)
            self.assertEqual(summary['mean_absolute_mean_difference'], 0)
            self.assertEqual(summary['mean_absolute_median_difference'], 0)


if __name__ == '__main__':
    unittest.main()
