import csv
import tempfile
import unittest
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image

from main import (analyze_image, create_qc_image, generate_grid, load_tiff,
                  make_display_image, prepare_rois, save_results,
                  calculate_foreground_statistics, calculate_threshold_percentages)


class AnalysisTests(unittest.TestCase):
    def test_bilinear_skew(self):
        corners = np.array([[15, 15], [700, 30], [30, 700], [730, 740]])
        grid = generate_grid(corners)
        np.testing.assert_allclose(grid[0, 0], corners[0])
        np.testing.assert_allclose(grid[0, -1], corners[1])
        np.testing.assert_allclose(grid[-1, 0], corners[2])
        np.testing.assert_allclose(grid[-1, -1], corners[3])
        np.testing.assert_allclose(grid[17, 17], corners.mean(axis=0))

    def test_statistics(self):
        roi = np.arange(121, dtype=np.uint16).reshape(11, 11) + 1000
        stats = calculate_foreground_statistics(roi)
        self.assertEqual(stats['Mean'], 1060)
        self.assertEqual(stats['Total_Intensity'], 128260)
        self.assertAlmostEqual(stats['SD'], np.sqrt(121*122/12))
        percentages = calculate_threshold_percentages(roi, {'BG_Mean': 1050, 'BG_SD': 10})
        self.assertAlmostEqual(percentages['Pct_BG_1SD'], 60 / 121 * 100)
        self.assertAlmostEqual(percentages['Pct_BG_2SD'], 50 / 121 * 100)

    def test_full_outputs_and_neighbor_exclusion(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            grid = generate_grid([[10.5, 10.5], [418.5, 10.5], [10.5, 418.5], [418.5, 418.5]])
            raw = np.full((430, 430), 100, dtype=np.uint16)
            boxes, occupied = prepare_rois(grid, raw.shape)
            raw[occupied] = 60000
            tifffile.imwrite(directory / 'chip.tif', raw)
            loaded = load_tiff(directory / 'chip.tif')
            display = make_display_image(loaded)
            np.testing.assert_array_equal(raw, loaded)
            self.assertEqual(loaded.dtype, np.uint16)
            self.assertFalse(loaded.flags.writeable)
            self.assertFalse(np.shares_memory(display, loaded))
            results, measured = analyze_image(loaded, grid)
            self.assertEqual(len(results), 1225)
            self.assertEqual(results[-1]['Spot'], 'R35C35')
            self.assertEqual(results[0]['Center_X'], 11)
            for row in results:
                self.assertEqual(row['Mean'], 60000)
                self.assertEqual(row['Total_Intensity'], 60000 * 121)
                self.assertEqual(row['BG_Mean'], 100)
                self.assertEqual(row['BG_SD'], 0)
                self.assertEqual(row['Pct_BG_2SD'], 100)
            save_results(results, directory / 'results.csv')
            with open(directory / 'results.csv', encoding='utf-8-sig') as handle:
                self.assertEqual(len(list(csv.DictReader(handle))), 1225)
            create_qc_image(display, measured, directory / 'QC.png')
            with Image.open(directory / 'QC.png') as image:
                self.assertEqual(image.size, (430, 430))
                for x0, y0, x1, y1 in boxes:
                    self.assertEqual(image.getpixel((x0, y0)), (0, 255, 0))
                    self.assertEqual(image.getpixel((x1-1, y1-1)), (0, 255, 0))

    def test_invalid_geometry_and_boundaries(self):
        with self.assertRaises(ValueError):
            generate_grid([[10, 10], [700, 10], [700, 700], [10, 700]])
        with self.assertRaises(ValueError):
            prepare_rois(generate_grid([[0, 0], [700, 0], [0, 700], [700, 700]]), (710, 710))
        with self.assertRaises(ValueError):
            prepare_rois(generate_grid([[10, 10], [350, 10], [10, 350], [350, 350]]), (370, 370))


if __name__ == '__main__':
    unittest.main()
