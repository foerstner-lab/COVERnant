import unittest
import pandas as pd
from pandas.util.testing import assert_frame_equal
from covernantlib.extract import CoverageExtractor


class TestCoverageExtractor(unittest.TestCase):

    def setUp(self):
        self.coverage_extractor = CoverageExtractor(
            coordinate_file="", coverage_file="", output_prefix="",
            flip_reverse_strand=None, matrix_alignment=None,
            window_size=1, step_size=1)

    def test_init_coverage_dataframe_1(self):
        df = self.coverage_extractor._init_coverage_dataframe(3)
        self.assertListEqual(
            list(df.columns),
            ['Replicon', 'Start', 'End', 'Strand', 1, 2, 3])

    def test_init_coverage_dataframe_2(self):
        df = self.coverage_extractor._init_coverage_dataframe(
            3, [["test_rep", 1, 4, "+", 4, 5, 6]])
        self.assertListEqual(
            list(df.columns),
            ['Replicon', 'Start', 'End', 'Strand', 1, 2, 3])

    def test_generate_coverage_matrix(self):
        self.coverage_extractor._coverage_lists = [[9, 9, 9, 9]]
        self.coverage_extractor._coordinates = [
            {"replicon": "test", "start": 1, "end": 5, "strand": "+"}]
        self.coverage_extractor.generate_coverage_matrix()
        # assert_frame_equal(
        #     self.coverage_extractor._coverage_df,
        #     pd.DataFrame(
        #         [["test", 1, 5, "+", 9, 9, 9, 9]],
        #         columns=['Replicon', 'Start', 'End', 'Strand',
        #                  1, 2, 3, 4]))

    def test_align_coverages_1(self):
        """Alignment left"""
        self.coverage_extractor._matrix_alignment = "left"
        self.assertEqual(
            self.coverage_extractor._align_coverages([1.0, 2.0, 3.0], 10),
            [1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            
    def test_align_coverages_2(self):
        """Alignment center"""
        self.coverage_extractor._matrix_alignment = "center"
        self.assertEqual(
            self.coverage_extractor._align_coverages([1.0, 2.0, 3.0], 10),
            [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 0.0, 0.0, 0.0])

    def test_align_coverages_3(self):
        """Alignment right"""
        self.coverage_extractor._matrix_alignment = "right"
        self.assertEqual(
            self.coverage_extractor._align_coverages([1.0, 2.0, 3.0], 10),
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0])
