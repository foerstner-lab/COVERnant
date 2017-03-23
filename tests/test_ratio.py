import unittest
from covernantlib.ratio import CoverageRatioCalculator
from argparse import Namespace
import numpy as np


class TestCoverageRatioCalculator(unittest.TestCase):

    def setUp(self):
        self._cov_ratio_calculator = CoverageRatioCalculator(
            Namespace(denominator_bam_file="", numerator_bam_file="",
                      output_prefix="", window_size=1, step_size=1,
                      factor=1, factor_numerator=None, factor_denominator=None,
                      keep_zero_coverage=True, denominator_name="",
                      numerator_name="", ratio_name="", paired_end=False))

    def test_calc_averaged_coverages_1(self):
        averaged_coverages = (
            self._cov_ratio_calculator._calc_averaged_coverages(
                {"chrom1": [1.0, 1.0, 1.0, 2.0, 2.0]}))
        assert list(averaged_coverages["chrom1"]) == [1.0, 1.0, 1.0, 2.0, 2.0]

    def test_sliding_windows_average_1(self):
        self._cov_ratio_calculator._window_size = 1
        self._cov_ratio_calculator.step_size = 1
        assert list(self._cov_ratio_calculator._sliding_windows_average(
            [1.0, 1.0, 1.0, 2.0, 2.0])) == [1.0, 1.0, 1.0, 2.0, 2.0]

    def test_sliding_windows_average_2(self):
        self._cov_ratio_calculator._window_size = 3
        self._cov_ratio_calculator.step_size = 1
        assert list(self._cov_ratio_calculator._sliding_windows_average(
            [1.0, 1.0, 1.0, 4.0, 1.0])) == [0.0, 1.0, 2.0, 2.0, 0.0]
        
    def test_sliding_windows_average_3(self):
        self._cov_ratio_calculator._window_size = 1
        self._cov_ratio_calculator._step_size = 2
        assert list(self._cov_ratio_calculator._sliding_windows_average(
            [1.0, 1.0, 1.0, 4.0, 1.0])) == [1.0, 0, 1.0, 0.0, 1.0]
        
    def test_ratio_1(self):
        assert self._cov_ratio_calculator._ratio(10.0, 5.0) == 2

    def test_ratio_2(self):
        assert self._cov_ratio_calculator._ratio(5.0, 10.0) == 0.5

    def test_ratio_3(self):
        assert self._cov_ratio_calculator._ratio(0.0, 10.0) == 0

    def test_ratio_4(self):
        assert self._cov_ratio_calculator._ratio(10.0, 0.0) == 0
        
    def test_calc_normalization_factors_1(self):
        self._cov_ratio_calculator.no_of_mapped_bases_numerator = 10
        self._cov_ratio_calculator.no_of_mapped_bases_denominator = 5
        self._cov_ratio_calculator.calc_normalization_factors()
        assert self._cov_ratio_calculator._numerator_normalization_factor == (
            100000.0)
        assert (
            self._cov_ratio_calculator._denominator_normalization_factor == (
                200000.0))
        
    def test_calc_normalization_factors_2(self):
        """In case numerator and denominator factors are given the two number
        are used for the ration calculation.

        """
        self._cov_ratio_calculator._numerator_factor_given = 12
        self._cov_ratio_calculator._denominator_factor_given = 6
        self._cov_ratio_calculator.calc_normalization_factors()
        assert self._cov_ratio_calculator._numerator_normalization_factor == 12
        assert (
            self._cov_ratio_calculator._denominator_normalization_factor == 6)

    def test_calc_normalization_factors_3(self):
        self._cov_ratio_calculator._numerator_factor_given = 12
        self._cov_ratio_calculator._denominator_factor_given = 6
        self._cov_ratio_calculator.calc_normalization_factors()
        assert self._cov_ratio_calculator._numerator_normalization_factor == 12
        assert (
            self._cov_ratio_calculator._denominator_normalization_factor == 6)
        
    def test_compare_coverages_1(self):
        assert self._cov_ratio_calculator._compare_coverages(
            "chrom1", {"chrom1": [2, 4, 8]}, {"chrom1": [4, 8, 48]}) == [
                2, 2, 6]

    def test_normalize_coverages_1(self):
        # Denominator
        self._cov_ratio_calculator._denominator_normalization_factor = 5
        self._cov_ratio_calculator.coverage_denominator = {
            "chrom": np.array([1.0, 1.0, 1.0, 1.0, 1.0])}
        self._cov_ratio_calculator.coverage_denominator_forward = {
            "chrom": np.array([1.0, 1.0, 1.0, 1.0, 1.0])}
        self._cov_ratio_calculator.coverage_denominator_reverse = {
            "chrom": np.array([1.0, 1.0, 1.0, 1.0, 1.0])}
        # Numerator
        self._cov_ratio_calculator._numerator_normalization_factor = 10
        self._cov_ratio_calculator.coverage_numerator = {
            "chrom": np.array([1.0, 1.0, 1.0, 1.0, 1.0])}
        self._cov_ratio_calculator.coverage_numerator_forward = {
            "chrom": np.array([1.0, 1.0, 1.0, 1.0, 1.0])}
        self._cov_ratio_calculator.coverage_numerator_reverse = {
            "chrom": np.array([1.0, 1.0, 1.0, 1.0, 1.0])}
        self._cov_ratio_calculator.normalize_coverages()
        assert list(
            self._cov_ratio_calculator.coverage_denominator_normalized[
                "chrom"]) == [5.0, 5.0, 5.0, 5.0, 5.0]
        assert list(
            self._cov_ratio_calculator.coverage_denominator_reverse_normalized[
                "chrom"]) == [5.0, 5.0, 5.0, 5.0, 5.0]
        assert list(
            self._cov_ratio_calculator.coverage_denominator_forward_normalized[
                "chrom"]) == [5.0, 5.0, 5.0, 5.0, 5.0]
        assert list(self._cov_ratio_calculator.coverage_numerator_normalized[
            "chrom"]) == [10.0, 10.0, 10.0, 10.0, 10.0]
        assert list(
            self._cov_ratio_calculator.coverage_numerator_reverse_normalized[
                "chrom"]) == [10.0, 10.0, 10.0, 10.0, 10.0]
        assert list(
            self._cov_ratio_calculator.coverage_numerator_forward_normalized[
                "chrom"]) == [10.0, 10.0, 10.0, 10.0, 10.0]
