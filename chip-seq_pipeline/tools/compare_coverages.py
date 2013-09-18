#!/usr/bin/env python

__description__ = ""
__author__ = "Konrad Foerstner <konrad@foerstner.org>"
__copyright__ = "2013 by Konrad Foerstner <konrad@foerstner.org>"
__license__ = "ISC license"
__email__ = "konrad@foerstner.org"
__version__ = ""

import argparse
import sys
sys.path.append("..")
from libs.coveragecalculator import CoverageCalculator
import pysam

def main():
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument("bam_file_control")
    parser.add_argument("bam_file_chip")
    parser.add_argument("--output", "-o", dest="output_prefix",
                        default="chip_seq.wig", required=False)
    parser.add_argument("--window_size", type=int, default=None,
                        help="Window size for sliding window average calculation.")
    parser.add_argument("--step_size", type=int, default=1,
                        help="Step size for sliding window average calculation."
                        " Default is 1.")
    parser.add_argument("--factor", type=float, default=None,
                        help="A factor the final ratio is multiplied with.")
    parser.add_argument("--keep_zero_coverage", default=False, action="store_true",
                        help="Also write positions with that have coverage of 0. "
                        "Default is to discard those.")

    # TODO
    # parser.add_argument("--pseudocount",
    #                    help="Use pseudocounts to initate the coverage vector.")
    # parser.add_argument("--uniquely", help="Use uniquely mapped reads only")
    # parser.add_argument("--zero_div", help="What to do if the divisor is zero - 'zero', "minus", "addone"")
    # 'zero' - the value becomes zero
    # 'minus' - the value becomes
    args = parser.parse_args()
    coverage_comparer = CoverageComparer(
        args.bam_file_control, args.bam_file_chip, args.output_prefix,
        args.window_size, args.step_size, args.factor, args.keep_zero_coverage)
    coverage_comparer.calc_coverages()
    coverage_comparer.print_no_aligned_reads()
    coverage_comparer.calc_combined_factor()
    coverage_comparer.write_chip_and_control_wiggle_files()
    coverage_comparer.compare()
    coverage_comparer.write_ratio_wiggle_file()

    # multi smallest or 1M
    # print total number of mapped reads

class CoverageComparer(object):

    def __init__(self, bam_file_control, bam_file_chip, output_prefix,
                 window_size, step_size, factor, keep_zero_coverage):
        self._bam_file_control = bam_file_control
        self._bam_file_chip = bam_file_chip
        self._output_prefix = output_prefix
        if window_size % 2 == 0:
            sys.stderr.write("Error. Window size must be an odd number!\n")
            sys.exit(2)
        self._window_size = window_size
        self._step_size = step_size
        self._factor = factor
        self._combine_factor = None
        self._keep_zero_coverage = keep_zero_coverage

    def calc_coverages(self):
        self._print_file_names()
        self.no_of_mapped_reads_chip = self._count_no_of_mapped_reads(
            self._bam_file_chip)
        self.no_of_mapped_reads_control = self._count_no_of_mapped_reads(
            self._bam_file_control)
        self.coverage_control = self._prepare_coverage(self._bam_file_control)
        self.coverage_chip = self._prepare_coverage(self._bam_file_chip)

    def print_no_aligned_reads(self):
        print("Number of mapped reads in reference sample: %s" % 
              self.no_of_mapped_reads_control)
        print("Number of mapped reads in ChIP-Seq sample: %s" % 
              self.no_of_mapped_reads_chip)

    def calc_combined_factor(self):
        min_no_of_mapped_reads = float(min([
                    self.no_of_mapped_reads_chip, 
                    self.no_of_mapped_reads_control]))
        if self._factor != None:
            self._combine_factor = min_no_of_mapped_reads * self._factor
            print("Multiplication factor: %s (%s * %s)" % (
                    self._combine_factor, min_no_of_mapped_reads, self._factor))
        else:
            self._combine_factor = min_no_of_mapped_reads
            print("Multiplication factor: %s (min. number of aligned reads)" % (
                    self._combine_factor))

    def write_chip_and_control_wiggle_files(self):
        self._write_wiggle(
            self._calc_averaged_coverages(self.coverage_control), 
            "control",
            self._combine_factor/float(self.no_of_mapped_reads_control))
        self._write_wiggle(
            self._calc_averaged_coverages(self.coverage_chip), 
            "chip",
            self._combine_factor/float(self.no_of_mapped_reads_chip))

    def _calc_averaged_coverages(self, coverages):
        averaged_coverages = {}
        for element, element_coverages in coverages.items():
            averaged_coverages[element] = self._sliding_windows_average(
                element_coverages)
        return averaged_coverages

    def write_ratio_wiggle_file(self):
        self._write_wiggle(
            self.elements_and_coverage_ratios, "ratio", self._combine_factor)

    def _write_wiggle(self, elements_and_coverages, name, factor):
        output_fh = open("%s-%s.wig" % (self._output_prefix, name), "w")
        output_fh.write("track type=wiggle_0 name=\"ChipSeq %s\"\n" % (name))
        for element in sorted(elements_and_coverages.keys()):
            output_fh.write("variableStep chrom=%s span=1\n" % (element))
            # Remove position with as coverage of 0. pos is increased
            # by 1 as a translation from a 0-based sysem (Python list)
            # to a 1 based system (wiggle) takes place.
            output_fh.write(
                "\n".join(
                    ["%s %s" % (pos + 1, float(coverage) * factor)
                     for pos, coverage in
                     self._pos_and_coverages(elements_and_coverages[element])])
                + "\n")
        output_fh.close()

    def _pos_and_coverages(self, coverages):
        if self._keep_zero_coverage is False:
            for pos, coverage in filter(
                lambda pos_and_cov: pos_and_cov[1] != 0.0,
                enumerate(coverages)):
                yield (pos, coverage)
        else:
            for pos, coverage in enumerate(coverages):
                if pos % self._step_size == 2:
                    yield (pos, coverage)

    def compare(self):
        self.elements_and_coverage_ratios = {}
        for element, coverages in self.coverage_control.items():
            self.elements_and_coverage_ratios[
                element] = self._compare_coverages(element)

    def _prepare_coverage(self, bam_file):
        ref_seq_and_coverages = {}
        coverage_calculator = CoverageCalculator(
            read_count_splitting=False, uniqueley_aligned_only=False,
            first_base_only=False)
        for ref_seq, coverages in coverage_calculator.ref_seq_and_coverages(bam_file):
            assert len(coverages["forward"]) == len(coverages["reverse"])
            # Sum up the coverage of the forward and reverse strand
            summed_coverage = [
                abs(cov_for) + abs(cor_rev) for cov_for, cor_rev 
                in zip(coverages["forward"], coverages["reverse"])]
            ref_seq_and_coverages[ref_seq] = summed_coverage
        return ref_seq_and_coverages

    def _compare_coverages(self, element):
        cur_cov_control = self.coverage_control[element]
        cur_cov_chip = self.coverage_chip[element]
        if len(cur_cov_control) != len(cur_cov_chip):
            sys.stderr.write("Error! Different number of nucleotides.\n")
            sys.exit(2)
        if self._window_size != None:
            cur_cov_control = self._sliding_windows_average(cur_cov_control)
            cur_cov_chip = self._sliding_windows_average(cur_cov_chip)
        # Calculate the ratio of Chip data to control data
        coverage_ratios = [
            self._ratio(
                float(chip) / float(self.no_of_mapped_reads_chip),
                float(con) / float(self.no_of_mapped_reads_control))
                for chip, con in zip(cur_cov_chip, cur_cov_control)]
        return(coverage_ratios)

    def _print_file_names(self):
        print("Performing ChIP-Seq analysis")
        print("- Input file with ChIP-Seq data: %s" % self._bam_file_chip)
        print("- Input file with reference data: %s" % self._bam_file_control)
        print("- Output file prefix : %s" % self._output_prefix)

    def _sliding_windows_average(self, coverages):
        averaged_coverages = [0] * len(coverages)
        for pos in range(int(self._window_size/2),
                         len(coverages) - int(self._window_size/2),
                         self._step_size):
            averaged_coverages[pos] = (
                sum(coverages[pos-int(self._window_size/2):
                              pos+int(self._window_size/2)+1])
                / self._window_size)
        return(averaged_coverages)

    def _count_no_of_mapped_reads(self, bam_file):
        reads = {}

        with pysam.Samfile(bam_file, "rb") as bam_fh:
            for read in bam_fh.fetch():
                reads[read.qname] = 1
        return(len(reads))

    def _ratio(self, mult, div):
        try:
            return(mult/div)
        except:
            return(0.0)

if __name__ == "__main__":
   main()