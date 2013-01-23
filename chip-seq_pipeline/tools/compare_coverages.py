#!/usr/bin/env python

__description__ = ""
__author__ = "Konrad Foerstner <konrad@foerstner.org>"
__copyright__ = "2012 by Konrad Foerstner <konrad@foerstner.org>"
__license__ = "ISC license"
__email__ = "konrad@foerstner.org"
__version__ = ""

import argparse
import sys
sys.path.append("..")
from libs.coveragecreator import CoverageCreator
from libs.sam import SamParser

def main():
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument("bam_file_control")
    parser.add_argument("bam_file_chip")
    parser.add_argument("--output", "-o", dest="output_file",
                        default="chip_seq.wig", required=False)
    parser.add_argument("--window_size", type=int, default=None,
                        help="Window size for sliding window average calculation.")
    parser.add_argument("--step_size", type=int, default=1,
                        help="Step size for sliding window average calculation."
                        " Default is 1.")
    parser.add_argument("--factor", type=float, default=None,
                        help="A factor the final ratio is multiplied with.")
    # TODO
    # parser.add_argument("--pseudocount",
    #                    help="Use pseudocounts to initate the coverage vector.")
    # parser.add_argument("--uniquely", help="Use uniquely mapped reads only")
    # parser.add_argument("--zero_div", help="What to do if the divisor is zero - 'zero', "minus", "addone"")
    # 'zero' - the value becomes zero
    # 'minus' - the value becomes
    args = parser.parse_args()
    coverage_comparer = CoverageComparer(
        args.bam_file_control, args.bam_file_chip, args.output_file,
        args.window_size, args.step_size, args.factor)
    coverage_comparer.compare()
    coverage_comparer.write_to_wiggle_files()

class CoverageComparer(object):

    def __init__(self, bam_file_control, bam_file_chip, output_file,
                 window_size, step_size, factor):
        self._bam_file_control = bam_file_control
        self._bam_file_chip = bam_file_chip
        self._output_file = output_file
        if window_size % 2 == 0:
            sys.stderr.write("Error. Window size must be an odd number!\n")
            sys.exit(2)
        self._window_size = window_size
        self._step_size = step_size
        self._factor = factor

    def _calc_coverage(self, bam_file):
        coverage_creator = CoverageCreator()
        coverage_creator.init_coverage_lists(bam_file)
        coverage_creator.count_coverage(bam_file)

    def compare(self):
        self._print_file_names()
        no_of_mapped_reads_chip = self._count_no_of_mapped_reads(
            self._bam_file_chip)
        no_of_mapped_reads_control = self._count_no_of_mapped_reads(
            self._bam_file_control)
        coverage_control = self._calc_coverage(self._bam_file_control)
        coverage_chip = self._calc_coverage(self._bam_file_chip)
        self.elements_and_coverage_ratios = {}
        for element in coverage_control.elements_and_coverages.keys():
            cur_cov_control = coverage_control.elements_and_coverages[
                element]
            cur_cov_chip = coverage_chip.elements_and_coverages[element]
            if len(cur_cov_control) != len(cur_cov_chip):
                sys.stderr.write("Error! Different number of nucleotides.\n")
                sys.exit(2)
            if self._window_size != None:
                cur_cov_control = self._sliding_windows_average(
                    cur_cov_control)
                cur_cov_chip = self._sliding_windows_average(
                    cur_cov_chip)
            # Calculate the ratio of Chip data to control data
            self.elements_and_coverage_ratios[element] = [
                self._ratio(
                    float(chip) / float(no_of_mapped_reads_chip),
                    float(con) / float(no_of_mapped_reads_control))
                    for chip, con in zip(cur_cov_chip, cur_cov_control)]
            # Multiply the ratio by a factor if given
            if self._factor != None:
                self.elements_and_coverage_ratios[element] = [
                    ratio * self._factor
                    for ratio in self.elements_and_coverage_ratios[element]]

    def _print_file_names(self):
        print("Performing ChIP-Seq analysis")
        print("- Input file with ChIP-Seq data: %s" % self._bam_file_chip)
        print("- Input file with reference data: %s" % self._bam_file_control)
        print("- Output file: %s" % self._output_file)

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
        sam_parser = SamParser()
        for entry in sam_parser.entries_bam(bam_file):
            reads[entry.query_id] = 1
        return(len(reads))

    def _ratio(self, mult, div):
        try:
            return(mult/div)
        except:
            return(0.0)

    def write_to_wiggle_files(self):
        output_fh = open("%s" % (self._output_file), "w")
        output_fh.write("track type=wiggle_0 name=\"%s\"\n" % (
            "ChipSeq")) # TODO: give better name
        for element in sorted(self.elements_and_coverage_ratios.keys()):
            output_fh.write("variableStep chrom=%s span=1\n" % (element))
            # Filter values of 0 and multiply other the remaining
            # ones by the given factor. pos is increased by 1 as a
            # translation from a 0-based sysem (Python list) to a
            # 1 based system (wiggle) takes place.
            output_fh.write(
                "\n".join(
                    ["%s %s" % (pos + 1, coverage)
                     for pos, coverage in
                     filter(lambda pos_and_cov: pos_and_cov[1] != 0.0,
                            enumerate(self.elements_and_coverage_ratios[element]))]) + "\n")
        output_fh.close()

if __name__ == "__main__":
   main()
