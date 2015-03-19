#!/usr/bin/env python
"""
FUNCTION: 

USAGE: 

Copyright (c) 2015, Konrad Foerstner <konrad@foerstner.org>

Permission to use, copy, modify, and/or distribute this software for
any purpose with or without fee is hereby granted, provided that the
above copyright notice and this permission notice appear in all
copies.

THE SOFTWARE IS PROVIDED 'AS IS' AND THE AUTHOR DISCLAIMS ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE
AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL
DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER
TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.

"""
__description__ = ""
__author__ = "Konrad Foerstner <konrad@foerstner.org>"
__copyright__ = "2015 by Konrad Foerstner <konrad@foerstner.org>"
__license__ = "ISC license"
__email__ = "konrad@foerstner.org"
__version__ = ""

import argparse
import csv
from wiggle import WiggleParser
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument("position_file")
    parser.add_argument("coverage_file")
    parser.add_argument("--output_prefix", default="output")
    parser.add_argument("--pos_range", type=int)

    args = parser.parse_args()

    coverage_combiner = CoverageCombiner(
        args.position_file, args.coverage_file, args.pos_range, 
        args.output_prefix)
    coverage_combiner.read_wiggle_file()
    coverage_combiner.read_position_file()
    coverage_combiner.calc_averages()

    
class CoverageCombiner(object):
        
    def __init__(self, position_file, coverage_file, pos_range, output_prefix):
        self._position_file = position_file
        self._coverage_file = coverage_file
        self._pos_range = pos_range
        self._output_prefix = output_prefix
        self._replicons_and_coverages = {}
            
    def read_wiggle_file(self):
        wiggle_parser = WiggleParser()
        for entry in wiggle_parser.entries(open(self._coverage_file)):
            # Return dictionaries of the position and absolute coverage values
            self._replicons_and_coverages[entry.replicon] = dict(
                [(pos, abs(value)) for pos, value in entry.pos_value_pairs])
                
    def read_position_file(self):
        self._all_coverage_lists = []
        for row in csv.reader(open(self._position_file), delimiter="\t"):
            if row[0].startswith("#"):
                continue
            pos_1, pos_2 = [int(pos) for pos in row[1:3]]
            assert(pos_1 <= pos_2)
            pos_range = (pos_2 - pos_1) / 2
            positions_and_coverage = self._replicons_and_coverages[row[0]]
            mid_pos = pos_1 + pos_range
            positions = range(mid_pos - pos_range, mid_pos + pos_range)
            coverages = [positions_and_coverage.get(
                pos, 0.0) for pos in positions]
            self._all_coverage_lists.append(coverages)
        
    def calc_averages(self):
        max_list_length = max(
            [len(coverages) for coverages in self._all_coverage_lists])
        matrix = []
        for coverage_list in self._all_coverage_lists:
            no_of_gaps_to_add_to_sides = (
                max_list_length - len(coverage_list))/2
            coverage_list = [
                0.0] * no_of_gaps_to_add_to_sides + coverage_list + [
                    0.0] * no_of_gaps_to_add_to_sides
            matrix.append(coverage_list)
        matrix = np.array(matrix)
        output_fh = open("{}_combined.csv".format(self._output_prefix), "w")
        output_fh.write("\t".join([
            "Pos", "Mean", "Median", "Standard deviation"]) + "\n")
        for pos in range(max_list_length):
            column = matrix[:, pos]
            shifted_pos = pos - (max_list_length/2)
            if np.sum(column) == 0.0:
                continue
            output_fh.write("\t".join([str(cell) for cell in [
                shifted_pos,
                np.mean(column),
                np.median(column),
                np.std(column)]]) + "\n")
        output_fh.close()

if __name__ == "__main__":
    main()
