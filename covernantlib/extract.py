import csv
from datetime import datetime
import sys
import numpy as np
from covernantlib.wiggle import WiggleParser


def extract(args):
    coverage_extractor = CoverageExtractor(
        args.position_file, args.coverage_file, args.output_prefix,
        args.flip_reverse_strand)
    log("Reading wiggle file")
    coverage_extractor.read_wiggle_file()
    log("Reading position files")
    coverage_extractor.read_position_file()
    log("Calculating averages")
    coverage_extractor.calc_averages()


def log(msg):
    sys.stderr.write("{} - {}\n".format(
        datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S'), msg))


class CoverageExtractor(object):
    """ """
    def __init__(self, position_file, coverage_file, output_prefix,
                 flip_reverse_strand):
        self._position_file = position_file
        self._coverage_file = coverage_file
        self._output_prefix = output_prefix
        self._flip_reverse_strand = flip_reverse_strand
        self._replicons_and_coverages = {}
        self._all_coverage_lists = []

    def read_wiggle_file(self):
        wiggle_parser = WiggleParser()
        for entry in wiggle_parser.entries(open(self._coverage_file)):
            # Return dictionaries of the position and absolute coverage values
            self._replicons_and_coverages[entry.replicon] = dict(
                [(pos, abs(value)) for pos, value in entry.pos_value_pairs])

    def read_position_file(self):
        for row in csv.reader(open(self._position_file), delimiter="\t"):
            if row[0].startswith("#"):
                continue
            pos_1, pos_2 = [int(pos) for pos in row[1:3]]
            assert pos_1 <= pos_2
            strand = row[3]
            assert strand in ["+", "-"]
            pos_range = round((pos_2 - pos_1) / 2)
            positions_and_coverage = self._replicons_and_coverages[row[0]]
            mid_pos = pos_1 + pos_range
            positions = range(mid_pos - pos_range, mid_pos + pos_range)
            coverages = [positions_and_coverage.get(
                pos, 0.0) for pos in positions]
            if self._flip_reverse_strand and strand == "-":
                coverages = coverages[::-1]
            self._all_coverage_lists.append(coverages)

    def calc_averages(self):
        max_list_length = max(
            [len(coverages) for coverages in self._all_coverage_lists])
        matrix = []
        for coverage_list in self._all_coverage_lists:
            no_of_gaps_to_add_to_sides = round((
                max_list_length - len(coverage_list))/2)
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
