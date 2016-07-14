import csv
from datetime import datetime
import sys
import math
import numpy as np
import pandas as pd
from covernantlib.wiggle import WiggleParser


def extract(args):
    coverage_extractor = CoverageExtractor(
        args.coordinate_file, args.coverage_file, args.output_prefix,
        args.flip_reverse_strand, args.matrix_alignment, args.window_size,
        args.step_size)
    log("Reading wiggle file")
    coverage_extractor.read_wiggle_file()
    log("Reading coordinate files")
    coverage_extractor.read_coordinate_file()
    log("Generate coverage matrix")
    coverage_extractor.extract_coverages()
    coverage_extractor.generate_coverage_matrix()
    log("Writing matrix")
    coverage_extractor.write_matrix_to_file()
    log("Calculating and writing averages")
    coverage_extractor.calc_averages()

    
def log(msg):
    sys.stderr.write("{} - {}\n".format(
        datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S'), msg))


class CoverageExtractor(object):
    """ """
    def __init__(self, coordinate_file, coverage_file, output_prefix,
                 flip_reverse_strand, matrix_alignment, window_size,
                 step_size):
        self._coordinate_file = coordinate_file
        self._coverage_file = coverage_file
        self._output_prefix = output_prefix
        self._flip_reverse_strand = flip_reverse_strand
        self._matrix_alignment = matrix_alignment
        self._replicons_and_coverages = {}
        self._coordinates = []
        self._coverage_lists = []
        self._window_size = window_size
        self._step_size = step_size

    def read_wiggle_file(self):
        wiggle_parser = WiggleParser()
        for entry in wiggle_parser.entries(open(self._coverage_file)):
            # Return dictionaries of the coordinate and absolute
            # coverage values
            self._replicons_and_coverages[entry.replicon] = dict(
                [(pos, abs(value)) for pos, value in entry.pos_value_pairs])

    def read_coordinate_file(self):
        for row in csv.reader(open(self._coordinate_file), delimiter="\t"):
            if row[0].startswith("#"):
                continue
            start, end = [int(pos) for pos in row[1:3]]
            assert start <= end
            strand = row[3]
            assert strand in ["+", "-"]
            self._coordinates.append(
                {"replicon": row[0], "start": start,
                 "end": end, "strand": strand})

    def extract_coverages(self):
        for coordinate in self._coordinates:
            self._check_coordinate_validity(coordinate)
            coverages = [
                self._replicons_and_coverages[coordinate["replicon"]].get(
                    pos, 0.0)
                for pos in range(coordinate["start"], coordinate["end"])
            ]
            if self._flip_reverse_strand and coordinate["strand"] == "-":
                coverages = coverages[::-1]
            self._coverage_lists.append(coverages)

    def _check_coordinate_validity(self, coordinate):
        pass
            
    def generate_coverage_matrix(self):
        max_range = max(
            [len(coverages) for coverages in self._coverage_lists])
        self._coverage_df = self._init_coverage_dataframe(max_range)
        for coordinate, coverages in zip(
                self._coordinates, self._coverage_lists):
            aligned_coverages = self._align_coverages(coverages, max_range)
            if self._window_size == 1 and self._step_size == 1:
                row = [coordinate[key] for key in
                       ["replicon", "start", "end",
                        "strand"]] + aligned_coverages
            else:
                row = [coordinate[key] for key in
                       ["replicon", "start", "end", "strand"]]
                for pos in range(int(self._window_size/2) + 1,
                                 max_range - int(self._window_size/2) + 1,
                                 self._step_size):
                    start = pos-int(self._window_size/2)
                    end = pos+int(self._window_size/2)
                    coverages = aligned_coverages[start:end]
                    row.append(np.mean(coverages))
            row_df = self._init_coverage_dataframe(max_range, [row])
            self._coverage_df = self._coverage_df.append(
                row_df, ignore_index=True)

    def _align_coverages(self, coverages, max_range):
        filling = np.nan
        aligned_coverages = list(np.full(max_range, filling))
        if self._matrix_alignment == "left":
            aligned_coverages[:len(coverages)] = coverages
        elif self._matrix_alignment == "center":
            margin = round((max_range - len(coverages)) / 2)
            aligned_coverages[margin:margin+len(coverages)] = coverages
        elif self._matrix_alignment == "right":
            left_margin = max_range - len(coverages)
            aligned_coverages[
                left_margin:left_margin+len(coverages)] = coverages
        return aligned_coverages
            
    def _init_coverage_dataframe(self, max_range, data=None):
        """
        - max_range
        - data - list of list containig data
        """
        return pd.DataFrame(
            data, columns=["Replicon", "Start", "End", "Strand"] +
            list(range(int(self._window_size/2) + 1,
                       max_range - int(self._window_size/2) + 1,
                       self._step_size)))

    def write_matrix_to_file(self):
        self._coverage_df.to_csv(
            "{}_matrix.csv".format(self._output_prefix),
            sep="\t", index=False, na_rep="NaN")

    def calc_averages(self):
        average_df = pd.DataFrame(columns=["Value"] +
                                  list(self._coverage_df.columns)[4:])
        average_df["Value"] = [
            "Median (all values)",
            "Mean (all values)",
            "Standard deviation (all values)",
            "Median (zeros removed)",
            "Mean (zeros removed)",
            "Standard deviation (zeros removed)",
            "Median (NAs removed)",
            "Mean (NAs removed)",
            "Standard deviation (NAs removed)",
            "Median (zeroes and NAs removed)",
            "Mean (zeroes and NAs removed)",
            "Standard deviation (zeroes and NAs removed)"
            ]
        for column_name in list(self._coverage_df.columns)[4:]:
            column = self._coverage_df[column_name]
            column_without_zeros = np.array(list(
                filter(lambda cov: cov != 0.0, column)))
            column_without_na = np.array(list(
                filter(lambda cov: not math.isnan(cov), column)))
            column_without_zeros_and_na = np.array(list(filter(
                lambda cov: not (math.isnan(cov) or cov == 0.0), column)))
            average_df[column_name] = [
                column.median(),
                column.mean(),
                column.std(),
                self._save_column_operation(np.median, column_without_zeros),
                self._save_column_operation(np.mean, column_without_zeros),
                self._save_column_operation(np.std, column_without_zeros),
                self._save_column_operation(np.median, column_without_na),
                self._save_column_operation(np.mean, column_without_na),
                self._save_column_operation(np.std, column_without_na),
                self._save_column_operation(
                    np.median, column_without_zeros_and_na),
                self._save_column_operation(
                    np.mean, column_without_zeros_and_na),
                self._save_column_operation(
                    np.std, column_without_zeros_and_na)
            ]
        average_df.to_csv(
            "{}_averages.csv".format(self._output_prefix),
            sep="\t", index=False)

    def _save_column_operation(self, function, column):
        if len(column) > 0:
            return function(column)
        else:
            return "-"
