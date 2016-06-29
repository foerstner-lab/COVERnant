import pybedtools
from collections import defaultdict
import numpy as np


def bed_to_wiggle(input_file, output_prefix, window_size, step_size):
    for strand in ["+", "-"]:
        _generate_coverage_for_strand(
            strand, input_file, output_prefix, window_size, step_size)


def _generate_coverage_for_strand(strand, input_file, output_prefix,
                                  window_size, step_size):
    replicons_and_coverages = init_coverage_dict(input_file)
    entries = pybedtools.BedTool(input_file)
    for replicon, coverages in replicons_and_coverages.items():
        for window_center in range(
                int(window_size/2),
                len(coverages)-int(window_size/2), step_size):
            interval = pybedtools.create_interval_from_list(
                [replicon, window_center-int(window_size/2),
                 window_center+int(window_size/2), "tmp_window", ".", strand])
            factor = 1
            if strand == "-":
                factor = -1
            replicons_and_coverages[replicon][
                window_center] = factor * entries.count_hits(
                    interval, same_strand=True)

    output_fh = open("{}_{}_strand.wig".format(
        output_prefix, strand), "w")
    output_fh.write("track type=wiggle_0 name=\"coverage {} strand\"\n".format(
        strand))
    
    for replicon in sorted(replicons_and_coverages.keys()):
        output_fh.write("variableStep chrom=%s span=1\n" % (replicon))
        # Remove position with as coverage of 0. pos is increased
        # by 1 as a translation from a 0-based sysem (Python list)
        # to a 1 based system (wiggle) takes place.
        output_fh.write(
            "\n".join(
                ["%s %s" % (pos + 1, coverage)
                 for pos, coverage in enumerate(
                         replicons_and_coverages[replicon])]) + "\n")
    output_fh.close()


def init_coverage_dict(input_file):
    replicons_and_max_end_pos = defaultdict(int)
    for entry in pybedtools.BedTool(input_file):
        if entry.end > replicons_and_max_end_pos[entry.chrom]:
            replicons_and_max_end_pos[entry.chrom] = entry.end
    replicons_and_coverages = {}
    for replicon, end_pos in replicons_and_max_end_pos.items():
        replicons_and_coverages[replicon] = np.zeros(end_pos)
    return replicons_and_coverages
