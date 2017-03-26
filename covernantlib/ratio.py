import pysam
import sys
import numpy as np
from covernantlib.coveragecalculator import CoverageCalculator


def calc_ratio(args):
    coverage_ratio_calculator = CoverageRatioCalculator(args)
    coverage_ratio_calculator.set_names()
    coverage_ratio_calculator.calc_coverages()
    coverage_ratio_calculator.print_no_aligned_bases()
    coverage_ratio_calculator.calc_normalization_factors()
    coverage_ratio_calculator.normalize_coverages()
    coverage_ratio_calculator.write_numerator_and_denominator_wiggle_files()
    coverage_ratio_calculator.calc_coverage_ratio()
    coverage_ratio_calculator.write_ratio_wiggle_files()


class CoverageRatioCalculator(object):
    """Calculate the coverages and ratio of two input BAM files

    Input: Two BAM files
       - One later used for ratio calculation as numerator
       - One later used for ratio calculation as denominator

    Output - 9 wiggle files:
       - Three wiggle file sets:
         - Numerator
           - One wiggle file forward strand (CPB)
           - One wiggle file reverse strand (CPB)
           - One wiggle file both strands combined (CPB)
         - Denominator
           - One wiggle file forward strand (CPB)
           - One wiggle file reverse strand (CPB)
           - One wiggle file both strands combined (CPB)
         - Ratio:
           - Ratio of the normalized numerator coverages and the
             normalized denominator coverages forward strand
           - Ratio of the normalized numerator coverages and the
             normalized denominator coverages reverse strand
           - Ratio of the normalized numerator coverages and the
             normalized denominator coverages both strand combined

    Normalization:

       - Defaul: Count the number of valid alignments (one read can
         lead to several alignments depending on the mapper/downstream
         processing) and use this for the normalization. Counts per
         billion - "counts" are read for single end libraries and
         fragments for paired end libraries. Even for the strand
         specific coverages the total number of alignments are used.

       - Alternatively the user can specificy factors that will be
         used instead of the calculated number of alignments. This
         might be useful e.g. if spike-ins are used.

    Main parameters:
       - Window size and step size for the calculation of the coverage
         averages inside of a sliding window.
    """

    def __init__(self, args):
        self._denominator_bam_file = args.denominator_bam_file
        self._numerator_bam_file = args.numerator_bam_file
        self._output_prefix = args.output_prefix
        if args.window_size is not None and args.window_size % 2 == 0:
            sys.stderr.write("Error. Window size must be an odd number!\n")
            sys.exit(2)
        self._window_size = args.window_size
        self._step_size = args.step_size
        self._numerator_factor_given = args.factor_numerator
        self._denominator_factor_given = args.factor_denominator
        self._keep_zero_coverage = args.keep_zero_coverage
        self._denominator_name = args.denominator_name
        self._numerator_name = args.numerator_name
        self._ratio_name = args.ratio_name
        self._paired_end = args.paired_end

    def calc_coverages(self):
        (self.coverage_denominator,
         self.coverage_denominator_forward,
         self.coverage_denominator_reverse,
         self.no_of_mapped_bases_denominator) = self._prepare_coverage(
             self._denominator_bam_file)
        (self.coverage_numerator,
         self.coverage_numerator_forward,
         self.coverage_numerator_reverse,
         self.no_of_mapped_bases_numerator) = self._prepare_coverage(
             self._numerator_bam_file)

    def print_no_aligned_bases(self):
        print("Number of mapped bases in numerator sample total: %s" %
              self.no_of_mapped_bases_numerator)
        print("Number of used bases in denominator sample - "
              "total: %s" % self.no_of_mapped_bases_denominator)        

    def calc_normalization_factors(self):
        """Calculate the normalization factor based on the number of aligned
        nucleotides of the two input libraries and generate counts per
        billion (i.e. read counts for single end; fragment counts for
        paired end).
        """
        # Default behavior if no factor is given - use the number of
        # alignments for the normalization
        if not (self._numerator_factor_given is not None
                or self._denominator_factor_given is not None):
            self._numerator_normalization_factor = 1000000000.0/float(
                self.no_of_mapped_bases_numerator)
            self._denominator_normalization_factor = 1000000000.0/float(
                self.no_of_mapped_bases_denominator)
            print("CPB factor numerator: %s (1B/%s)" % (
                self._numerator_normalization_factor,
                self.no_of_mapped_bases_numerator))
            print("CPB factor denominator: %s (1B/%s)" % (
                self._denominator_normalization_factor,
                self.no_of_mapped_bases_denominator))
        # In case one or both factors are set:
        else:
            if self._numerator_factor_given is None:
                self._numerator_factor_given = 1.0
            if self._denominator_factor_given is None:
                self._denominator_factor_given = 1.0
            self._numerator_normalization_factor = self._numerator_factor_given
            self._denominator_normalization_factor = (
                self._denominator_factor_given)
            print("Factor numerator: %s (manually set)" % (
                self._numerator_normalization_factor))
            print("Factor denominator: %s (manually set)" % (
                self._denominator_normalization_factor))

    def normalize_coverages(self):
        self.coverage_denominator_normalized = (
            self._normalize_coverages(
                self.coverage_denominator,
                self._denominator_normalization_factor))
        self.coverage_denominator_forward_normalized = (
            self._normalize_coverages(
                self.coverage_denominator_forward,
                self._denominator_normalization_factor))
        self.coverage_denominator_reverse_normalized = (
            self._normalize_coverages(
                self.coverage_denominator_reverse,
                self._denominator_normalization_factor))
        self.coverage_numerator_normalized = (
            self._normalize_coverages(
                self.coverage_numerator,
                self._numerator_normalization_factor))
        self.coverage_numerator_forward_normalized = (
            self._normalize_coverages(
                self.coverage_numerator_forward,
                self._numerator_normalization_factor))
        self.coverage_numerator_reverse_normalized = (
            self._normalize_coverages(
                self.coverage_numerator_reverse,
                self._numerator_normalization_factor))

    def _normalize_coverages(self, coverages, factor):
        normalized_coverages = {}
        for element, element_coverages in coverages.items():
            normalized_coverages[element] = element_coverages * factor
        return normalized_coverages
            
    def set_names(self):
        if self._denominator_name is None:
            self._denominator_name = self._file_name_to_name(
                self._denominator_bam_file)
        if self._numerator_name is None:
            self._numerator_name = self._file_name_to_name(
                self._numerator_bam_file)
        if self._ratio_name is None:
            self._ratio_name = "{}_vs_{}".format(
                self._file_name_to_name(self._denominator_bam_file),
                self._file_name_to_name(self._numerator_bam_file))
        print(self._denominator_name, self._numerator_name, self._ratio_name)

    def _file_name_to_name(self, file_name):
        name = file_name.split("/")[-1].replace(".bam", "")
        return name

    def write_numerator_and_denominator_wiggle_files(self):
        """Write 6 wiggles files. 3 for denominator and for nominator - total
        coverage, coverage forward strand and coverage reverse
        strand. At this point the normalization by count per billion
        takes place.
        """
        self._write_wiggle(
            self._calc_averaged_coverages(
                self.coverage_denominator_normalized),
            "denominator_{}".format(self._denominator_name))
        self._write_wiggle(
            self._calc_averaged_coverages(
                self.coverage_denominator_forward_normalized),
            "denominator_{}_forward".format(self._denominator_name))
        self._write_wiggle(
            self._calc_averaged_coverages(
                self.coverage_denominator_reverse_normalized),
            "denominator_{}_reverse".format(self._denominator_name))
        self._write_wiggle(
            self._calc_averaged_coverages(self.coverage_numerator_normalized),
            "numerator_{}".format(self._numerator_name))
        self._write_wiggle(
            self._calc_averaged_coverages(
                self.coverage_numerator_forward_normalized),
            "numerator_{}_forward".format(self._numerator_name))
        self._write_wiggle(
            self._calc_averaged_coverages(
                self.coverage_numerator_reverse_normalized),
            "numerator_{}_reverse".format(self._numerator_name))

    def _calc_averaged_coverages(self, coverages):
        averaged_coverages = {}
        if self._window_size is None:
            self._window_size = 1
        for element, element_coverages in coverages.items():
            averaged_coverages[element] = self._sliding_windows_average(
                element_coverages)
        return averaged_coverages

    def write_ratio_wiggle_files(self):
        self._write_wiggle(
            self.elements_and_coverage_ratios,
            "ratio_{}".format(self._ratio_name))
        self._write_wiggle(
            self.elements_and_coverage_ratios_forward,
            "ratio_{}_forward".format(self._ratio_name))
        self._write_wiggle(
            self.elements_and_coverage_ratios_reverse,
            "ratio_{}_reverse".format(self._ratio_name))

    def _write_wiggle(self, elements_and_coverages, name):
        output_fh = open("%s%s.wig" % (self._output_prefix, name), "w")
        output_fh.write("track type=wiggle_0 name=\"%s%s\"\n" % (
            self._output_prefix.split("/")[-1], name))
        for element in sorted(elements_and_coverages.keys()):
            output_fh.write("variableStep chrom=%s span=1\n" % (element))
            # Remove position with as coverage of 0. pos is increased
            # by 1 as a translation from a 0-based sysem (Python list)
            # to a 1 based system (wiggle) takes place.
            output_fh.write(
                "\n".join(
                    ["%s %s" % (pos + 1, float(coverage))
                     for pos, coverage in
                     self._pos_and_coverages(elements_and_coverages[element])])
                + "\n")
        output_fh.close()

    def _pos_and_coverages(self, coverages):
        if not self._keep_zero_coverage:
            # As the values between the nucleotides selected by the
            # step size are zero only the selected nucleotides will be
            # returned due to this filter.
            for pos, coverage in filter(
                lambda pos_and_cov: pos_and_cov[1] != 0.0,
                    enumerate(coverages)):
                yield (pos, coverage)
        else:
            # Here the step size has to be used explicitly to return
            # the selected nucleotide
            for pos in range(int(self._window_size/2),
                             len(coverages) - int(self._window_size/2),
                             self._step_size):
                yield (pos, coverages[pos])

    def calc_coverage_ratio(self):
        self.elements_and_coverage_ratios = {}
        for element in self.coverage_denominator_normalized.keys():
            self.elements_and_coverage_ratios[
                element] = self._compare_coverages(
                    element, self.coverage_denominator_normalized,
                    self.coverage_numerator_normalized)
        self.elements_and_coverage_ratios_forward = {}
        for element in self.coverage_denominator_forward_normalized.keys():
            self.elements_and_coverage_ratios_forward[
                element] = self._compare_coverages(
                    element, self.coverage_denominator_forward_normalized,
                self.coverage_numerator_forward_normalized)
        self.elements_and_coverage_ratios_reverse = {}
        for element in self.coverage_denominator_reverse_normalized.keys():
            self.elements_and_coverage_ratios_reverse[
                element] = self._compare_coverages(
                    element, self.coverage_denominator_reverse_normalized,
                    self.coverage_numerator_reverse_normalized)

    def _compare_coverages(self, element, coverage_denominator,
                           coverage_numerator):
        cur_cov_denominator = coverage_denominator[element]
        cur_cov_numerator = coverage_numerator[element]
        if len(cur_cov_denominator) != len(cur_cov_numerator):
            sys.stderr.write("Error! Different number of nucleotides.\n")
            sys.exit(2)
        if self._window_size is not None:
            cur_cov_denominator = self._sliding_windows_average(
                cur_cov_denominator)
            cur_cov_numerator = self._sliding_windows_average(
                cur_cov_numerator)
        # Calculate the ratio of numerator data to denominator data
        coverage_ratios = [
            self._ratio(
                float(numerator), float(denominator))
            for numerator, denominator in zip(
                    cur_cov_numerator, cur_cov_denominator)]
        return(coverage_ratios)

    def _prepare_coverage(self, bam_file):
        """We generate three sets of coverages. One for the forward strand
        only, one for the reverse strand only, and one for the sum of
        both. The coverage returned are raw i.e. are not normalized to
        the number of used alignment.
        """
        ref_seq_and_coverages_sum = {}
        ref_seq_and_coverages_forward = {}
        ref_seq_and_coverages_reverse = {}
        coverage_calculator = CoverageCalculator(paired_end=self._paired_end)
        for ref_seq, coverages in coverage_calculator.ref_seq_and_coverages(
                bam_file):
            assert len(coverages["forward"]) == len(coverages["reverse"])
            # Sum up the coverage of the forward and reverse strand
            summed_coverage = coverages["forward"] + coverages["reverse"]
            ref_seq_and_coverages_sum[ref_seq] = summed_coverage
            ref_seq_and_coverages_forward[ref_seq] = coverages["forward"]
            ref_seq_and_coverages_reverse[ref_seq] = coverages["reverse"]
        return (ref_seq_and_coverages_sum, ref_seq_and_coverages_forward,
                ref_seq_and_coverages_reverse,
                coverage_calculator.no_of_used_bases)

    def _print_file_names(self):
        print("Performing ChIP-Seq analysis")
        print("- Input file with numerator data: %s" %
              self._numerator_bam_file)
        print("- Input file with denominator data: %s" %
              self._denominator_bam_file)
        print("- Output file prefix: %s" % self._output_prefix)

    def _sliding_windows_average(self, coverages):
        averaged_coverages = np.zeros(len(coverages))
        for pos in range(int(self._window_size/2),
                         len(coverages) - int(self._window_size/2),
                         self._step_size):
            averaged_coverages[pos] = (
                sum(coverages[pos-int(self._window_size/2):
                              pos+int(self._window_size/2)+1])
                / self._window_size)
        return(averaged_coverages)

    def _ratio(self, mult, div):
        try:
            return(mult/div)
        except:
            return(0.0)
