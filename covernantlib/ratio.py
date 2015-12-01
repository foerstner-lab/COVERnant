import pysam
import sys
from covernantlib.coveragecalculator import CoverageCalculator


def calc_ratio(args):
    coverage_ratio_calculator = CoverageRatioCalculator(args)
    coverage_ratio_calculator.set_names()
    coverage_ratio_calculator.calc_coverages()
    coverage_ratio_calculator.print_no_aligned_reads()
    coverage_ratio_calculator.calc_combined_factor()
    coverage_ratio_calculator.write_numerator_and_denominator_wiggle_files()
    coverage_ratio_calculator.compare()
    coverage_ratio_calculator.write_ratio_wiggle_files()


class CoverageRatioCalculator(object):

    def __init__(self, args):
        self._denominator_bam_file = args.denominator_bam_file
        self._numerator_bam_file = args.numerator_bam_file
        self._output_prefix = args.output_prefix
        if args.window_size is not None and args.window_size % 2 == 0:
            sys.stderr.write("Error. Window size must be an odd number!\n")
            sys.exit(2)
        self._window_size = args.window_size
        self._step_size = args.step_size
        self._factor = args.factor
        self._keep_zero_coverage = args.keep_zero_coverage
        self._denominator_name = args.denominator_name
        self._numerator_name = args.numerator_name
        self._ratio_name = args.ratio_name
        self._paired_end = args.paired_end

    def calc_coverages(self):
        self._print_file_names()
        (self.no_of_mapped_reads_numerator,
        self.no_of_mapped_reads_numerator_forward,
         self.no_of_mapped_reads_numerator_reverse
        ) = self._count_no_of_mapped_reads(
            self._numerator_bam_file)
        (self.no_of_mapped_reads_denominator,
         self.no_of_mapped_reads_denominator_forward,
         self.no_of_mapped_reads_denominator_reverse
        ) = self._count_no_of_mapped_reads(
            self._denominator_bam_file)
        (self.coverage_denominator,
         self.coverage_denominator_forward,
         self.coverage_denominator_reverse) = self._prepare_coverage(
             self._denominator_bam_file)
        (self.coverage_numerator,
         self.coverage_numerator_forward,
         self.coverage_numerator_reverse) = self._prepare_coverage(
             self._numerator_bam_file)

    def print_no_aligned_reads(self):
        print("Number of mapped reads in denominator sample - "
              "total: %s" % self.no_of_mapped_reads_denominator)
        print("Number of mapped reads in denominator sample - "
              "forward strand: %s" %
              self.no_of_mapped_reads_denominator_forward)
        print("Number of mapped reads in denominator sample - "
              "reverse strand: %s" %
              self.no_of_mapped_reads_denominator_reverse)
        print("Number of mapped reads in numerator sample total: %s" %
              self.no_of_mapped_reads_numerator)
        print("Number of mapped reads in numerator sample - "
              "forward strand: %s" % self.no_of_mapped_reads_numerator_forward)
        print("Number of mapped reads in numerator sample - "
              "reverse strand: %s" % self.no_of_mapped_reads_numerator_reverse)

    def calc_combined_factor(self):
        self._numerator_rpm_factor = 1000000.0/float(
            self.no_of_mapped_reads_numerator)
        self._denominator_rpm_factor = 1000000.0/float(
            self.no_of_mapped_reads_denominator)
        self._ratio_factor = float(self.no_of_mapped_reads_numerator)/float(
            self.no_of_mapped_reads_denominator)
        if self._factor is not None:
            self._ratio_factor *= self._factor
            print("Ratio factor: %s (%s/%s * %s) " % (
                self._ratio_factor, self.no_of_mapped_reads_numerator,
                self.no_of_mapped_reads_denominator, self._factor))
        else:
            print("Ratio factor: %s (%s/%s) " % (
                self._ratio_factor, self.no_of_mapped_reads_numerator,
                self.no_of_mapped_reads_denominator,))
        print("RPM factor numerator: %s (1M/%s)" % (
            self._numerator_rpm_factor, self.no_of_mapped_reads_numerator))
        print("RPM factor denominator: %s (1M/%s)" % (
            self._denominator_rpm_factor, self.no_of_mapped_reads_denominator))

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
        self._write_wiggle(
            self._calc_averaged_coverages(self.coverage_denominator),
            "denominator_{}".format(self._denominator_name),
            self._denominator_rpm_factor)
        self._write_wiggle(
            self._calc_averaged_coverages(self.coverage_denominator_forward),
            "denominator_{}_forward".format(self._denominator_name),
            self._denominator_rpm_factor)
        self._write_wiggle(
            self._calc_averaged_coverages(self.coverage_numerator_forward),
            "numerator_{}_forward".format(self._numerator_name),
            self._numerator_rpm_factor)
        self._write_wiggle(
            self._calc_averaged_coverages(self.coverage_numerator),
            "numerator_{}".format(self._numerator_name),
            self._numerator_rpm_factor)
        self._write_wiggle(
            self._calc_averaged_coverages(self.coverage_denominator_reverse),
            "denominator_{}_reverse".format(self._denominator_name),
            self._denominator_rpm_factor)
        self._write_wiggle(
            self._calc_averaged_coverages(self.coverage_numerator_reverse),
            "numerator_reverse".format(self._numerator_name),
            self._numerator_rpm_factor)

    def _calc_averaged_coverages(self, coverages):
        averaged_coverages = {}
        # TODO - CHECK
        if self._window_size is None:
            self._window_size = 1
        for element, element_coverages in coverages.items():
            averaged_coverages[element] = self._sliding_windows_average(
                element_coverages)
        return averaged_coverages

    def write_ratio_wiggle_files(self):
        self._write_wiggle(
            self.elements_and_coverage_ratios,
            "ratio_{}".format(self._ratio_name), self._ratio_factor)
        self._write_wiggle(
            self.elements_and_coverage_ratios_forward,
            "ratio_{}_forward".format(self._ratio_name),
            self._ratio_factor)
        self._write_wiggle(
            self.elements_and_coverage_ratios_reverse,
            "ratio_{}_reverse".format(self._ratio_name),
            self._ratio_factor)

    def _write_wiggle(self, elements_and_coverages, name, factor):
        output_fh = open("%s-%s.wig" % (self._output_prefix, name), "w")
        output_fh.write("track type=wiggle_0 name=\"%s_%s\"\n" % (
            self._output_prefix.split("/")[-1], name))
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

    def compare(self):
        self.elements_and_coverage_ratios = {}
        for element, coverages in self.coverage_denominator.items():
            self.elements_and_coverage_ratios[
                element] = self._compare_coverages(
                    element, self.coverage_denominator,
                    self.coverage_numerator)
        self.elements_and_coverage_ratios_forward = {}
        for element, coverages in self.coverage_denominator_forward.items():
            self.elements_and_coverage_ratios_forward[
                element] = self._compare_coverages(
                    element, self.coverage_denominator_forward,
                    self.coverage_numerator_forward)
        self.elements_and_coverage_ratios_reverse = {}
        for element, coverages in self.coverage_denominator_reverse.items():
            self.elements_and_coverage_ratios_reverse[
                element] = self._compare_coverages(
                    element, self.coverage_denominator_reverse,
                    self.coverage_numerator_reverse)

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
                float(numerator) / float(self.no_of_mapped_reads_numerator),
                float(con) / float(self.no_of_mapped_reads_denominator))
            for numerator, con in zip(cur_cov_numerator, cur_cov_denominator)]
        return(coverage_ratios)

    def _prepare_coverage(self, bam_file):
        ref_seq_and_coverages_sum = {}
        ref_seq_and_coverages_forward = {}
        ref_seq_and_coverages_reverse = {}
        coverage_calculator = CoverageCalculator(
            read_count_splitting=False, uniqueley_aligned_only=False,
            first_base_only=False, paired_end=self._paired_end)
        for ref_seq, coverages in coverage_calculator.ref_seq_and_coverages(
                bam_file):
            assert len(coverages["forward"]) == len(coverages["reverse"])
            # Sum up the coverage of the forward and reverse strand
            summed_coverage = [
                abs(cov_for) + abs(cor_rev) for cov_for, cor_rev
                in zip(coverages["forward"], coverages["reverse"])]
            ref_seq_and_coverages_sum[ref_seq] = summed_coverage
            ref_seq_and_coverages_forward[ref_seq] = [
                abs(cor_forw) for cor_forw in coverages["forward"]]
            ref_seq_and_coverages_reverse[ref_seq] = [
                abs(cor_rev) for cor_rev in coverages["reverse"]]
        print("Number of used alignments: {}".format(
            coverage_calculator.used_alignmets))
        # print("Number of discarded alignments: {}".format(
        #     coverage_calculator.discared_alignments))
        return (ref_seq_and_coverages_sum, ref_seq_and_coverages_forward,
                ref_seq_and_coverages_reverse)

    def _print_file_names(self):
        print("Performing ChIP-Seq analysis")
        print("- Input file with numerator data: %s" %
              self._numerator_bam_file)
        print("- Input file with denominator data: %s" %
              self._denominator_bam_file)
        print("- Output file prefix: %s" % self._output_prefix)

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
        reads_forward = set()
        reads_reverse = set()
        with pysam.Samfile(bam_file, "rb") as bam_fh:
            for read in bam_fh.fetch():
                if read.is_reverse is False:
                    reads_forward.add(read.qname)
                else:
                    reads_reverse.add(read.qname)
        return(len(reads_forward) + len(reads_reverse), len(reads_forward),
               len(reads_reverse))

    def _ratio(self, mult, div):
        try:
            return(mult/div)
        except:
            return(0.0)
