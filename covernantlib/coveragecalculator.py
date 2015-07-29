import pysam


class CoverageCalculator(object):

    def __init__(self, read_count_splitting=True, uniqueley_aligned_only=False,
                 first_base_only=False, paired_end=False):
        self._read_count_splitting = read_count_splitting
        self._uniqueley_aligned_only = uniqueley_aligned_only
        self._first_base_only = first_base_only
        self._coverage_add_function = self._select_coverage_add_function()
        self._coverages = {}
        self._paired_end = paired_end
        self.used_alignmets = 0

    def ref_seq_and_coverages(self, bam_path):
        bam = self._open_bam_file(bam_path)
        for ref_seq, length in zip(bam.references, bam.lengths):
            self._init_coverage_list(length)
            self._calc_coverage(ref_seq, bam)
            yield(ref_seq, self._coverages)

    def _init_coverage_list(self, length):
        for strand in ["forward", "reverse"]:
            self._coverages[strand] = [0.0] * length

    def _calc_coverage(self, ref_seq, bam):
        if not self._paired_end:
            self._calc_coverage_single_end(ref_seq, bam)
        else:
            self._calc_coverage_paired_end(ref_seq, bam)
            
    def _calc_coverage_single_end(self, ref_seq, bam):
        for alignment in bam.fetch(ref_seq):
            self._add_coverage_of_single_end_reads(alignment)
            
    def _add_coverage_of_single_end_reads(self, alignment):
        number_of_hits = self._number_of_hits(alignment)
        if self._uniqueley_aligned_only is True and number_of_hits != 1:
            return
        # Note: No translation from bam coordinates to python
        # list coorindates is needed.
        start = alignment.pos
        end = alignment.aend
        increment = self._calc_increment_values(number_of_hits)
        if end is None or start is None:
            return
        self.used_alignmets += 1
        self._coverage_add_function(alignment, increment, start, end)

    def _calc_coverage_paired_end(self, ref_seq, bam):
        """Can handle single and paired end reads. In case of paired end
        reads the reads are treates as one fragment."""
        read_pairs_by_qname = {}
        for alignment in bam.fetch(ref_seq):
            if not alignment.is_proper_pair:
                self._add_coverage_of_single_end_reads(alignment)
            else:
                # Collect all read pairs - this might be memory intensive!
                index = 0
                if alignment.is_read2:
                    index = 1
                read_pairs_by_qname.setdefault(
                    alignment.qname, [None, None])
                read_pairs_by_qname[alignment.qname][index] = alignment
        for pair in read_pairs_by_qname.values():
            # TO FIX - raises Segementation fault:
            # if None in pair:
            #    continue
            read_1, read_2 = pair
            try:
                start = min(read_1.pos, read_1.aend, read_2.pos, read_2.aend)
            except AttributeError:
                continue
            end = max(read_1.pos, read_1.aend, read_2.pos, read_2.aend)
            number_of_hits = self._number_of_hits(read_1)
            increment = self._calc_increment_values(number_of_hits)
            self._coverage_add_function(read_1, increment, start, end)
            self.used_alignmets += 1

    def _number_of_hits(self, alingment):
            try:
                return dict(alingment.tags)["NH"]
            except KeyError:
                return 1

    def _calc_increment_values(self, number_of_hits):
        # Normalize coverage increment by number of read alignments
        # per read
        if self._read_count_splitting is True:
            return 1.0 / float(number_of_hits)
        else:
            return 1.0
            
    def _select_coverage_add_function(self):
        if self._first_base_only is False:
            return self._add_whole_alignment_coverage
        else:
            return self._add_first_base_coverage

    def _open_bam_file(self, bam_file):
        return pysam.AlignmentFile(bam_file)

    def _close_bam_fh(self, bam_fh):
        bam_fh.close()

    def _add_whole_alignment_coverage(self, alignment, increment, start, end):
        if alignment.is_reverse is False:
            self._coverages["forward"][start:end] = [
                coverage + increment for coverage in
                self._coverages["forward"][start:end]]
        else:
            self._coverages["reverse"][start:end] = [
                coverage - increment for coverage in
                self._coverages["reverse"][start:end]]

    def _add_first_base_coverage(self, alignment, increment, start, end):
        if alignment.is_reverse is False:
            self._coverages["forward"][start] = self._coverages[
                "forward"][start] + increment
        else:
            self._coverages["reverse"][end-1] = self._coverages[
                "reverse"][end-1] - increment

