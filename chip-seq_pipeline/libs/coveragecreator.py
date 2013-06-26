from libs.sam import SamParser

class CoverageCreator(object):

    def __init__(self, samtools_bin="samtools"):
        self._sam_parser = SamParser(samtools_bin)
        self.elements_and_coverages = {}

    def init_coverage_lists(self, bam_file):
        for ref_seq, length in self._sam_parser.ref_seq_ids_and_lengths_bam(
            bam_file).items():
            self.elements_and_coverages[ref_seq] = [0.0] * length

    def count_coverage(self, bam_file, read_count_splitting=True,
                       uniqueley_mapped_only=False):
        for entry in self._sam_parser.entries_bam(bam_file):
            if uniqueley_mapped_only and entry.number_of_hits_as_int != 1:
                continue
            # Here a translation from 1-based system (SAM) to a
            # 0-based system (python lists) takes place. Due to this
            # each position is decreased by one. To cover the full
            # range of the end postion would need to be increased by
            # one. The substraction and addition result in a change of
            # zero.
            start = entry.start - 1
            end = entry.end
            # Normalize coverage increment by number of read mappings
            # per read
            if read_count_splitting:
                increment = 1.0 / float(entry.number_of_hits_as_int)
            else:
                increment = 1.0
            self.elements_and_coverages[entry.reference][
                start:end] = [
                    coverage + increment for coverage in
                    self.elements_and_coverages[entry.reference][
                        start:end]]
