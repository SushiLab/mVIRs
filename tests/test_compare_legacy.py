import unittest
from mVIRs.oprs import insertize_bamfile_by_name, _generate_paired_alignments
from mVIRs.alignment_processing import (
    mapped_reads_from_alignment,
    get_paired_alignments
)

class TestLegacy(unittest.TestCase):
    def setUp(self):
        # this part is different because of data structures used
        bam_file = "/nfs/cds-peta/exports/biol_micro_cds_gr_sunagawa/scratch/vbezshapkin/mVIRs/tests/data/ERR4552622_100k_mVIRs.sam"
        self.orig = [x for x in insertize_bamfile_by_name(bam_file)]
        self.mapped = mapped_reads_from_alignment(bam_file)

    def test_paired_alignments(self):
        self.p_alignments = dict([x for x in _generate_paired_alignments(self.orig)])
        self.paired = get_paired_alignments(self.mapped)
        self.assertEqual(self.p_alignments, self.paired)
