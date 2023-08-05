import unittest

from mVIRs.alignment_processing import mapped_reads_from_alignment

class TestAlignmentProcessing(unittest.TestCase):
    def setUp(self) -> None:
        self.bam_file = "tests/data/small.sam"

    def test_mapped_reads_from_alignment(self):
        mapped = mapped_reads_from_alignment(self.bam_file)
        self.assertEqual(len(mapped["ERR4552622.32"]["R1"]), 2)
        self.assertEqual(len(mapped["ERR4552622.13"]["R2"]), 1)
