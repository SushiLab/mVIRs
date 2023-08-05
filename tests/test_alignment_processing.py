import unittest

from mVIRs.alignment_processing import (
    mapped_reads_from_alignment,
    _calculate_orientation
)

class TestMappedReadsFromAlignment(unittest.TestCase):
    def setUp(self) -> None:
        self.bam_file = "tests/data/small.sam"

    def test_mapped_reads_from_alignment(self):
        mapped = mapped_reads_from_alignment(self.bam_file)
        self.assertEqual(len(mapped["ERR4552622.32"]["R1"]), 2)
        self.assertEqual(len(mapped["ERR4552622.13"]["R2"]), 1)

class TestCalculateOrientation(unittest.TestCase):
    def setUp(self) -> None:
        # same direction
        self.same = (False, 1, False, 5)
        # the first read is reversed
        self.opr_reverse_r1 = (True, 1, False, 5)
        # the second read is reversed
        self.opr_reverse_r2 = (False, 5, True, 1)
        # inward oriented == pairedend
        self.paired_end_1 = (True, 5, False, 1)
        self.paired_end_2 = (False, 1, True, 5)



    def test_calculate_orientation(self):
        self.assertEqual(_calculate_orientation(*self.same), "SAME")
        self.assertEqual(_calculate_orientation(*self.opr_reverse_r1), "OPR")
        self.assertEqual(_calculate_orientation(*self.opr_reverse_r2), "OPR")
        self.assertEqual(_calculate_orientation(*self.paired_end_1), "PAIREDEND")
        self.assertEqual(_calculate_orientation(*self.paired_end_2), "PAIREDEND")