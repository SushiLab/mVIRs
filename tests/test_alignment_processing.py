import unittest
from tempfile import TemporaryDirectory
import io

from mVIRs.containers import MappedRead, Mapping
from mVIRs.alignment_processing import (
    mapped_reads_from_alignment,
    _calculate_orientation,
    get_paired_alignments,
    _estimate_insert_size,
    find_oprs,
    find_clipped_reads
)

class TestMappedReadsFromAlignment(unittest.TestCase):
    def setUp(self) -> None:
        self.small_bam = "tests/data/small.sam"

    def test_mapped_reads_from_alignment(self):
        mapped = mapped_reads_from_alignment(self.small_bam)
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

class TestGetPairedReads(unittest.TestCase):
    def setUp(self) -> None:
        self.read_name = 'ERR4552622.32'
        self.test_case = MappedRead(name=self.read_name)
        self.test_case["R1"] = Mapping(False, "SalmonellaLT2", 0, 151,
                                       151, [(0, 151)], [(0, 151)])

        self.test_case["R2"] = Mapping(True, "SalmonellaLT2", 2399129, 2399276,
                                       147, [(0, 147)], [(2399129, 2399276)])

    def test_get_paired_reads(self):
        pe_aln = get_paired_alignments({self.read_name: self.test_case})[self.read_name]
        self.assertEqual(pe_aln[0].iss, 2399276)
        self.assertEqual(pe_aln[0].orientation, 'PAIREDEND')

class TestEstimateInsertSizes(unittest.TestCase):

    def test_estimate_insert_size(self):
        inserts = [1, 2, 10, 11, 11, 12, 14, 18]
        self.assertEqual(_estimate_insert_size(inserts, sd=2), (0, 19, 11))

class TestFindOprs(unittest.TestCase):
    def setUp(self) -> None:
        self.bam_file = "tests/data/ERR4552622_100k_mVIRs.bam"
        self.expected_oprs = "tests/expected/ERR4552622_100k_mVIRs.oprs"

    def test_find_oprs(self):
        with TemporaryDirectory() as tmp:
            output_oprs = tmp + "/test.oprs"
            find_oprs(self.bam_file, output_oprs,
                      min_coverage=0.8, min_alength=45)

            with io.open(output_oprs) as out, io.open(self.expected_oprs) as expected:
                self.assertListEqual(list(out), list(expected))

class TestFindClipped(unittest.TestCase):
    def setUp(self):
        self.bam_file = "tests/data/ERR4552622_100k_mVIRs.bam"
        self.expected_clipped = "tests/expected/ERR4552622_100k_mVIRs.clipped"

    def test_find_clipped(self):
        mapped = mapped_reads_from_alignment(self.bam_file, min_coverage=0, min_alength=0,
                                             extended_info=True)
        with TemporaryDirectory() as tmp:
            output_clipped = tmp + "/test.clipped"
            find_clipped_reads(mapped, output_clipped)

            with io.open(output_clipped) as out, io.open(self.expected_clipped) as expected:
                self.assertListEqual(list(out), list(expected))
