from pysam.libcalignmentfile import AlignmentFile
import unittest
from mVIRs.containers import MappedRead, Mapping

class TestMapping(unittest.TestCase):
    def setUp(self):
        self.bam_file = AlignmentFile("tests/data/small.sam", "rb")
        self.aln1 = next(self.bam_file)
        self.aln2 = next(self.bam_file)
        self.aln3 = next(self.bam_file)
        self.aln4 = next(self.bam_file)

    def test_mapping_extended_false(self):
        mapping = Mapping(self.aln1, extended=False)
        self.assertEqual(mapping.rev, False)
        self.assertEqual(mapping.ref, "SalmonellaLT2")
        self.assertEqual(mapping.rstart, 0)
        self.assertEqual(mapping.rend, 151)
        self.assertEqual(mapping.score, 151)
        self.assertEqual(mapping.blocks, None)
        self.assertEqual(mapping.cigartuples, None)

    def test_mapping_extended_true(self):
        mapping = Mapping(self.aln2, extended=True)
        self.assertEqual(mapping.rev, True)
        self.assertEqual(mapping.ref, "SalmonellaLT2")
        self.assertEqual(mapping.rstart, 1)
        self.assertEqual(mapping.rend, 148)
        self.assertEqual(mapping.score, 147)
        self.assertEqual(mapping.blocks, [(1, 148)])
        self.assertEqual(mapping.cigartuples, [(0, 147)])

    def test_mapped_read(self):
        readname = self.aln1.qname.rsplit('/', 1)[0]
        read = MappedRead(name=readname)
        for aligned_seg in [self.aln1, self.aln2,
                            self.aln3, self.aln4]:

            data_tmp = aligned_seg.qname.rsplit('/', 1)
            readname = data_tmp[0]
            orientation = data_tmp[1]
            read[orientation] = Mapping(aligned_seg, extended=False)

        self.assertEqual(read.name, "ERR4552622.32")
        self.assertEqual(len(read.mapping["R1"]), 2)
        self.assertEqual(len(read.mapping["R2"]), 2)




