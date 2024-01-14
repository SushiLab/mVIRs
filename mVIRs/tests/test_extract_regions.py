import unittest
from tempfile import TemporaryDirectory
import io
from time import time

from mVIRs.extract_regions import extract_regions, read_clipped_file as read_clipped
from ..extract_utils import read_clipped_file

class TestClippedFile(unittest.TestCase):
    def setUp(self) -> None:
        self.clipped_file = "mVIRs/tests/expected/ERR4552622_100k_mVIRs.clipped"

    def test_read_clipped_file(self):
        start = time()
        clipped_dict2, softclipped_pos2 = read_clipped(self.clipped_file)
        print(f"read_clipped took {time() - start} seconds")

        start = time()
        clipped_dict, softclipped_pos = read_clipped_file(self.clipped_file)
        print(f"Cython read_clipped_file took {time() - start} seconds")



# class TestExtractRegions(unittest.TestCase):

#     def setUp(self) -> None:
#         self.opr_file = "tests/expected/ERR4552622_100k_mVIRs.oprs"
#         self.clipped_file = "tests/expected/ERR4552622_100k_mVIRs.clipped"
#         self.reference = "tests/data/np_salmoLT2.fasta.gz"
#         self.expected_fasta = "tests/expected/ERR4552622_100k_mVIRs.fasta"

#     def test_extract_regions(self):
#        with TemporaryDirectory() as tmp:
#            output_fasta = tmp + "/test.fasta"
#            extract_regions(self.clipped_file, self.opr_file,  self.reference, output_fasta,
#                            min_length=4000, max_length=800_000, allow_complete_scaffolds=True)

#            with io.open(output_fasta) as out, io.open(self.expected_fasta) as expected:
#                self.assertListEqual(list(out), list(expected))


