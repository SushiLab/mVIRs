import unittest
from tempfile import TemporaryDirectory
import io

from mVIRs.extract_regions import extract_regions


class TestExtractRegions(unittest.TestCase):

    def setUp(self) -> None:
        self.opr_file = "tests/expected/ERR4552622_100k_mVIRs.oprs"
        self.clipped_file = "tests/expected/ERR4552622_100k_mVIRs.clipped"
        self.reference = "tests/data/np_salmoLT2.fasta.gz"
        self.expected_fasta = "tests/expected/ERR4552622_100k_mVIRs.fasta"

    def test_extract_regions(self):
       with TemporaryDirectory() as tmp:
           output_fasta = tmp + "/test.fasta"
           extract_regions(self.clipped_file, self.opr_file,  self.reference, output_fasta,
                           min_length=4000, max_length=800_000, allow_complete_scaffolds=True)

           with io.open(output_fasta) as out, io.open(self.expected_fasta) as expected:
               self.assertListEqual(list(out), list(expected))
