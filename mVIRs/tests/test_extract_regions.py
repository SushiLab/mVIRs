import unittest
from tempfile import TemporaryDirectory
import io
from time import time

from mVIRs.extract_regions import denoise_softclips
from ..extract_utils import read_clipped_file, read_oprs_file as oprs_cython

# class TestClippedFile(unittest.TestCase):
#     def setUp(self) -> None:
#         self.clipped_file = "mVIRs/tests/expected/ERR4552622_100k_mVIRs.clipped"


# class TestOprsFile(unittest.TestCase):
#     def setUp(self) -> None:
#         self.opr_file = "mVIRs/tests/expected/ERR4552622_100k_mVIRs.oprs"
#     def test_oprs_parsing(self):
#         start = time()
#         read_oprs_file(self.opr_file)
#         print(f"Python: {time() - start}")

#         start = time()
#         oprs_cython(self.opr_file)
#         print(f"Cython: {time() - start}")

#         start = time()
#         read_oprs_file(self.opr_file)
#         print(f"Python: {time() - start}")


class TestDenoiseSoftclips(unittest.TestCase):
    def test_empty_input(self):
        self.assertEqual(denoise_softclips({}, 10), {})

    def test_single_position(self):
        soft_clipped_positions = {(10, 'A'): 2}
        expected_output = {(10, 'A'): 2}
        self.assertEqual(denoise_softclips(soft_clipped_positions, 10), expected_output)

    def test_multiple_positions_same_scaffold(self):
        soft_clipped_positions = {(10, 'A'): 2, (15, 'A'): 3, (20, 'A'): 1}
        expected_output = {(10, 'A'): 2, (15, 'A'): 4}
        self.assertEqual(denoise_softclips(soft_clipped_positions, 5), expected_output)

    def test_multiple_positions_different_scaffolds(self):
        soft_clipped_positions = {(10, 'A'): 2, (15, 'B'): 3, (20, 'C'): 1}
        expected_output = {(10, 'A'): 2, (15, 'B'): 3, (20, 'C'): 1}
        self.assertEqual(denoise_softclips(soft_clipped_positions, 5), expected_output)

    def test_multiple_positions_overlapping_ranges(self):
        soft_clipped_positions = {(10, 'A'): 2, (15, 'A'): 3, (20, 'A'): 1, (25, 'A'): 2}
        expected_output = {(15, 'A'): 8}
        self.assertEqual(denoise_softclips(soft_clipped_positions, 10), expected_output)

    def test_multiple_positions_multiple_winners(self):
        soft_clipped_positions = {(10, 'A'): 2, (15, 'A'): 3, (20, 'A'): 1, (25, 'A'): 2, (30, 'A'): 2}
        expected_output = {(15, 'A'): 8, (30, 'A'): 2}
        self.assertEqual(denoise_softclips(soft_clipped_positions, 10), expected_output)


