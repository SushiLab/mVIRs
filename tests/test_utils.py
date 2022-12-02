import pytest
import filecmp
from pathlib import Path

from mVIRs import (
    find_clipped_reads,
    find_oprs,
    extract_regions 
)


@pytest.fixture
def bam_file():
    return "tests/data/ERR4552622_100k_mVIRs.bam"

@pytest.fixture
def reference_fasta():
    return "tests/data/np_salmoLT2.fasta.gz"

# @pytest.fixture
# def test_dir():
#     path = Path("tests/data/test_out")
#     path.mkdir(exist_ok=True)
#     return path

@pytest.fixture(scope="session")
def test_dir(tmp_path_factory):
    test_dir = tmp_path_factory.mktemp("tmp")
    return test_dir

@pytest.fixture
def out_prefix(bam_file):
    return bam_file.rsplit("/", 1)[1].replace(".bam", "")

@pytest.fixture
def clipped_file(out_prefix, test_dir):
    clipped_filepath = test_dir / (out_prefix + ".clipped")
    return clipped_filepath

@pytest.fixture
def oprs_file(out_prefix, test_dir):
    oprs_filepath = test_dir / (out_prefix + ".oprs")
    return oprs_filepath

@pytest.fixture
def output_fasta(out_prefix, test_dir):
    output_fasta = test_dir / (out_prefix + ".fasta")
    return output_fasta

# def test_load_fasta():
#     small_fasta = "tests/data/small.fasta"
#     loaded_fasta = load_fasta(small_fasta)
#     assert loaded_fasta == {"SalmonellaLT2": "CGGG",
#                             "SalmonellaLT2_plasmid": "ACAG"}

def test_find_clipped_reads(bam_file, clipped_file):
    # Find clipped reads
    find_clipped_reads(str(bam_file), str(clipped_file))
    clipped_expected = "tests/expected/ERR4552622_100k_mVIRs.clipped"   
    assert filecmp.cmp(clipped_file, clipped_expected, shallow=False)


def test_find_oprs(bam_file, oprs_file):
    find_oprs(str(bam_file), str(oprs_file), min_coverage=0.8, min_alength=45)
    oprs_expected = "tests/expected/ERR4552622_100k_mVIRs.oprs"
    assert filecmp.cmp(oprs_file, oprs_expected, shallow=False)


def test_extract_regions(clipped_file, oprs_file, reference_fasta, output_fasta):
    extract_regions(str(clipped_file), 
                    str(oprs_file), 
                    str(reference_fasta), 
                    str(output_fasta),
                    minmvirlength=1000, 
                    maxmvirlength=1000000, 
                    allow_complete_scaffolds=True)
    output_fasta_expected = "tests/expected/ERR4552622_100k_mVIRs.fasta"    
    assert filecmp.cmp(output_fasta, output_fasta_expected, shallow=False)






    



    