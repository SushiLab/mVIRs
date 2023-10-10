import os
import unittest
import urllib.request
import tempfile
import filecmp

from mVIRs.alignment import index_genome
from mVIRs.mvirs import _execute_oprs

class TestRegression(unittest.TestCase):

    def test_full_pipeline(self):

        required_files = dict(
        r1_file = 'https://sunagawalab.ethz.ch/share/MVIRS_TEST//ERR4552622_100k_1.fastq.gz',
        r2_file = 'https://sunagawalab.ethz.ch/share/MVIRS_TEST//ERR4552622_100k_2.fastq.gz',
        reference_genome = 'https://sunagawalab.ethz.ch/share/MVIRS_TEST//np_salmoLT2.fasta.gz'
        )

        with tempfile.TemporaryDirectory() as tmp:
            for name in required_files:
                url = required_files[name]
                localpath = os.path.join(tmp, name.split('/')[-1])
                urllib.request.urlretrieve(url, localpath)
                required_files[name] = localpath

            # index genome
            bwa_ref_name = index_genome(required_files['reference_genome'], tmp)
            out_bam_file = os.path.join(tmp, 'ERR4552622_100k_mVIRs.bam')
            opr_file = os.path.join(tmp, 'ERR4552622_100k_mVIRs.oprs')
            clipped_file = os.path.join(tmp, 'ERR4552622_100k_mVIRs.clipped')
            output_fasta_file = os.path.join(tmp, 'ERR4552622_100k_mVIRs.fasta')

            _execute_oprs(required_files['r1_file'], required_files['r2_file'],
                          out_bam_file, bwa_ref_name, required_files['reference_genome'],
                          opr_file, clipped_file, output_fasta_file)

            filecmp.cmp(opr_file, 'tests/expected/ERR4552622_100k_mVIRs.oprs')
            filecmp.cmp(clipped_file, 'tests/expected/ERR4552622_100k_mVIRs.clipped')
            filecmp.cmp(output_fasta_file, 'tests/expected/ERR4552622_100k_mVIRs.fasta')