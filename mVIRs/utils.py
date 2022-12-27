import logging
import sys
from pysam import FastxFile


def startup():
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=logging.INFO)
    logging.info('Starting mVIRs')


def shutdown(status=0):
    logging.info('Finishing mVIRs')
    sys.exit(status)


def check_sequences(r1_filepath, r2_filepath, n_headers=250):
    if r1_filepath == r2_filepath:
        logging.error(f'Input read files can not be the same file. Quitting')
        shutdown(1)

    r1_fasta, r2_fasta = FastxFile(r1_filepath), FastxFile(r2_filepath)
    counter = 0
    for f_entry, r_entry in zip(r1_fasta, r2_fasta):
        if f_entry.name != r_entry.name:
            logging.error(f'Names of input reads do not match. ({h1} != {h2}). Check if read files belong together. Quitting')
            shutdown(1)
        counter += 1
        if counter == n_headers:
            break