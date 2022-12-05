import logging
import sys
from mVIRs import read_seq_file


def startup():
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=logging.INFO)
    logging.info('Starting mVIRs')


def shutdown(status=0):
    logging.info('Finishing mVIRs')
    sys.exit(status)


def check_sequences(r1_file, r2_file):
    if r1_file == r2_file:
        logging.error(f'Input read files can not be the same file. Quitting')
        shutdown(1)

    seq_headers_r1 = read_seq_file(r1_file)
    seq_headers_r2 = read_seq_file(r2_file)
    for h1, h2 in zip(seq_headers_r1, seq_headers_r2):
        if h1 != h2:
            logging.error(f'Names of input reads do not match. ({h1} != {h2}). Check if read files belong together. Quitting')
            shutdown(1)