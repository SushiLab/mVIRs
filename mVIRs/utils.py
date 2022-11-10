import gzip
import logging
import sys


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


def read_seq_file(seq_file):
    lines = []
    if seq_file.endswith('gz'):
        with gzip.open(seq_file, 'rt') as handle:
            for line in handle:
                lines.append(line.strip())
                if len(lines) == 1000:
                    break
    else:
        with open(seq_file) as handle:
            for line in handle:
                lines.append(line.strip())
                if len(lines) == 1000:
                    break

    modulo = 2
    if lines[0].startswith('@'):
        modulo = 4
    seq_headers = []
    for cnt, line in enumerate(lines):
        if cnt % modulo == 0:
            seq_headers.append(line.split()[0])
    return seq_headers


def get_read_orientation(rev: bool) -> str:
    """
    Helper function to print if an alignments is forward or reverse
    :param rev:
    :return:
    """
    if rev:
        return 'reverse'
    else:
        return 'forward'


def load_fasta(sequence_file):
    '''
    Read a fasta file and put it into a dictionary
    :param sequence_file:
    :return:
    '''
    if sequence_file.endswith('.gz'):
        handle = gzip.open(sequence_file, 'rt')
    else:
        handle = open(sequence_file)
    sequences = {}
    current_header = None
    for line in handle:
        line = line.strip().split()[0]
        if len(line) == 0:
            continue
        if line.startswith('>'):
            line = line[1:]
            current_header = line
            sequences[current_header] = []
        else:
            sequences[current_header].append(line)

    handle.close()
    sequences2 = {}
    for header, sequence in sequences.items():
        tmp_seq = ''.join(sequence)
        sequences2[header] = tmp_seq
    return sequences2