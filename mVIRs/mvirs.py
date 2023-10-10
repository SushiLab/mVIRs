import argparse
import logging
import os
import pathlib
import sys

from .alignment import index_genome, align
from .utils import check_sequences, shutdown, startup

from mVIRs import add_info
from mVIRs.oprs import (
    find_clipped_reads,
    find_oprs,
)

from mVIRs.extract_regions import extract_regions





def _execute_oprs(forward_read_file, reversed_read_file, out_bam_file, bwa_ref_name, reference_genome,
                   opr_file, clipped_file, output_fasta_file, minlength_report=4000, maxlength_report=800000, allow_fl_report=True, threads=1):
    min_coverage = 0.8
    min_alength = 45

    ### CHECKS START
    if minlength_report < 0:
        raise argparse.ArgumentTypeError('-ml has to be >0'.format(minlength_report))
        shutdown(1)

    if maxlength_report < 0:
        raise argparse.ArgumentTypeError('-ML has to be >0'.format(maxlength_report))
        shutdown(1)
    if threads <= 0:
        raise argparse.ArgumentTypeError('Number of threads has to be > 0'.format(threads))

    if not pathlib.Path(forward_read_file).exists() or not pathlib.Path(forward_read_file).is_file():
        logging.error(f'The forward reads file does not exist: {forward_read_file}. Quitting')
        shutdown(1)
    if not pathlib.Path(reversed_read_file).exists() or not pathlib.Path(reversed_read_file).is_file():
        logging.error(f'The reverse reads file does not exist: {reversed_read_file}. Quitting')
        shutdown(1)

    check_sequences(forward_read_file, reversed_read_file)

    ### CHECKS END


    align(forward_read_file, reversed_read_file, bwa_ref_name, out_bam_file, threads)
    find_clipped_reads(out_bam_file, clipped_file)
    find_oprs(out_bam_file, opr_file, min_coverage, min_alength)
    extract_regions(clipped_file, opr_file, reference_genome, output_fasta_file,
                    min_length=minlength_report, max_length=maxlength_report,
                    allow_complete_scaffolds=allow_fl_report)

    # remove bwa index files
    for suffix in ['.amb', '.ann', '.bwt', '.pac', '.sa']:
        os.remove(bwa_ref_name.with_suffix(suffix))

class CapitalisedHelpFormatter(argparse.HelpFormatter):
    def add_usage(self, usage, actions, groups, prefix=None):
        if prefix is None:
            prefix = ''
        return super(CapitalisedHelpFormatter, self).add_usage(usage, actions, groups, prefix)

def check_positive(value):
    if value <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return value

def oprs():

    parser = argparse.ArgumentParser(description='Align reads against a reference and find OPRs and IPRs.',
                                     usage=f'''{add_info}

Usage: mvirs oprs [options]

    Input:
        -f  FILE   Forward reads file. FastA/Q. Can be gzipped. [Required]
        -r  FILE   Reverse reads file. FastA/Q. Can be gzipped. [Required]
        -db FILE   Reference database file (prefix) created by mvirs index. [Required]

    Output:
        -o  PATH   Output folder. [Required]

    Options:
        -t  INT    Number of threads. [1]
        -ml INT    Minimum sequence length for extraction.. [4000]
        -ML INT    Maximum sequence length for extraction.. [800000]
        -m         Allow full contigs/scaffolds/chromosomes to be reported
                   (When OPRs and clipped reads are found at the start and
                   end of contigs/scaffolds/
    ''', formatter_class=CapitalisedHelpFormatter,add_help=False)
    parser.add_argument('-ref', action='store', help='Reference FASTA file. Can be gzipped', required=True)
    parser.add_argument('-f', action='store', help='Forward reads file. Can be gzipped', required=True, dest='forward')
    parser.add_argument('-r', action='store', help='Reverse reads file. Can be gzipped', required=True, dest='reverse')

    parser.add_argument('-ml', action='store', type=int, help='Minimum length to extract. [4000]', default=4000, required=False, dest='ml')
    parser.add_argument('-ML', action='store', type=int, help='Maximum length to extract. [800000]', default=800000, required=False, dest='ML')
    parser.add_argument('-m', action='store_true', help='Allow full scaffolds to be reported', required=False, dest='afs')
    parser.add_argument('-o', action='store', help='Output folder', required=True, dest='output')
    parser.add_argument('-t', action='store', dest='threads',help='Number of threads to use. (Default = 1)',type=int, default=1)


    try:
        args = parser.parse_args(sys.argv[2:])
    except:
        shutdown(1)

    genome_file = args.ref
    if not pathlib.Path(genome_file).exists() or not pathlib.Path(genome_file).is_file():
        raise FileNotFoundError(f'The input file for bwa index building does not exist: {genome_file}. Quitting')

    forward_read_file = args.forward
    reversed_read_file = args.reverse

    min_length = args.ml
    max_length = args.ML
    allow_fl_report = args.afs
    threads = args.threads

    bwa_ref_name = index_genome(genome_file, args.output)
    filename = pathlib.Path(bwa_ref_name)

    out_file_prefix = os.path.join(args.output, filename)
    out_bam_file = out_file_prefix + '.bam'
    opr_file = pathlib.Path(out_file_prefix + '.oprs')
    clipped_file = pathlib.Path(out_file_prefix + '.clipped')
    output_fasta_file = pathlib.Path(out_file_prefix + '.fasta')

    # make output folder
    output_folder = pathlib.Path(args.output)
    output_folder.mkdir(parents=True, exist_ok=True)

    _execute_oprs(forward_read_file, reversed_read_file,
                  str(out_bam_file), str(bwa_ref_name),
                  str(genome_file),
                  str(opr_file),
                  str(clipped_file), str(output_fasta_file),
                  min_length, max_length, allow_fl_report, threads)


def main():
    startup()

    parser = argparse.ArgumentParser(
        description='Bioinformatic toolkit for finding prophages in sequencing data',
        usage=f'''{add_info}

Usage: mvirs <command> [options]
Command:

    oprs    align reads against reference and used clipped
            alignment positions and OPRs to extract potential
            prophages

    ''', formatter_class=CapitalisedHelpFormatter,add_help=False)

    parser.add_argument('command', help='Subcommand to run: oprs')

    args = parser.parse_args(sys.argv[1:2])

    if args.command == 'oprs':
        oprs()
    else:
        print('Unrecognized command')
        parser.print_usage()
        shutdown(1)
    shutdown(0)


if __name__ == '__main__':
    main()