import argparse
import glob
import logging
import os
import pathlib
import subprocess
import sys
import urllib.request
import urllib.response

from .alignment import align
from .index import index_genome
from .utils import check_sequences, shutdown, startup

from mVIRs import ( 
    find_clipped_reads,
    find_oprs,
    extract_regions
) 

VERSION = '1.1.1'


def _execute_oprs(forward_read_file, reversed_read_file, out_bam_file, bwa_ref_name, opr_file, clipped_file, output_fasta_file, minlength_report=4000, maxlength_report=800000, allow_fl_report=True, threads=1):
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
        shutdown(1)


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
    extract_regions(clipped_file, opr_file, bwa_ref_name, output_fasta_file,
                    minmvirlength=minlength_report, maxmvirlength=maxlength_report, allow_complete_scaffolds=allow_fl_report)

class CapitalisedHelpFormatter(argparse.HelpFormatter):
    def add_usage(self, usage, actions, groups, prefix=None):
        if prefix is None:
            prefix = ''
        return super(CapitalisedHelpFormatter, self).add_usage(usage, actions, groups, prefix)


def index():
    parser = argparse.ArgumentParser(description='Generates the BWA index that is required for the oprs command.', usage=f'''
Program: mVIRs - Localisation of inducible prophages using NGS data
Version: {VERSION}
Reference: Zünd, Ruscheweyh, et al. 
High throughput sequencing provides exact genomic locations of inducible 
prophages and accurate phage-to-host ratios in gut microbial strains. 
Microbiome (2021). doi:10.1186/s40168-021-01033-w    

Usage: mvirs index [options]

    Input:
        -f  FILE   Reference FASTA file. Can be gzipped. [Required]
        -o  PATH   Output folder [Required]
    ''', formatter_class=CapitalisedHelpFormatter,add_help=False)
    parser.add_argument('-f', action='store', help='Input FastA or FastQ file for index building. Gzipped input allowed.', required=True)
    parser.add_argument('-o', action='store', help='Output folder for bwa index files.', required=False)
    try:
        args = parser.parse_args(sys.argv[2:])
    except:
        shutdown(1)
    
    seq_file = args.f
    if not pathlib.Path(seq_file).exists() or not pathlib.Path(seq_file).is_file():
        logging.error(f'The input file for bwa index building does not exist: {seq_file}. Quitting')
        shutdown(1)
    
    index_genome(seq_file, args.o)



def oprs():

    parser = argparse.ArgumentParser(description='Align reads against a reference and find OPRs and IPRs.', usage=f'''
Program: mVIRs - Localisation of inducible prophages using NGS data
Version: {VERSION}
Reference: Zünd, Ruscheweyh, et al. 
High throughput sequencing provides exact genomic locations of inducible 
prophages and accurate phage-to-host ratios in gut microbial strains. 
Microbiome (2021). doi:10.1186/s40168-021-01033-w    

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
    parser.add_argument('-f', action='store', help='Forward reads file. Can be gzipped', required=True, dest='forward')
    parser.add_argument('-r', action='store', help='Reverse reads file. Can be gzipped', required=True, dest='reverse')
    parser.add_argument('-db', action='store', help='BWA reference. Has to be created upfront', required=True, dest='db')

    parser.add_argument('-ml', action='store', type=int, help='Minimum length to extract. [4000]', default=4000, required=False, dest='ml')
    parser.add_argument('-ML', action='store', type=int, help='Maximum length to extract. [800000]', default=800000, required=False, dest='ML')
    parser.add_argument('-m', action='store_true', help='Allow full scaffolds to be reported', required=False, dest='afs') 
    parser.add_argument('-o', action='store', help='Output folder', required=True, dest='output')
    parser.add_argument('-t', action='store', dest='threads',help='Number of threads to use. (Default = 1)',type=int, default=1)
    
    
    try:
        args = parser.parse_args(sys.argv[2:])
    except:
        shutdown(1)

    forward_read_file = args.forward
    reversed_read_file = args.reverse

    bwa_ref_folder = args.db
    minlength_report = args.ml
    maxlength_report = args.ML
    allow_fl_report = args.afs
    threads = args.threads

    required_index_ext = ['.bwt', '.pac', '.ann', '.amb', '.sa']

    files = glob.glob(bwa_ref_folder + "/*")

    bwa_ref_name = None
    for file in files:
       
        if file.endswith(required_index_ext[0]):
            filename =  file.split("/")[-1].rsplit(".", 1)[0]
            bwa_ref_name = os.path.join(bwa_ref_folder, filename)
            
    
    if not bwa_ref_name: 
        logging.error(f'Could not find BWA index files in {bwa_ref_folder}. Quitting.')
        shutdown(1)
    
    out_file_prefix = os.path.join(args.output, filename)
    out_bam_file = out_file_prefix + '.bam'
    opr_file = pathlib.Path(out_file_prefix + '.oprs')
    clipped_file = pathlib.Path(out_file_prefix + '.clipped')
    output_fasta_file = pathlib.Path(out_file_prefix + '.fasta')
    
    required_index_files = [bwa_ref_name + suffix for suffix in required_index_ext]
    is_reference_missing = False
    for rif in required_index_files:
        if not pathlib.Path(rif).exists() or not pathlib.Path(rif).is_file():
            is_reference_missing = True
    if is_reference_missing:
        logging.error(f'Could not find BWA index files in {bwa_ref_folder}. Quitting')
        shutdown(1)
     
    # make output folder
    output_folder = pathlib.Path(args.output)
    output_folder.mkdir(parents=True, exist_ok=True)

    _execute_oprs(forward_read_file, reversed_read_file, out_bam_file, bwa_ref_name, opr_file, clipped_file, output_fasta_file, minlength_report, maxlength_report, allow_fl_report, threads)


def test():
    parser = argparse.ArgumentParser(description='Run mVIRs on a public dataset', usage=f'''
Program: mVIRs - Localisation of inducible prophages using NGS data
Version: {VERSION}
Reference: Zünd, Ruscheweyh, et al. 
High throughput sequencing provides exact genomic locations of inducible 
prophages and accurate phage-to-host ratios in gut microbial strains. 
Microbiome (2021). doi:10.1186/s40168-021-01033-w    
Usage: mvirs test [options]

    Input:
        -o  PATH   Output folder [Required]
    ''', formatter_class=CapitalisedHelpFormatter,add_help=False)
    parser.add_argument('-o', action='store', help='Output folder', required=True)
    try:
        args = parser.parse_args(sys.argv[2:])
    except:
        shutdown(1)
    output_folder = pathlib.Path(args.o)
    # make output folder

    output_folder.mkdir(parents=True, exist_ok=True)
    #download reads/ref
    logging.info('Downloading reference and read files (22MB)')
    remote_r1_file = 'https://sunagawalab.ethz.ch/share/MVIRS_TEST//ERR4552622_100k_1.fastq.gz'
    local_r1_file = str(output_folder) + '/ERR4552622_100k_1.fastq.gz'
    remote_r2_file = 'https://sunagawalab.ethz.ch/share/MVIRS_TEST//ERR4552622_100k_2.fastq.gz'
    local_r2_file = str(output_folder)+ '/ERR4552622_100k_2.fastq.gz'
    remote_reference_file = 'https://sunagawalab.ethz.ch/share/MVIRS_TEST//np_salmoLT2.fasta.gz'
    local_reference_file = str(output_folder) + '/np_salmoLT2.fasta.gz'
    urllib.request.urlretrieve(remote_r1_file, local_r1_file)
    urllib.request.urlretrieve(remote_r2_file, local_r2_file)
    urllib.request.urlretrieve(remote_reference_file, local_reference_file)
    logging.info('Finished downloading')
    logging.info('Building mVIRs index')
    # build index
    command = f'bwa index {local_reference_file}'

    try:
        returncode: int = subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        raise Exception(e)
    if returncode != 0:
        raise(Exception(f'Command: {command} failed with return code {returncode}'))
        shutdown(1)
    logging.info(f'Successfully built index on {local_reference_file}')
    # run mvirs
    logging.info('Run mVIRs test')
    forward_read_file = local_r1_file
    reversed_read_file = local_r2_file

    out_bam_file = str(output_folder) + '/ERR4552622_100k_mVIRs.bam'
    bwa_ref_name = local_reference_file
    opr_file = str(output_folder) + '/ERR4552622_100k_mVIRs.oprs'
    clipped_file = str(output_folder) + '/ERR4552622_100k_mVIRs.clipped'
    output_fasta_file = str(output_folder) + '/ERR4552622_100k_mVIRs.fasta'
    _execute_oprs(forward_read_file, reversed_read_file, out_bam_file, bwa_ref_name, opr_file, clipped_file, output_fasta_file)


def main():
    startup()

    parser = argparse.ArgumentParser(
        description='Bioinformatic toolkit for finding prophages in sequencing data', usage=f'''
        
Program: mVIRs - Localisation of inducible prophages using NGS data
Version: {VERSION}
Reference: Zünd, Ruscheweyh, et al. 
High throughput sequencing provides exact genomic locations of inducible 
prophages and accurate phage-to-host ratios in gut microbial strains. 
Microbiome (2021). doi:10.1186/s40168-021-01033-w    

Usage: mvirs <command> [options]
Command:

    index   create index files for reference used in the 
            mvirs oprs routine
            
    oprs    align reads against reference and used clipped
            alignment positions and OPRs to extract potential
            prophages
            
    test    run mVIRs for a public dataset

    ''', formatter_class=CapitalisedHelpFormatter,add_help=False)

    parser.add_argument('command', help='Subcommand to run: oprs')

    args = parser.parse_args(sys.argv[1:2])

    if args.command == 'oprs':
        oprs()
    elif args.command == 'index':
        index()
    elif args.command == 'test':
        test()
    else:
        print('Unrecognized command')
        parser.print_usage()
        shutdown(1)
    shutdown(0)


if __name__ == '__main__':
    main()