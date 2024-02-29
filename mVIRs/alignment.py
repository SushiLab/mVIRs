import logging
import pathlib
import subprocess

import pysam

from .utils import shutdown


def index_genome(seq_file: str, output_folder: str) -> None:
    """
    Builds a bwa index for a given genome file.

    Args:
        seq_file (str): The path to the genome file.
        output_folder (str): The output folder for the index files.

    Returns:
        index_path (str): The path to the index file.
    """
    logging.info(f'Start building bwa index on {seq_file}')
    # make output folder
    out_folder = pathlib.Path(output_folder)
    out_folder.mkdir(parents=True, exist_ok=True)
    index_path = (pathlib.Path(output_folder) / seq_file.split("/")[-1]).absolute()
    command = f'bwa index {seq_file} -p {index_path}'

    try:
        returncode: int = subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        raise Exception(e)
    if returncode != 0:
        raise(Exception(f'Command: {command} failed with return code {returncode}'))

    logging.info(f'Successfully built index on {seq_file}')

    return index_path



def add_tags(alignedSegment: pysam.AlignedSegment) -> pysam.AlignedSegment:
    """
    Takes an AlignedSegment and add percent identity and alignment length as tags
    alignment length = MID
    mismatches = NM
    percent identity = (MID - NM) / MID
    The percent identity is a value between 0.0 and 1.0
    If the segment is unmapped then it is returned as with a percent identity of 0
    and an alignment length of 0.

    Args:
        alignedSegment (pysam.AlignedSegment): The aligned segment to add tags to.

    Returns:
        pysam.AlignedSegment: The aligned segment with the added tags.
    """

    # Assuming that if the id tag is present that the other tags are also there.
    if alignedSegment.has_tag('id'):
        return alignedSegment

    alnlength = sum(alignedSegment.get_cigar_stats()[0][0:3])

    query_covered_bases = sum(alignedSegment.get_cigar_stats()[0][0:2])

    query_length = alignedSegment.infer_read_length()
    mismatches = alignedSegment.get_tag('NM')
    percid = (alnlength - mismatches) / float(alnlength)
    qcov = query_covered_bases / float(query_length)
    alignedSegment.set_tag('id', percid, 'f')
    alignedSegment.set_tag('qc', qcov, 'f')
    alignedSegment.set_tag('al', alnlength, 'i')
    return alignedSegment


def align(forward_read_file: str,
          reversed_read_file: str,
          bwa_ref_name: str,
          out_bam_file: str,
          threads: int,
          min_percid: float = 0.97,
          min_coverage: float = 0,
          min_alength: int = 0):
    """
    Takes two paired end fastq/fasta files and aligns them against a reference genome and reports sorted and filtered alignments.
    Prerequisites:
    1. R1 and R2 reads files need to have the same read names for the same insert. This can be problematic with data downloaded from SRA
    2. The bwa needs to be contructed beforehands. E.g. You want to align against the reference genome.fasta. Then you need to call "bwa index genome.fasta".
    Then you can provide genome.fasta as parameter for -r
    Execution:
    1. Align R1/R2 read files against the reference
    2. Filter alignments by 97% identity, read coverage >=80%, alignmentlength >= 45 and remove unmapped alignments
    3. Each readname from R1 file gets an /R1 tag to its readname. /R2 is added for R2 alignments.
    4. This creates a temporary bam file. This bam file is sorted with samtools by name to produce the final bam file.

    Args:
        forward_read_file (str): The path to the forward read file.
        reversed_read_file (str): The path to the reversed read file.
        bwa_ref_name (str): The path to the reference genome index file.
        out_bam_file (str): The path to the output bam file.
        threads (int): The number of threads to use for the alignment.
        min_percid (float): The minimum percent identity of the alignment to keep.
        min_coverage (float): The minimum coverage of the alignment to keep.
        min_alength (int): The minimum alignment length to keep.

    Returns:
        None
    """



    logging.info('Start alignment step')
    logging.info('\tInput files:')

    temp_bam_file = out_bam_file + '_temp.bam'
    temp_bam_file_handle = None
    logging.info(f'Executing BWA alignment:')
    for readsfile, orientation in [(forward_read_file, '/R1'), (reversed_read_file, '/R2')]:

        command = f'bwa mem -a -t {threads} {bwa_ref_name} {readsfile} | samtools view -h -F 4 -'
        logging.info(f'\tCommand executed {command}')

        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)

        # https://github.com/pysam-developers/pysam/issues/939
        save = pysam.set_verbosity(0)
        in_bam_file_handle = pysam.AlignmentFile(process.stdout, 'rb')
        if not temp_bam_file_handle:
            temp_bam_file_handle = pysam.AlignmentFile(temp_bam_file, "wb", template=in_bam_file_handle)
        pysam.set_verbosity(save)

        for record in in_bam_file_handle:
            if record.is_unmapped:
                continue
            else:
                updated_record = add_tags(record)
                if updated_record.get_tag('al') >= min_alength and updated_record.get_tag('qc') >= min_coverage and updated_record.get_tag('id') >= min_percid:
                    record.qname = ''.join([record.qname, orientation])
                    temp_bam_file_handle.write(record)

        process.stdout.close()
        return_code = process.wait()
        if return_code != 0:
            logging.error(f'BWA command failed with return code {return_code}')
            shutdown(1)

    in_bam_file_handle.close()
    temp_bam_file_handle.close()


    logging.info(f'Executing samtools sort:')
    command = f'samtools sort -n -m 4G -@ {threads} -o {out_bam_file} {temp_bam_file}'
    logging.info(f'\tCommand executed {command}')
    try:
        returncode = subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        raise Exception(e)
    if returncode != 0:
        raise(Exception(f'Command: {command} failed with return code {returncode}'))

    pathlib.Path(temp_bam_file).unlink()
    logging.info('Finished alignment step')