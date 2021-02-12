import pysam
import pathlib
import collections
import statistics
from typing import List, Dict, Tuple, Counter, DefaultDict, Iterable, IO, Generator
import logging
import argparse
import sys
import subprocess
import gzip



def add_tags(alignedSegment: pysam.AlignedSegment) -> pysam.AlignedSegment:
    """ Takes an AlignedSegment and add percent identity and alignment length as tags
    alignment length = MID
    mismatches = NM
    percent identity = (MID - NM) / MID
    The percent identity is a value between 0.0 and 1.0
    If the segment is unmapped then it is returned as with a percent identity of 0
    and an alignment length of 0.
    :param alignedSegment: The pysam AlignedSegment object
    :return: alignedSegment: The updated pysam AlignedSegment object
    """

    # Assuming that if the id tag is present that the other tags are also there.
    if alignedSegment.has_tag('id'):
        return alignedSegment
    if alignedSegment.is_unmapped:
        alignedSegment.set_tag('id', 0.0, 'f')
        alignedSegment.set_tag('al', 0, 'i')
        alignedSegment.set_tag('qc', 0.0, 'f')
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

def align(forward_read_file, reversed_read_file, bwa_ref_name, threads, out_bam_file, min_coverage, min_percid, min_alength):
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

    :return:
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
        in_bam_file_handle = pysam.AlignmentFile(process.stdout, 'rb')
        if not temp_bam_file_handle:
            temp_bam_file_handle = pysam.AlignmentFile(temp_bam_file, "wb", template=in_bam_file_handle)

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
        returncode: int = subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        raise Exception(e)
    if returncode != 0:
        raise(Exception(f'Command: {command} failed with return code {returncode}'))
        shutdown(1)

    pathlib.Path(temp_bam_file).unlink()
    logging.info('Finished alignment step')
















PAlignment = collections.namedtuple('PAlignment',
                                    'iss, ref revr1 revr2 score startr1 endr1 startr2 endr2 orientation')
SAMLine = collections.namedtuple('SAMLine', 'rev ref rstart rend score')


def _calc_orientation(revr1: bool, revr2: bool, posr1: int, posr2: int) -> str:
    """
    Based on orientation of both reads and their alignment start positions, estimate if a read is PE/SAME/OPR
    :param revr1:
    :param revr2:
    :param posr1:
    :param posr2:
    :return: PAIREDEND, SAME, OPR
    """
    orientation: str = 'PAIREDEND'
    if revr1 == revr2:
        orientation = 'SAME'
        return orientation
    fwpos: int = 0
    revpos: int = 0
    if revr1:
        revpos = posr1
        fwpos = posr2
    else:
        revpos = posr2
        fwpos = posr1
    if revpos < fwpos:
        orientation = 'OPR'
    return orientation


def _generate_paired_alignments(insert2alignments: Dict[str, Dict[str, List[SAMLine]]],
                                cutoff_bestscore: float = 0.95) -> Generator[
    Tuple[str, DefaultDict[str, List[PAlignment]]], None, None]:
    """
    Iterate over all inserts and generate pairs with a cross product. Pairs with best score are returned.
    Careful --> throws away single end mappers. Considers only paired reads


    :param insert2alignments:
    :return:
    """

    cnt_singleend: int = 0
    cnt_pairedend: int = 0
    query: str
    alignments: Dict[str, List[SAMLine]]
    for query, alignments in insert2alignments.items():
        if len(alignments) < 2:
            cnt_singleend += 1
            continue
        else:
            cnt_pairedend += 1
            matches: DefaultDict[str, List[PAlignment]] = collections.defaultdict(list)
            samematches: List[PAlignment] = []
            bestscore: int = 0
            alnr1: List[SAMLine]
            alnr2: List[SAMLine]

            for alnr1 in alignments['R1']:
                for alnr2 in alignments['R2']:
                    if alnr1.ref != alnr2.ref:  # if both reads align to 2 different references we continue
                        continue
                    orientation: str = _calc_orientation(alnr1.rev, alnr2.rev, alnr1.rstart, alnr2.rstart)
                    isp = abs(alnr1.rstart - alnr2.rstart)
                    positions = [alnr1.rstart, alnr2.rstart, alnr1.rend, alnr2.rend]
                    iss = abs(max(positions) - min(positions))
                    if orientation == 'SAME':
                        samematches.append(
                            PAlignment(iss=iss, revr1=alnr1.rev, revr2=alnr2.rev,
                                       score=alnr1.score + alnr2.score, startr1=alnr1.rstart, endr1=alnr1.rend, startr2=alnr2.rstart, endr2=alnr2.rend,
                                       orientation=orientation, ref=alnr1.ref))
                        continue
                    score: int = alnr1.score + alnr2.score
                    if bestscore < score:
                        bestscore = score
                    matches[alnr1.score + alnr2.score].append(
                        PAlignment(iss=iss, revr1=alnr1.rev, revr2=alnr2.rev, score=score,
                                   startr1=alnr1.rstart, endr1=alnr1.rend, startr2=alnr2.rstart,
                                   endr2=alnr2.rend, orientation=orientation, ref=alnr1.ref))


            minscore: int = int(bestscore * cutoff_bestscore)
            scoreinsertsizesortedmatches: DefaultDict[str, List[PAlignment]] = []
            score: int
            for score in sorted(matches.keys(), reverse=True):
                if score < minscore:
                    continue
                match: List[PAlignment] = matches[score]
                tmpmatch: Iterable[PAlignment] = sorted(match, key=lambda x: x.iss)
                scoreinsertsizesortedmatches = scoreinsertsizesortedmatches + tmpmatch
            samematch: PAlignment
            for samematch in samematches:
                if samematch.score >= minscore:
                    scoreinsertsizesortedmatches.append(samematch)
            if len(scoreinsertsizesortedmatches) == 0:
                continue
            yield (query, scoreinsertsizesortedmatches)


def _read_bam_file(bam_file: pathlib.Path, max_sam_lines: int = -1) -> Dict[str, Dict[str, List[SAMLine]]]:
    """
    Reads a bam file into memory and splits by R1/R2/S. Each insert is either singleton or paired.
    Singletons are removed
    :param bamfile:
    :return: Dictionary with inserts and the location where reads align
    """
    logging.info('Start reading alignments from file:\t{}'.format(bam_file))
    bf: pysam.libcalignmentfile.AlignmentFile = pysam.AlignmentFile(str(bam_file), "rb")

    insert2alignments: Dict[str, Dict[str, List[SAMLine]]] = collections.defaultdict(
        lambda: collections.defaultdict(list))
    cnt: int
    alignment: pysam.AlignedSegment
    for cnt, alignment in enumerate(bf):
        if cnt % 1000000 == 0:
            logging.info('Alignments read:\t{}\tInserts found:\t{}'.format(format(cnt, ',d'), format(len(insert2alignments), ',d')))
        qname: str = alignment.qname
        splits: List[str] = qname.rsplit('/', 1)
        qname: str = splits[0]

        if splits[1] == 'S':  # Remove singletons
            continue
        ascore: int = alignment.get_tag('AS')

        orientation: str = splits[1]
        reverse: bool = True if alignment.is_reverse else False
        refname: str = alignment.reference_name
        refstart: int = alignment.reference_start
        refend: int = alignment.reference_end


        insert2alignments[qname][orientation].append(
            SAMLine(rev=reverse, ref=refname, rstart=refstart, rend=refend, score=ascore))
        if max_sam_lines <= cnt:
            break
    logging.info('Alignments read:\t{}\tInserts found:\t{}'.format(format(cnt, ',d'), format(len(insert2alignments), ',d')))

    bf.close()
    return insert2alignments


def _estimate_insert_size(insert2alignments: Dict[str, Dict[str, List[SAMLine]]]) -> Tuple[int, int, int]:
    """
    Goes through UNIQUE PAIRED END alignments and estimates the median insert size.
    :param insert2alignments:
    :return:
    """

    cnter: Counter = collections.Counter()
    insertsizes: List[int] = []

    alignments: List[PAlignment]
    for (_, alignments) in _generate_paired_alignments(insert2alignments):
        if len(alignments) == 1 and alignments[0].orientation == 'PAIREDEND':
            insertsizes.append(alignments[0].iss)

    so: List[int] = insertsizes
    minvaltmp: int = -1
    maxvaltmp: int = -1
    iteration = 0
    while True:
        iteration += 1
        logging.info(
            'Estimating insert size by mean/stdev convergence. Iteration:\t{}. Min insert size:\t{}. Max insert size:\t{}. Mean insert size:\t{}. Stdev insert size:\t{}'.format(
                iteration, format(min(so), ',d'), format(max(so), ',d'), format(int(statistics.mean(so)), ',d'), format(int(statistics.stdev(so)), ',d')))
        tmp: List[int] = []
        mean: int = int(statistics.mean(so))
        stdev: int = int(statistics.pstdev(so))
        print(stdev)
        minval: int = mean - 7 * stdev
        if minval < 0:
            minval = 0
        maxval: int = mean + 7 * stdev
        val: int
        for val in so:
            if minval > val:
                continue
            if maxval < val:
                continue
            tmp.append(val)
        so: List[int] = tmp

        if minvaltmp == minval and maxval == maxvaltmp:
            break

        minvaltmp: int = minval
        maxvaltmp: int = maxval


    # For mapping Illumina short-insert reads to the human genome, x is about 6-7 sigma away from the mean
    # http://bio-bwa.sourceforge.net/bwa.shtml section Estimating Insert Size Distribution

    minvaltmp: int = mean - 7 * int(statistics.pstdev(so))

    if minvaltmp != 0:
        minvaltmp = 0

    maxvaltmp: int = mean + 7 * int(statistics.pstdev(so))
    logging.info(
        'Estimating insert size by mean/stdev convergence. Iteration:\t{}. Min insert size:\t{}. Max insert size:\t{}. Mean insert size:\t{}. Stdev insert size:\t{}'.format(
            iteration, format(minvaltmp, ',d'), format(maxvaltmp, ',d'), format(int(statistics.mean(so)), ',d'),
            format(int(statistics.stdev(so)), ',d')))

    return minvaltmp, maxvaltmp, int(statistics.median(so))


def _calc_primary_paired_alignment(insert2alignments: Dict[str, Dict[str, List[SAMLine]]], insertsize: int,
                                   minreasonable_insertsize: int, maxreasonable_insertsize: int) -> Generator[
    Tuple[str, List[PAlignment], bool], None, None]:
    """
    Go through all paired alignments and find one that has a reasonable good score and is close to the insert size.
    1. if there is a PE alignment with reasonable insert size --> take that one (and all other equally good PE with reasonable IS)
    2. Otherwise take the best one
    3. If there are multiple best ones return all and state add the false flag to indicate that alignment could not be resolved to a single best alignment

    :param insert2alignments:
    :param insertsize:
    :return:
    """

    qname: str
    alignments: List[PAlignment]

    for (qname, alignments) in _generate_paired_alignments(insert2alignments):


        if len(alignments) == 1:  # This is a unique alignment. No filtering required.
            yield (qname, alignments, True)
        else:
            # Check if there are any alignments with PE and reasonable insert size
            pe_candidates: List[PAlignment] = []
            alignment: PAlignment
            for alignment in alignments:
                if alignment.orientation == 'PAIREDEND' and alignment.iss >= minreasonable_insertsize and alignment.iss <= maxreasonable_insertsize:
                    pe_candidates.append(alignment)
            if len(pe_candidates) != 0:
                # if there are alignments with PE and reasonable insertsize then pick the one(s) with the best score
                bestscore: int = -1
                pe_bestalns: List[PAlignment] = []
                alignment: PAlignment
                for alignment in sorted(pe_candidates, key=lambda x: x.score, reverse=True):
                    if bestscore == -1:
                        bestscore = alignment.score
                    if alignment.score == bestscore:
                        pe_bestalns.append(alignment)
                # and yield the best one(s)
                if len(pe_bestalns) != 0:
                    yield (qname, pe_bestalns, True)
                    continue
            else:
                yield (qname, alignments, False)
                continue



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













def startup():
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=logging.INFO)
    logging.info('Starting mVIRs')


def shutdown(status=0):
    logging.info('Finishing mVIRs')
    sys.exit(status)


def find(bam_file, opr_file) -> None:
    """
    Find OPRs in aligned inserts. This includes a couple steps:

    1. Read the original bam file into memory and pair --> insert2alignments. Singleton inserts are tossed.

    2. Estimate the insert size from the unique paired mappers. We do that to get a reasonable insert size for filtering
    potentially good (paired-end) or strange (OPR) alignments.

    3. calculate the best primary alignments.

    3.1  if there is a unique best alignment --> take that one
    3.2. elif there is paired-end alignment with reasonable insert size and within 5% of the best score --> take that one
    3.3. elif there is an OPR that is the best alignment for this insert --> take that one
    3.4  elif take every best mapper

    The output is a file that reports for each insert for which we can identify a best mapper and has either an:
    - an unreasonable insert size
    - an unnormal orientation (OPR <---- ---->) or (SAME ---> --->)

    :param bam_file:
    :param opr_file:
    :return:
    """

    logging.info('Start OPR finding step')
    logging.info('Input BAM File:\t{}'.format(bam_file))
    logging.info('Output OPR File:\t{}'.format(opr_file))


    # START DEBUG PARAMETERS

    max_sam_lines = 10000000000
    # END DEBUG PARAMETERS

    insert2alignments: Dict[str, Dict[str, List[SAMLine]]] = _read_bam_file(bam_file, max_sam_lines)
    logging.info('Start estimating insert size from paired end alignments.')
    minreasonable_insertsize, maxreasonable_insertsize, estimated_insertsize = _estimate_insert_size(insert2alignments)
    logging.info('Finished estimating insert size from paired end alignments.')
    logging.info('Min reasonable insert size:\t{}'.format(minreasonable_insertsize))
    logging.info('Max reasonable insert size:\t{}'.format(maxreasonable_insertsize))





    with open(opr_file, 'w') as handle:
        handle.write(
            '#MIN_REASONABLE_INSERTSIZE={}\n#MAX_REASONABLE_INSERTSIZE={}\n#READNAME\tREFERENCE\tINSERT_SIZE\tR1_ORIENTATION\tR2_ORIENTATION\tBWA_SCORE\tR1_START\tR2_START\tR1_ALNLENGTH\tR2_ALNLENGTH\tREAD_ORIENTATION\n'.format(
                minreasonable_insertsize, maxreasonable_insertsize))



        template: str = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'


        logging.info('Start screening for OPRs and Paired-End inserts with unreasonable insert size.')
        alncnt = 0
        for alncnt, (query, alignments, found) in enumerate(
                _calc_primary_paired_alignment(insert2alignments, estimated_insertsize, minreasonable_insertsize,
                                               maxreasonable_insertsize)):
            if found:
                for cnt, alignment in enumerate(alignments):
                    printstring = template.format(query, alignment.ref, alignment.iss,
                                                  get_read_orientation(alignment.revr1),
                                                  get_read_orientation(alignment.revr2), alignment.score, alignment.startr1,
                                                  alignment.startr2, alignment.endr1 - alignment.startr1, alignment.endr2 - alignment.startr2,
                                                  alignment.orientation)
                    if alignment.orientation == 'PAIREDEND' and alignment.iss >= minreasonable_insertsize and alignment.iss <= maxreasonable_insertsize:  # pe with reasonable insert size
                        x = 0
                    else:  # pe with unreasonable insert size and opr/same
                        printstring = printstring.replace('PAIREDEND', 'IPR')
                        handle.write(printstring)
            else:  # undefined
                continue

            if alncnt % 100000 == 0:
                logging.info('Paired inserts processed:\t{} / {}'.format(format(alncnt, ',d'), format(len(insert2alignments), ',d')))
        logging.info('Paired inserts processed:\t{} / {}'.format(format(alncnt, ',d'), format(len(insert2alignments), ',d')))
        logging.info('Singleton inserts tossed:\t{} / {}'.format(format(len(insert2alignments) - alncnt, ',d'), format(len(insert2alignments), ',d')))
        logging.info('Finished screening for OPRs and Paired-End inserts with unreasonable insert size.')



def read_seq_file(seq_file):
    lines = []
    if seq_file.endswith('gz'):
        with gzip.open(seq_file, 'rt') as handle:
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


def oprs():

    parser = argparse.ArgumentParser(description='Align paired reads against a reference database and find outward orientated paired reads (OPRs).', prog='mvirs oprs')
    parser.add_argument('i1', action='store', help='Forward reads file. FastQ or FastA files supported. Input can be gzipped')
    parser.add_argument('i2', action='store', help='Reverse reads file. FastQ or FastA files supported. Input can be gzipped')
    parser.add_argument('r', action='store', help='BWA reference. Has to be created upfront using the mvirs oprs command.')
    parser.add_argument('b', action='store', help='Output BAM file. File with all filtered alignments created by aligning forward and reverse reads against the reference database.')
    parser.add_argument('o', action='store',help='Output OPR file. File with all OPRs and IPRs found in the alignment file.')
    parser.add_argument('-t', action='store', dest='threads',help='Number of threads to use. (Default = 1)',type=int, default=1)

    try:
        args = parser.parse_args(sys.argv[2:])
    except:
        shutdown(1)



    forward_read_file = args.i1
    reversed_read_file = args.i2

    out_bam_file = args.b
    bwa_ref_name = args.r
    opr_file = pathlib.Path(args.o)

    min_percid = 0.97
    remove_unmapped = True
    min_coverage = 0.8
    min_alength = 45

    threads = args.threads

    if not pathlib.Path(forward_read_file).exists() or not pathlib.Path(forward_read_file).is_file():
        logging.error(f'The forward reads file does not exist: {forward_read_file}. Quitting')
        shutdown(1)
    if not pathlib.Path(reversed_read_file).exists() or not pathlib.Path(reversed_read_file).is_file():
        logging.error(f'The reverse reads file does not exist: {reversed_read_file}. Quitting')
        shutdown(1)

    # check the names of sequences match

    required_index_files = [bwa_ref_name + suffix for suffix in ['.bwt', '.pac', '.ann', '.amb', '.sa']]
    is_reference_missing = False
    for rif in required_index_files:
        if not pathlib.Path(rif).exists() or not pathlib.Path(rif).is_file():
            is_reference_missing = True
    if is_reference_missing:
        logging.error(f'The bwa index files are missing. Please rerun mvirs index. Quitting')

    check_sequences(forward_read_file, reversed_read_file)


    if threads <= 0:
        raise argparse.ArgumentTypeError('Number of threads has to be >0'.format(threads))
        shutdown(1)

    align(forward_read_file, reversed_read_file, bwa_ref_name, threads, out_bam_file, min_coverage, min_percid,
          min_alength)
    logging.info('\n\n\n\n')
    find(out_bam_file, opr_file)


def index():
    parser = argparse.ArgumentParser(description='Generates the BWA index that is required for the oprs command.', prog='mvirs index')
    parser.add_argument('r', action='store', help='Input FastA or FastQ file for index building. Gzipped input allowed.')
    try:
        args = parser.parse_args(sys.argv[2:])
    except:
        shutdown(1)
    seq_file = args.r
    if not pathlib.Path(seq_file).exists() or not pathlib.Path(seq_file).is_file():
        logging.error(f'The input file for bwa index building does not exist: {seq_file}. Quitting')
        shutdown(1)
    logging.info(f'Start building bwa index on {seq_file}')
    command = f'bwa index {seq_file}'

    try:
        returncode: int = subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        raise Exception(e)
    if returncode != 0:
        raise(Exception(f'Command: {command} failed with return code {returncode}'))
        shutdown(1)
    logging.info(f'Successfully built index on {seq_file}')
    shutdown(0)


def main():
    startup()

    parser = argparse.ArgumentParser(
        description='Bioinformatic toolkit for finding prophages in sequencing data', usage='''mvirs <command> [<args>]

    Command options
        oprs    align reads to a reference and report outward orientated alignments (OPRs)
        index   index the reference database using bwa
    ''')
    parser.add_argument('command', help='Subcommand to run: oprs | index')

    args = parser.parse_args(sys.argv[1:2])

    if args.command == 'oprs':
        oprs()
    if args.command == 'index':
        index()
    else:
        print('Unrecognized command')
        parser.print_usage()
        shutdown(1)
    shutdown(0)



if __name__ == '__main__':
    main()
