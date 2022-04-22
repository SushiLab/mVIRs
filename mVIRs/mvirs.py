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
import urllib.request
import urllib.response


VERSION = '1.1.1'

debug = False

PYSAM_BAM_CSOFT_CLIP = 4
PYSAM_BAM_CHARD_CLIP = 5

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


# def align_individual(read_file, orientation, bwa_ref_name, threads, out_bam_file, out_map_file, min_coverage, min_percid, min_alength):
#
#
#     temp_bam_file = out_bam_file + '_temp.bam'
#     logging.info(f'Executing BWA alignment:')
#
#     output_map_file = open(out_map_file, 'w')
#
#     command = f'bwa mem -a -t {threads} {bwa_ref_name} {read_file} | samtools view -h -F 4 -'
#     logging.info(f'\tCommand executed {command}')
#
#     process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
#     in_bam_file_handle = pysam.AlignmentFile(process.stdout, 'rb')
#     temp_bam_file_handle = pysam.AlignmentFile(temp_bam_file, "wb", template=in_bam_file_handle)
#
#     for record in in_bam_file_handle:
#         if record.is_unmapped:
#             continue
#         else:
#             updated_record = add_tags(record)
#             if updated_record.get_tag('al') >= min_alength and updated_record.get_tag('qc') >= min_coverage and updated_record.get_tag('id') >= min_percid:
#
#
#                 # START CHANGE
#                 qname: str = updated_record.qname
#                 ascore: int = updated_record.get_tag('AS')
#                 reverse: bool = True if updated_record.is_reverse else False
#                 refname: str = updated_record.reference_name
#                 refstart: int = updated_record.reference_start
#                 refend: int = updated_record.reference_end
#                 # blocks: list = updated_record.get_blocks()
#                 # cigartuples: list = updated_record.cigartuples
#                 alignment_length: int = updated_record.get_tag('al')
#                 alignment_qcov: float = updated_record.get_tag('qc')
#                 alignment_id: float = updated_record.get_tag('id')
#                 output_map_file.write(f'{qname}\t{reverse}\t{refname}\t{refstart}\t{refend}\t{ascore}\t{alignment_length}\t{alignment_qcov}\t{alignment_id}\n')
#                 # END CHANGE
#                 record.qname = ''.join([record.qname, '/R', orientation])
#                 temp_bam_file_handle.write(record)
#
#     process.stdout.close()
#     return_code = process.wait()
#     if return_code != 0:
#         logging.error(f'BWA command failed with return code {return_code}')
#         shutdown(1)
#
#     in_bam_file_handle.close()
#     temp_bam_file_handle.close()
#     output_map_file.close()
#
#
#     logging.info(f'Executing samtools sort:')
#     command = f'samtools sort -n -m 4G -@ {threads} -o {out_bam_file} {temp_bam_file}'
#     logging.info(f'\tCommand executed {command}')
#     try:
#         returncode: int = subprocess.check_call(command, shell=True)
#     except subprocess.CalledProcessError as e:
#         raise Exception(e)
#     if returncode != 0:
#         raise(Exception(f'Command: {command} failed with return code {returncode}'))
#         shutdown(1)
#
#     pathlib.Path(temp_bam_file).unlink()



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
SAMLine = collections.namedtuple('SAMLine', 'rev ref rstart rend score cigartuples blocks')


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


def _generate_paired_alignments(insert2alignments,
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
    for query, alignments in insert2alignments:
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
                match = matches[score]
                tmpmatch: Iterable[PAlignment] = sorted(match, key=lambda x: x.iss)
                scoreinsertsizesortedmatches = scoreinsertsizesortedmatches + tmpmatch
            samematch: PAlignment
            for samematch in samematches:
                if samematch.score >= minscore:
                    scoreinsertsizesortedmatches.append(samematch)
            if len(scoreinsertsizesortedmatches) == 0:
                continue
            yield (query, scoreinsertsizesortedmatches)


def insertize_bamfile_by_name(bam_file: pathlib.Path, max_sam_lines: int = -1, min_coverage: float = 0.0, min_alength: int = 0, need_extended = True):
    alignments: pysam.libcalignmentfile.AlignmentFile = pysam.AlignmentFile(str(bam_file), "rb")

    current_name = None
    current_insert = collections.defaultdict(list)


    for alignment in alignments:
        data_tmp = alignment.qname.rsplit('/', 1)
        readname = data_tmp[0]
        orientation = data_tmp[1]
        ascore: int = alignment.get_tag('AS')
        reverse: bool = True if alignment.is_reverse else False
        refname: str = alignment.reference_name
        refstart: int = alignment.reference_start
        refend: int = alignment.reference_end
        blocks = None
        cigartuples = None
        if need_extended:
            blocks: list = alignment.get_blocks()
            cigartuples: list = alignment.cigartuples
        if not (alignment.get_tag('al') >= min_alength and alignment.get_tag('qc') >= min_coverage):
            continue
        samline = SAMLine(rev=reverse, ref=refname, rstart=refstart, rend=refend, score=ascore, cigartuples=cigartuples, blocks=blocks)

        if not current_name:
            current_name = readname
        if readname == current_name:
            current_insert[orientation].append(samline)
        else:
            yield current_name, current_insert
            current_name = readname
            current_insert = collections.defaultdict(list)
            current_insert[orientation].append(samline)



    if len(current_insert) != 0:
        yield current_name, current_insert


# def _read_bam_file_generator(bam_file: pathlib.Path, max_sam_lines: int = -1, min_coverage: float = 0.0, min_alength: int = 0, toss_singletons: bool = True) -> Dict[str, Dict[str, List[SAMLine]]]:
#     """
#     Reads a bam file into memory and splits by R1/R2/S. Each insert is either singleton or paired.
#     Singletons are removed
#     :param bamfile:
#     :return: Dictionary with inserts and the location where reads align
#     """
#     logging.info('Start reading alignments from file:\t{}'.format(bam_file))
#     bf: pysam.libcalignmentfile.AlignmentFile = pysam.AlignmentFile(str(bam_file), "rb")
#     total_alignments_in_insert2alignments: int = 0
#     insert2alignments: Dict[str, Dict[str, List[SAMLine]]] = collections.defaultdict(
#         lambda: collections.defaultdict(list))
#     cnt: int
#     alignment: pysam.AlignedSegment
#     for cnt, alignment in enumerate(bf):
#         if cnt % 100000 == 0:
#             logging.info('Alignments read:\t{}\tInserts found:\t{}'.format(format(cnt, ',d'), format(len(insert2alignments), ',d')))
#         qname: str = alignment.qname
#         splits: List[str] = qname.rsplit('/', 1)
#         qname: str = splits[0]
#
#         if splits[1] == 'S':  # Remove singletons
#             if toss_singletons:
#                 continue
#         ascore: int = alignment.get_tag('AS')
#
#         orientation: str = splits[1]
#         reverse: bool = True if alignment.is_reverse else False
#         refname: str = alignment.reference_name
#         refstart: int = alignment.reference_start
#         refend: int = alignment.reference_end
#         blocks: list = alignment.get_blocks()
#         cigartuples: list = alignment.cigartuples
#         if not (alignment.get_tag('al') >= min_alength and alignment.get_tag('qc') >= min_coverage):
#             continue
#
#
#         samline = SAMLine(rev=reverse, ref=refname, rstart=refstart, rend=refend, score=ascore, cigartuples=cigartuples, blocks=blocks)
#
#         insert2alignments[qname][orientation].append(samline)
#         total_alignments_in_insert2alignments += 1
#         if max_sam_lines <= cnt:
#             break
#     logging.info('Alignments read:\t{}\tInserts found:\t{}'.format(format(cnt, ',d'), format(len(insert2alignments), ',d')))
#
#     bf.close()
#     return insert2alignments, total_alignments_in_insert2alignments


def _estimate_insert_size(insert2alignments: Dict[str, Dict[str, List[SAMLine]]]) -> Tuple[int, int, int]:
    """
    Goes through UNIQUE PAIRED END alignments and estimates the median insert size.
    :param insert2alignments:
    :return:
    """

    cnter: Counter = collections.Counter()
    insertsizes: List[int] = []
    import pprint

    alignments: List[PAlignment]
    for (_, alignments) in _generate_paired_alignments(insert2alignments):

        if len(alignments) == 1 and alignments[0].orientation == 'PAIREDEND':
            insertsizes.append(alignments[0].iss)
            if len(insertsizes) > 20000:
                break

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
        #print(stdev)
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
                                   minreasonable_insertsize: int, maxreasonable_insertsize: int) -> Generator[Tuple[str, List[PAlignment], bool], None, None]:
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


def find_clipped_reads(bam_file, clipped_file) -> None:
    """
    Find clipped reads in a bam file. The goal of the clipped reads is to find start/end positions of
    activated prophages. This can be done by looking at cigar string that contain a S.
    If the S is before a M then this is a START position. If the S is after an M then that is a END
    position
    :param bam_file:
    :param clipped_file:

    :return:
    """
    logging.info('Start clipped reads finding step')
    logging.info('Input BAM File:\t{}'.format(bam_file))
    logging.info('Output Clipped File:\t{}'.format(clipped_file))


    # START DEBUG PARAMETERS
    if debug:
        max_sam_lines = 10000
    else:
        max_sam_lines = sys.maxsize


    # END DEBUG PARAMETERS
    min_coverage = 0.0
    min_alength = 0
    toss_singletons = False

    insert2alignments_generator = insertize_bamfile_by_name(bam_file, max_sam_lines, min_coverage, min_alength, need_extended=True)
    alignments_seen: int = 0
    alignments_softclipped: int = 0
    alignments_hardclipped: int = 0

    clipped_reads = collections.defaultdict(list)
    for insert, ori_2_alignments in insert2alignments_generator:
        for ori, alignments in ori_2_alignments.items():
            for alignment in alignments:
                alignments_seen += 1
                if alignments_seen % 500000 == 0:
                    logging.info(f'{format(alignments_seen, ",d")} - {format(alignments_softclipped, ",d")} - {format(alignments_hardclipped, ",d")} | Alignments - Softclips - Hardclips')
                if len(alignment.cigartuples) > 1:
                    softclipped: bool = True if len(set(filter(lambda cg: cg == PYSAM_BAM_CSOFT_CLIP, map(lambda cigar: cigar[0], alignment.cigartuples)))) == 1 else False
                    hardclipped: bool = True if len(set(filter(lambda cg: cg == PYSAM_BAM_CHARD_CLIP,map(lambda cigar: cigar[0],alignment.cigartuples)))) == 1 else False

                    if softclipped:
                        clip_location = [0, 0]
                        if alignment.cigartuples[0][0] == PYSAM_BAM_CSOFT_CLIP:
                            clip_location[0] = 1
                        if alignment.cigartuples[-1][0] == PYSAM_BAM_CSOFT_CLIP:
                            clip_location[1] = 1

                        if sum(clip_location) == 2:
                            # logging.info(f'Insert {insert}/{ori} mapped with softclips in front and end. Ignoring')
                            continue
                        alignments_softclipped += 1

                        if clip_location[0] == 1:
                            front = True
                            direction = '->'
                            startpos = alignment.blocks[0][0]
                        else:
                            end = True
                            direction = '<-'
                            startpos = alignment.blocks[-1][1]
                        clipped_reads[(insert, ori)].append(('S', direction, startpos, alignment.ref))

                    if hardclipped and not softclipped:
                        clip_location = [0, 0]
                        if alignment.cigartuples[0][0] == PYSAM_BAM_CHARD_CLIP:
                            clip_location[0] = 1
                        if alignment.cigartuples[-1][0] == PYSAM_BAM_CHARD_CLIP:
                            clip_location[1] = 1

                        if sum(clip_location) == 2:
                            #logging.info(f'Insert {insert}/{ori} mapped with hardclips in front and end. Ignoring')
                            continue

                        alignments_hardclipped += 1

                        if clip_location[0] == 1:
                            front = True
                            direction = '->'
                            startpos = alignment.blocks[0][0]
                        else:
                            end = True
                            direction = '<-'
                            startpos = alignment.blocks[-1][1]
                        clipped_reads[(insert, ori)].append(('H', direction, startpos, alignment.ref))

    logging.info(f'{format(alignments_seen, ",d")} - {format(alignments_softclipped, ",d")} - {format(alignments_hardclipped, ",d")} | Alignments - Softclips - Hardclips')
    logging.info(f'Keeping only paired soft/hard-clips and writing unfiltered soft and filtered hard-clips to {clipped_file}.')
    out_file = open(clipped_file, 'w')

    hardclips_written = 0
    softclips_written = 0
    out_file.write('#INSERT\tREADORIENTATION\tHARD/SOFTCLIP\tDIRECTION\tPOSITION\tSCAFFOLD\n')
    for (insert, ori), clipped_alignments in clipped_reads.items():
        softclipped = [aln for aln in clipped_alignments if aln[0] == 'S']
        hardclipped = [aln for aln in clipped_alignments if aln[0] == 'H']
        if len(softclipped) == 0: # A read with an alignment without softclip. Ignoring
            continue
        # there can be multiple hardclips
        # - Candidate hardclips have to face the opposite direction
        # - Candidate hardclips also need to face at each other. So the -> needs to have the smaller coordinate
        # There can me multiple hardclips even after filtering for direction
        softclipped_ori = softclipped[0][1]
        softclipped_coordinate = softclipped[0][2]
        softclipped_ref = softclipped[0][3]
        hardclipped_filtered = []
        for hardclip in hardclipped:
            hardclip_ori = hardclip[1]
            hardclip_coordinate = hardclip[2]
            hardclip_ref = hardclip[3]
            if hardclip_ref != softclipped_ref:
                continue
            if softclipped_ori == hardclip_ori:
                continue
            if softclipped_ori == '->' and softclipped_coordinate > hardclip_coordinate:
                continue
            else:
                hardclipped_filtered.append(hardclip)
        hardclips_written = hardclips_written + len(hardclipped_filtered)
        softclips_written = softclips_written + len(softclipped)
        for aln in (softclipped + hardclipped_filtered):
            out_file.write(f'{insert}\t{ori}\t{aln[0]}\t{aln[1]}\t{aln[2]}\t{aln[3]}\n')


    out_file.close()
    logging.info(f'Wrote {format(softclips_written, ",d")} soft and {format(hardclips_written, ",d")} hard-clips.')




def find_oprs(out_bam_file, opr_file, min_coverage, min_alength) -> None:
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
    logging.info('Input BAM File:\t{}'.format(out_bam_file))
    logging.info('Output OPR File:\t{}'.format(opr_file))


    # START DEBUG PARAMETERS
    if debug:
        max_sam_lines = 10000
    else:
        max_sam_lines = sys.maxsize
    # END DEBUG PARAMETERS

    insert2alignments = None
    logging.info('Start estimating insert size from paired end alignments.')
    minreasonable_insertsize, maxreasonable_insertsize, estimated_insertsize = _estimate_insert_size(insertize_bamfile_by_name(out_bam_file, max_sam_lines, min_coverage, min_alength, need_extended=False))
    logging.info('Finished estimating insert size from paired end alignments.')
    logging.info('Min reasonable insert size:\t{}'.format(minreasonable_insertsize))
    logging.info('Max reasonable insert size:\t{}'.format(maxreasonable_insertsize))





    with open(opr_file, 'w') as handle:
        handle.write(
            '#MIN_REASONABLE_INSERTSIZE={}\n#MAX_REASONABLE_INSERTSIZE={}\n#ESTIMATED_INSERTSIZE={}\n#READNAME\tREFERENCE\tINSERT_SIZE\tR1_ORIENTATION\tR2_ORIENTATION\tBWA_SCORE\tR1_START\tR2_START\tR1_ALNLENGTH\tR2_ALNLENGTH\tREAD_ORIENTATION\n'.format(
                minreasonable_insertsize, maxreasonable_insertsize, estimated_insertsize))



        template: str = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'

        insert2alignments_gen = insertize_bamfile_by_name(out_bam_file, max_sam_lines, min_coverage, min_alength, need_extended=False)
        logging.info('Start screening for OPRs and Paired-End inserts with unreasonable insert size.')
        alncnt = 0
        for alncnt, (query, alignments, found) in enumerate(
                _calc_primary_paired_alignment(insert2alignments_gen, estimated_insertsize, minreasonable_insertsize,
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

            if alncnt % 500000 == 0:
                logging.info('Paired inserts processed:\t{}'.format(format(alncnt, ',d')))
        logging.info('Paired inserts processed:\t{}'.format(format(alncnt, ',d')))
        logging.info('Finished screening for OPRs and Paired-End inserts with unreasonable insert size.')


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




def extract_regions(clipped_file, opr_file, reference_fasta_file, output_fasta_file, minmvirlength=1000, maxmvirlength=1000000, allow_complete_scaffolds=True):
    clipped_reads = collections.defaultdict(list)
    soft_clipped_positions = collections.Counter()
    soft_to_hardclip_pairs = collections.Counter() # (start, stop, scaffold) --> Count
    oprs_start_to_stop = collections.Counter() # (start, stop, scaffold) --> Count
    max_reasonable_insert_size = 0
    estimated_insert_size = 0

    softclip_range = 20



    # OPRS have an unprecise location of start and end position of virus but start and end are connected
    # Soft-Hard pairs are partly unprecise (But more precise then OPRs) but connect start and end
    # Soft-clips are precise and plenty but they don't connect start with end

    # 1. Create a start-end map with abundance from soft-hard pairs
    # 2. Add OPRs to the start-end map whereever they potentially fit +- 500 bp
    # 3. Get a more precise location with softclips

    logging.info('Finding potential viruses in the genome')

    logging.info('Reading reference fasta')
    reference_header_2_sequence = load_fasta(reference_fasta_file)

    with open(clipped_file) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            splits = line.strip().split()
            clipped_reads[(splits[0], splits[1])].append((splits[2], splits[3], int(splits[4]), splits[5]))
            if splits[2] == 'S':
                soft_clipped_positions[(int(splits[4]), splits[5])] += 1

    logging.info(f'Start finding start/end positions of non-continous alignment regions using clipped and OPR alignments')

    logging.info(f'Creating initial start/end positions using hard-soft alignment pairs')
    logging.info(f'Start denoising {len(soft_clipped_positions)} soft clipped positions from {sum(soft_clipped_positions.values())} reads.')
    updated_soft_clipped_positions = collections.Counter()

    for (softclip_position, softclip_scaffold), softclip_count in soft_clipped_positions.most_common():

        candidates = collections.Counter()
        for potential_pos in range(softclip_position-softclip_range, softclip_position+softclip_range):
            if (potential_pos, softclip_scaffold) in updated_soft_clipped_positions:
                candidates[(potential_pos, softclip_scaffold)] = updated_soft_clipped_positions[(potential_pos, softclip_scaffold)]


        if len(candidates) == 0:
            updated_soft_clipped_positions[(softclip_position, softclip_scaffold)] = softclip_count
        elif len(candidates) == 1:
            winner = list(candidates.keys())[0]
            updated_soft_clipped_positions[winner] += softclip_count
        else:
            winner = candidates.most_common()[0]
            updated_soft_clipped_positions[winner] += softclip_count
    cnt_tmp = len(soft_clipped_positions)
    sum_tmp = sum(soft_clipped_positions.values())
    logging.info(
        f'Denoising soft clipped reads finished. {len(updated_soft_clipped_positions)} ({int(len(updated_soft_clipped_positions) * 100.0 / cnt_tmp)}%) positions from {sum(updated_soft_clipped_positions.values())} ({int(sum(updated_soft_clipped_positions.values()) * 100.0 / sum_tmp)}%) reads were kept.')
    soft_clipped_positions = collections.Counter()
    for (softclip_position, softclip_scaffold), softclip_count in updated_soft_clipped_positions.most_common():
        if softclip_count > 1:
            soft_clipped_positions[(softclip_position, softclip_scaffold)] = softclip_count
    logging.info(f'Removing singleton positions finished. {len(soft_clipped_positions)} ({int(len(soft_clipped_positions)*100.0/cnt_tmp)}%) positions from {sum(soft_clipped_positions.values())} ({int(sum(soft_clipped_positions.values())*100.0/sum_tmp)}%) reads were kept.')






    '''
    Pairing hard/soft clips
    '''
    for insert, alignments in clipped_reads.items():
        softclipped = [aln for aln in alignments if aln[0] == 'S'][0]
        hardclipped = [aln for aln in alignments if aln[0] == 'H']
        if len(hardclipped) == 0:
            continue
        scaffold = softclipped[3]
        softclipped_pos = softclipped[2]
        for hardclip in hardclipped:
            hardclip_pos = hardclip[2]
            if hardclip_pos > softclipped_pos:
                soft_to_hardclip_pairs[(softclipped_pos, hardclip_pos, scaffold)] += 1
            else:
                soft_to_hardclip_pairs[(hardclip_pos, softclipped_pos,  scaffold)] += 1
    logging.info(f'Found {len(soft_to_hardclip_pairs)} hardclip-softclip split alignment pairs from {sum(soft_to_hardclip_pairs.values())} reads determining start/end positons.')


    with open(opr_file) as handle:
        for line in handle:
            if line.startswith('#'):
                if line.startswith('#MAX_REASONABLE_INSERTSIZE'):
                    max_reasonable_insert_size = int(line.split('=')[1])
                if line.startswith('#ESTIMATED_INSERTSIZE'):
                    estimated_insert_size = int(line.split('=')[1])
                continue
            splits = line.strip().split()
            if splits[-1] == 'OPR':
                pos1 = int(splits[6])
                pos2 = int(splits[7])
                scaffold = splits[1]
                if pos1 < pos2:
                    oprs_start_to_stop[(pos1, pos2, scaffold)] += 1
                else:
                    oprs_start_to_stop[(pos2, pos1, scaffold)] += 1

    logging.info(f'Adding OPR information from {len(oprs_start_to_stop)} positions and {sum(oprs_start_to_stop.values())} inserts.')

    opr_supported_start_ends = {}
    for (hsp1, hsp2, hsscaffold), hscnt in soft_to_hardclip_pairs.items():
        opr_supported_start_ends[(hsp1, hsp2, hsscaffold)] = (hscnt, 0)

    for (oprp1, oprp2, oprscaffold), oprcnt in oprs_start_to_stop.items():
        distances = {}
        for (hsp1, hsp2, hsscaffold), hscnt in soft_to_hardclip_pairs.items():
            if oprscaffold == hsscaffold:
                deltap1 = abs(hsp1 - oprp1)
                deltap2 = abs(hsp2 - oprp2)
                delta = deltap1 + deltap2
                if delta <= max_reasonable_insert_size:
                    distances[(hsp1, hsp2, hsscaffold, hscnt)] = delta
        winner = (oprp1, oprp2, oprscaffold, 0)
        if len(distances) > 0:
            max_hscnt = max(map(lambda distance: distance[3], distances))
            filtered_distances = {distance:insert_size for (distance,insert_size) in distances.items() if distance[3] == max_hscnt}
            if len(filtered_distances) == 1:
                winner = list(filtered_distances.items())[0][0]
            else:
                min_dev_from_insert_size = min([abs(isize - estimated_insert_size) for isize in filtered_distances.values()])
                double_filtered_distances = {distance: isize for (distance, isize) in filtered_distances.items() if abs(isize - estimated_insert_size) == min_dev_from_insert_size}
                winner = list(double_filtered_distances.items())[0][0]
        else:
            winner = (oprp1, oprp2, oprscaffold, 0)

        entry = opr_supported_start_ends.get((winner[0], winner[1], winner[2]), (0,0))
        entry = (entry[0], entry[1] + oprcnt)
        opr_supported_start_ends[(winner[0], winner[1], winner[2])] = entry




    logging.info(f'Added OPR information. Now working with {len(opr_supported_start_ends)} start/end combinations.')
    #minmvirlength=1000, maxmvirlength=1000000, allow_complete_scaffolds=True
    minsize = minmvirlength
    maxsize = maxmvirlength
    minoprcount = 1
    minhscount = 1
    mincombcount = 5
    


    # for (start, end, scaffold), (hs_cnt, opr_cnt) in sorted(opr_supported_start_ends.items(), key=lambda i: sum(i[1]),
    #                                                         reverse=True)[:10]:
    #     print((start, end, scaffold, end - start, hs_cnt, opr_cnt))
    
    logging.info(f'Filtering: \n\tby length {minsize} <= length <= {maxsize}. \n\t#OPRs >= {minoprcount}\n\t#HARD-SOFT >= {minhscount}\n\t#OPRS + #HARD-SOFT >= {mincombcount}\n\tAllow full scaffolds: {allow_complete_scaffolds}')
    updated_opr_supported_start_ends = {}
    for (start, end, scaffold), (hs_cnt, opr_cnt) in sorted(opr_supported_start_ends.items(), key=lambda i: sum(i[1]), reverse=True):
        # 1:
        length = end - start + 1
        
        if minsize >= length or length >= maxsize:
            continue
        # 2:
        if opr_cnt < minoprcount or hs_cnt < minhscount or (opr_cnt + hs_cnt) < mincombcount:
            continue
        # 3:
        scaffold_length = len(reference_header_2_sequence[scaffold])
        if not allow_complete_scaffolds:
            if scaffold_length * 0.99 <  length:
                continue
        updated_opr_supported_start_ends[(start, end, scaffold)] = (hs_cnt, opr_cnt)
    opr_supported_start_ends = updated_opr_supported_start_ends
    logging.info(f'Finished filtering. Working with {len(opr_supported_start_ends)} start/end combinations.')


    ref_opr_supported_start_ends = {}
    for (start, end, scaffold), (hs_cnt, opr_cnt) in sorted(opr_supported_start_ends.items(), key=lambda i: sum(i[1]),reverse=True):
        found = False
        for (refstart, refend, refscaffold), (ref_hs_cnt, ref_opr_cnt) in sorted(ref_opr_supported_start_ends.items(), key=lambda i: sum(i[1]), reverse=True):
            if refscaffold == scaffold:
                if refstart - softclip_range <= start <= refstart + softclip_range:
                    if refend - softclip_range <= end <= refend + softclip_range:
                        ref_opr_supported_start_ends[(refstart, refend, refscaffold)] = (hs_cnt + ref_hs_cnt, opr_cnt + ref_opr_cnt)
                        found = True
                        break
        if not found:
            ref_opr_supported_start_ends[(start, end, scaffold)] = (hs_cnt, opr_cnt)

    with open(output_fasta_file, 'w') as outhandle:
        for (start, end, scaffold), (hs_cnt, opr_cnt) in sorted(ref_opr_supported_start_ends.items(), key=lambda i: sum(i[1]),reverse=True):
            scaffold_sequence = reference_header_2_sequence[scaffold]
            scaffold_length = len(scaffold_sequence)
            subsequence = scaffold_sequence[start:end+1]
            scaffold_coverage = 100.0 * float(len(subsequence))/float(scaffold_length)
            scaffold_coverage = "{:0.6f}".format(scaffold_coverage)
            outhandle.write(f'>{scaffold}:{start}-{end}\tOPRs={opr_cnt}-HSs={hs_cnt}-SF={scaffold_coverage}\n{subsequence}\n')






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



class CapitalisedHelpFormatter(argparse.HelpFormatter):
    def add_usage(self, usage, actions, groups, prefix=None):
        if prefix is None:
            prefix = ''
        return super(CapitalisedHelpFormatter, self).add_usage(usage, actions, groups, prefix)


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
        -o  PATH   Prefix for output file. [Required]
    
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
    parser.add_argument('-o', action='store', help='Output prefix', required=True, dest='output')
    parser.add_argument('-t', action='store', dest='threads',help='Number of threads to use. (Default = 1)',type=int, default=1)

    try:
        args = parser.parse_args(sys.argv[2:])
    except:
        shutdown(1)



    forward_read_file = args.forward
    reversed_read_file = args.reverse

    out_bam_file = args.output + '.bam'
    bwa_ref_name = args.db
    opr_file = pathlib.Path(args.output + '.oprs')
    clipped_file = pathlib.Path(args.output + '.clipped')
    output_fasta_file = pathlib.Path(args.output + '.fasta')
    minlength_report = args.ml
    maxlength_report = args.ML
    allow_fl_report = args.afs
    threads = args.threads
    _execute_oprs(forward_read_file, reversed_read_file, out_bam_file, bwa_ref_name, opr_file, clipped_file, output_fasta_file, minlength_report, maxlength_report, allow_fl_report, threads)

def _execute_oprs(forward_read_file, reversed_read_file, out_bam_file, bwa_ref_name, opr_file, clipped_file, output_fasta_file, minlength_report=4000, maxlength_report=800000, allow_fl_report=True, threads=1):
    min_percid = 0.97
    remove_unmapped = True
    min_coverage = 0.8
    min_alength = 45


    align_minalength = 0
    align_mincoverage = 0.0





    ### CHECKS START
    if minlength_report < 0:
        raise argparse.ArgumentTypeError('-ml has to be >0'.format(minlength_report))
        shutdown(1)

    if maxlength_report < 0:
        raise argparse.ArgumentTypeError('-ML has to be >0'.format(maxlength_report))
        shutdown(1)
    if threads <= 0:
        raise argparse.ArgumentTypeError('Number of threads has to be >0'.format(threads))
        shutdown(1)



    if not pathlib.Path(forward_read_file).exists() or not pathlib.Path(forward_read_file).is_file():
        logging.error(f'The forward reads file does not exist: {forward_read_file}. Quitting')
        shutdown(1)
    if not pathlib.Path(reversed_read_file).exists() or not pathlib.Path(reversed_read_file).is_file():
        logging.error(f'The reverse reads file does not exist: {reversed_read_file}. Quitting')
        shutdown(1)


    required_index_files = [bwa_ref_name + suffix for suffix in ['.bwt', '.pac', '.ann', '.amb', '.sa']]
    is_reference_missing = False
    for rif in required_index_files:
        if not pathlib.Path(rif).exists() or not pathlib.Path(rif).is_file():
            is_reference_missing = True
    if is_reference_missing:
        logging.error(f'The bwa index files are missing. Please rerun mvirs index. Quitting')
        shutdown(1)

    check_sequences(forward_read_file, reversed_read_file)

    ### CHECKS END


    align(forward_read_file, reversed_read_file, bwa_ref_name, threads, out_bam_file, align_mincoverage, min_percid, align_minalength)
    find_clipped_reads(out_bam_file,clipped_file)
    find_oprs(out_bam_file, opr_file, min_coverage, min_alength)
    extract_regions(clipped_file, opr_file, bwa_ref_name, output_fasta_file,minmvirlength=minlength_report, maxmvirlength=maxlength_report, allow_complete_scaffolds=allow_fl_report)



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
    ''', formatter_class=CapitalisedHelpFormatter,add_help=False)
    parser.add_argument('-f', action='store', help='Input FastA or FastQ file for index building. Gzipped input allowed.', required=True)
    try:
        args = parser.parse_args(sys.argv[2:])
    except:
        shutdown(1)
    seq_file = args.f
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
