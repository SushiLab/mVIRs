import collections
import logging
import sys
import pysam
from typing import Dict, Generator, List, Tuple, DefaultDict
import statistics

debug = False

PYSAM_BAM_CSOFT_CLIP = 4
PYSAM_BAM_CHARD_CLIP = 5

PAlignment = collections.namedtuple('PAlignment',
                                    'iss ref revr1 revr2 score startr1 endr1 startr2 endr2 orientation')
SAMLine = collections.namedtuple('SAMLine', 'rev ref rstart rend score cigartuples blocks')

# immutable class
# class PAlignment:

#     def __init__(self, iss: bool, ref: str, revr1: bool, revr2: bool, score: int, startr1: int, endr1: int, startr2: int, endr2: int, orientation: str):
#         self.iss = iss
#         self.ref = ref
#         self.revr1 = revr1
#         self.revr2 = revr2
#         self.score = score
#         self.startr1 = startr1
#         self.endr1 = endr1
#         self.startr2 = startr2
#         self.endr2 = endr2
#         self.orientation = orientation


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

def insertize_bamfile_by_name(bam_file, max_sam_lines: int = -1, min_coverage: float = 0.0, min_alength: int = 0, need_extended = True):
    save = pysam.set_verbosity(0)
    alignments = pysam.AlignmentFile(str(bam_file), "rb")
    pysam.set_verbosity(save)

    current_name = None
    current_insert = collections.defaultdict(list)

    for alignment in alignments:
        # filter first, avoid unnecessary operations
        if not (alignment.get_tag('al') >= min_alength and alignment.get_tag('qc') >= min_coverage):
            continue
        # store data in data structures
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


def _estimate_insert_size(insert2alignments: Dict[str, Dict[str, List[SAMLine]]]) -> Tuple[int, int, int]:
    """
    Goes through UNIQUE PAIRED END alignments and estimates the median insert size.
    :param insert2alignments:
    :return:
    """
    insertsizes: List[int] = []

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