cimport cython
from collections import Counter, defaultdict, namedtuple
import logging
import statistics
import sys

from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment


PAlignment = namedtuple('PAlignment',
                        'iss, ref revr1 revr2 score startr1 endr1 startr2 endr2 orientation')
SAMLine = namedtuple('SAMLine', 'rev ref rstart rend score cigartuples blocks')


cdef inline str get_read_orientation(bint rev):
    cdef str orientation
    if rev:
        orientation = "reverse"
    else:
        orientation = "forward" 
    
    return orientation


cdef inline str _calc_orientation(bint revr1, bint revr2, int posr1, int posr2):
    """
    Based on orientation of both reads and their alignment start positions, estimate if a read is PE/SAME/OPR
    :param revr1:
    :param revr2:
    :param posr1:
    :param posr2:
    :return: PAIREDEND, SAME, OPR
    """
    cdef str orientation = 'PAIREDEND'
    cdef int fwpos = 0
    cdef int revpos = 0
    
    if revr1 == revr2:
        orientation = 'SAME'
        return orientation
 
    if revr1:
        revpos = posr1
        fwpos = posr2
    else:
        revpos = posr2
        fwpos = posr1
    if revpos < fwpos:
        orientation = 'OPR'
    return orientation

cdef inline tuple extract_alignment_info(AlignedSegment alignment, bint need_extended):
    
    cdef list data_tmp = alignment.qname.rsplit('/', 1)
    cdef str readname = data_tmp[0]
    cdef str orientation = data_tmp[1]
    cdef int ascore = alignment.get_tag('AS')
    cdef bint reverse = True if alignment.is_reverse else False
    cdef str refname = alignment.reference_name
    cdef int refstart = alignment.reference_start
    cdef int refend = alignment.reference_end
    cdef list blocks = None
    cdef list cigartuples = None

    if need_extended:
        blocks = alignment.get_blocks()
        cigartuples = alignment.cigartuples
    samline = SAMLine(rev=reverse, ref=refname, rstart=refstart, 
                      rend=refend, score=ascore, cigartuples=cigartuples, blocks=blocks)
    
    return readname, orientation, samline


def insertize_bamfile_by_name(str bam_file,
                              double max_sam_lines = -1, 
                              double min_coverage = 0.0, 
                              int min_alength = 0, 
                              bint need_extended = 1):
    
    cdef AlignmentFile alignments = AlignmentFile(str(bam_file), "rb")  
    cdef AlignedSegment aligment   
    cdef str current_name = None
    cdef str readname, orientation
    current_insert = defaultdict(list)

    for alignment in alignments:
        if not (alignment.get_tag('al') >= min_alength and alignment.get_tag('qc') >= min_coverage):
            continue

        readname, orientation, samline = extract_alignment_info(alignment, True)

        if not current_name:
            current_name = readname
        if readname == current_name:
            current_insert[orientation].append(samline)
        else:
            yield current_name, current_insert
            current_name = readname
            current_insert: defaultdict = defaultdict(list)
            current_insert[orientation].append(samline)

    if len(current_insert) != 0:
        yield current_name, current_insert


def _generate_paired_alignments(insert2alignments,
                                float cutoff_bestscore = 0.95):
    """
    Iterate over all inserts and generate pairs with a cross product. Pairs with best score are returned.
    Careful --> throws away single end mappers. Considers only paired reads
    :param insert2alignments:
    :return:
    """
    cdef int cnt_singleend = 0
    cdef int cnt_pairedend = 0
    cdef str query
    cdef list samematches, positions
    cdef int bestscore, score, minscore, isp, iss
    cdef str orientation

    for query, alignments in insert2alignments:
        if len(alignments) < 2:
            cnt_singleend += 1
            continue
        else:
            cnt_pairedend += 1
            matches = defaultdict(list)
            samematches = []
            bestscore = 0


            for alnr1 in alignments['R1']:
                for alnr2 in alignments['R2']:
                    if alnr1.ref != alnr2.ref:  # if both reads align to 2 different references we continue
                        continue
                    orientation: str = _calc_orientation(alnr1.rev, alnr2.rev, alnr1.rstart, alnr2.rstart)
                    isp: int = abs(alnr1.rstart - alnr2.rstart)
                    positions: list = [alnr1.rstart, alnr2.rstart, alnr1.rend, alnr2.rend]
                    iss: int = abs(max(positions) - min(positions))
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


            minscore = int(bestscore * cutoff_bestscore)
            scoreinsertsizesortedmatches = []
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


cdef (int, int, int) _estimate_insert_size(insert2alignments):
    """
    Goes through UNIQUE PAIRED END alignments and estimates the median insert size.
    :param insert2alignments:
    :return:
    """

    cdef list insertsizes = []
    cdef list alignments = []

    for (_, alignments) in _generate_paired_alignments(insert2alignments):

        if len(alignments) == 1 and alignments[0].orientation == 'PAIREDEND':
            insertsizes.append(alignments[0].iss)
            if len(insertsizes) > 20000:
                break

    cdef list so = insertsizes

    cdef int minvaltmp = -1
    cdef int maxvaltmp = -1
    cdef int iteration = 0
    cdef list tmp
    cdef int mean, stdev, minval, maxval, val
    
    while True:
        iteration += 1
        logging.info(
            'Estimating insert size by mean/stdev convergence. Iteration:\t{}. Min insert size:\t{}. Max insert size:\t{}. Mean insert size:\t{}. Stdev insert size:\t{}'.format(
                iteration, format(min(so), ',d'), format(max(so), ',d'), format(int(statistics.mean(so)), ',d'), format(int(statistics.stdev(so)), ',d')))
        tmp = []
        mean = int(statistics.mean(so))
        stdev = int(statistics.pstdev(so))
        minval = mean - 7 * stdev
        if minval < 0:
            minval = 0
        maxval = mean + 7 * stdev
        for val in so:
            if minval > val:
                continue
            if maxval < val:
                continue
            tmp.append(val)
        so = tmp

        if minvaltmp == minval and maxval == maxvaltmp:
            break

        minvaltmp = minval
        maxvaltmp = maxval


    # For mapping Illumina short-insert reads to the human genome, x is about 6-7 sigma away from the mean
    # http://bio-bwa.sourceforge.net/bwa.shtml section Estimating Insert Size Distribution

    minvaltmp = mean - 7 * int(statistics.pstdev(so))

    if minvaltmp != 0:
        minvaltmp = 0

    maxvaltmp = mean + 7 * int(statistics.pstdev(so))
    logging.info(
        'Estimating insert size by mean/stdev convergence. Iteration:\t{}. Min insert size:\t{}. Max insert size:\t{}. Mean insert size:\t{}. Stdev insert size:\t{}'.format(
            iteration, format(minvaltmp, ',d'), format(maxvaltmp, ',d'), format(int(statistics.mean(so)), ',d'),
            format(int(statistics.stdev(so)), ',d')))

    return minvaltmp, maxvaltmp, int(statistics.median(so))

cdef inline int get_score(pa_alignment):
    return pa_alignment.score

def _calc_primary_paired_alignment(insert2alignments, 
                                   int insertsize,
                                   int minreasonable_insertsize, 
                                   int maxreasonable_insertsize):
    """
    Go through all paired alignments and find one that has a reasonable good score and is close to the insert size.
    1. if there is a PE alignment with reasonable insert size --> take that one (and all other equally good PE with reasonable IS)
    2. Otherwise take the best one
    3. If there are multiple best ones return all and state add the false flag to indicate that alignment could not be resolved to a single best alignment
    :param insert2alignments:
    :param insertsize:
    :return:
    """

    cdef str qname
    cdef list alignments
    cdef list pa_bestalns, pa_candidates

    for (qname, alignments) in _generate_paired_alignments(insert2alignments):
        if len(alignments) == 1:  # This is a unique alignment. No filtering required.
            yield (qname, alignments, True)
        else:
            # Check if there are any alignments with PE and reasonable insert size
            pe_candidates = []

            for alignment in alignments:
                if alignment.orientation == 'PAIREDEND' and alignment.iss >= minreasonable_insertsize and alignment.iss <= maxreasonable_insertsize:
                    pe_candidates.append(alignment)
            if len(pe_candidates) != 0:
                # if there are alignments with PE and reasonable insertsize then pick the one(s) with the best score
                bestscore = -1
                pe_bestalns = []

                for alignment in sorted(pe_candidates, key=get_score, reverse=True):
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

cdef (int, int) _write_clipped(clipped_reads,
                               str output_file):
    
    cdef str insert, ori

    cdef list softclipped, hardclipped
    cdef str softclipped_ori, softclipped_ref
    cdef int softclipped_coordinate
    
    cdef list hardclipped_filtered
    cdef str hardclipped_ori, hardclipped_ref
    cdef int hardclipped_coordinate

    cdef int hardclips_written = 0
    cdef int softclips_written = 0

    cdef list clipped_alignments
    cdef tuple aln

    out_file = open(output_file, 'wb')

    out_file.write(str.encode('#INSERT\tREADORIENTATION\tHARD/SOFTCLIP'
                              '\tDIRECTION\tPOSITION\tSCAFFOLD\n'))
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
            out_file.write(str.encode(f'{insert}\t{ori}\t{aln[0]}\t'
                                      f'{aln[1]}\t{aln[2]}\t{aln[3]}\n'))

    out_file.close()

    return softclips_written, hardclips_written



def find_clipped_reads(str bam_file, str clipped_file):
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

    cdef bint debug = False

    cdef int PYSAM_BAM_CSOFT_CLIP = 4
    cdef int PYSAM_BAM_CHARD_CLIP = 5
    # START DEBUG PARAMETERS
    if debug:
        max_sam_lines = 10000
    else:
        max_sam_lines = sys.maxsize


    # END DEBUG PARAMETERS
    cdef double min_coverage = 0.0
    cdef int min_alength = 0
    cdef bint toss_singletons = False
    cdef int alignments_seen = 0
    cdef int alignments_softclipped = 0
    cdef int alignments_hardclipped = 0
    cdef list clip_location
    cdef bint front, end 
    cdef str direction
    cdef bint softclipped, hardclipped

    insert2alignments_generator = insertize_bamfile_by_name(bam_file, max_sam_lines, 
                                                            min_coverage, min_alength, 
                                                            need_extended=True)

    clipped_reads = defaultdict(list)
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
    
    softclips_written, hardclips_written = _write_clipped(clipped_reads, clipped_file)

    logging.info(f'Wrote {format(softclips_written, ",d")} soft and {format(hardclips_written, ",d")} hard-clips.')


cpdef void find_oprs(str out_bam_file, 
                     str opr_file,
                     double min_coverage, 
                     int min_alength):
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

    cdef (int, int, int) insertsize_tuple
    cdef int minreasonable_insertsize, maxreasonable_insertsize, estimated_insertsize
    cdef int alncnt, cnt 
    cdef str template = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'
    cdef str printstring
    cdef bint debug, found
    
    # START DEBUG PARAMETERS
    if debug:
        max_sam_lines = 10000
    else:
        max_sam_lines = sys.maxsize
    # END DEBUG PARAMETERS


    insert2alignments = None
    logging.info('Start estimating insert size from paired end alignments.')
    insertsize_tuple = _estimate_insert_size(insertize_bamfile_by_name(out_bam_file, max_sam_lines, 
                                                                       min_coverage, min_alength, need_extended=False))
    minreasonable_insertsize, maxreasonable_insertsize, estimated_insertsize = insertsize_tuple
    logging.info('Finished estimating insert size from paired end alignments.')
    logging.info('Min reasonable insert size:\t{}'.format(minreasonable_insertsize))
    logging.info('Max reasonable insert size:\t{}'.format(maxreasonable_insertsize))

    with open(opr_file, 'wb') as handle:
        header = (f'#MIN_REASONABLE_INSERTSIZE={minreasonable_insertsize}\n'
                  f'#MAX_REASONABLE_INSERTSIZE={maxreasonable_insertsize}\n'
                  f'#ESTIMATED_INSERTSIZE={estimated_insertsize}\n'
                   '#READNAME\tREFERENCE\tINSERT_SIZE\t'
                   'R1_ORIENTATION\tR2_ORIENTATION\tBWA_SCORE\tR1_START\tR2_START\t'
                   'R1_ALNLENGTH\tR2_ALNLENGTH\tREAD_ORIENTATION\n')
        handle.write(str.encode(header))

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
                        handle.write(str.encode(printstring))
            else:  # undefined
                continue

            if alncnt % 500000 == 0:
            
                logging.info('Paired inserts processed:\t{}'.format(format(alncnt, ',d')))
        logging.info('Paired inserts processed:\t{}'.format(format(alncnt, ',d')))
        logging.info('Finished screening for OPRs and Paired-End inserts with unreasonable insert size.')
