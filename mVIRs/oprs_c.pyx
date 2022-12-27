from collections import defaultdict, namedtuple
import logging
import statistics
import sys 

from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from libcpp.string cimport string 
from libcpp.vector cimport vector

PAlignment = namedtuple('PAlignment',
                        'iss, ref revr1 revr2 score startr1 endr1 startr2 endr2 orientation')
SAMLine = namedtuple('SAMLine', 'rev ref rstart rend score cigartuples blocks')

# cdef class PAlignment: 
#     cdef readonly int iss
#     cdef readonly str ref
#     cdef readonly bint revr1 
#     cdef readonly bint revr2 
#     cdef readonly int score 
#     cdef readonly int startr1
#     cdef readonly int endr1  
#     cdef readonly int startr2
#     cdef readonly int endr2 
#     cdef readonly str orientation

#     def __init__(self, iss, ref, revr1, revr2, score, startr1, endr1, 
#                  startr2, endr2, orientation):
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


# cdef class SAMLine:
#     cdef readonly bint rev
#     cdef readonly str ref
#     cdef readonly int rstart
#     cdef readonly int rend
#     cdef readonly int score
#     cdef readonly list blocks 
#     cdef readonly list cigartuples

#     def __init__(self, rev, ref, rstart, rend, score, blocks, cigartuples):
#         self.rev = rev
#         self.ref = ref
#         self.rstart = rstart
#         self.rend = rend
#         self.score = score
#         self.blocks = blocks
#         self.cigartuples = cigartuples    


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
    cdef: 
        str orientation = 'PAIREDEND'
        int fwpos = 0
        int revpos = 0
    
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
    
    cdef: 
        list data_tmp = alignment.qname.rsplit('/', 1)
        str readname = data_tmp[0]
        str orientation = data_tmp[1]
        int ascore = alignment.get_tag('AS')
        bint reverse = True if alignment.is_reverse else False
        str refname = alignment.reference_name
        int refstart = alignment.reference_start
        int refend = alignment.reference_end
        list blocks = None
        list cigartuples = None

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
                                       score=alnr1.score + alnr2.score, 
                                       startr1=alnr1.rstart, endr1=alnr1.rend, 
                                       startr2=alnr2.rstart, endr2=alnr2.rend,
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
        tmp = []
        iteration += 1
        logging.info(
            'Estimating insert size by mean/stdev convergence. Iteration:\t{}. Min insert size:\t{}. Max insert size:\t{}. Mean insert size:\t{}. Stdev insert size:\t{}'.format(
                iteration, format(min(so), ',d'), format(max(so), ',d'), format(int(statistics.mean(so)), ',d'), format(int(statistics.stdev(so)), ',d')))
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

        insert2alignments_gen = insertize_bamfile_by_name(out_bam_file, max_sam_lines, min_coverage,
                                                          min_alength, need_extended=False)
        logging.info('Start screening for OPRs and Paired-End inserts with unreasonable insert size.')
        alncnt = 0
        for alncnt, (query, alignments, found) in enumerate(
                _calc_primary_paired_alignment(insert2alignments_gen, estimated_insertsize, 
                                               minreasonable_insertsize,
                                               maxreasonable_insertsize)):
            if found:
                for cnt, alignment in enumerate(alignments):
                    printstring = template.format(query, alignment.ref, alignment.iss,
                                                  get_read_orientation(alignment.revr1),
                                                  get_read_orientation(alignment.revr2), 
                                                  alignment.score, 
                                                  alignment.startr1,
                                                  alignment.startr2, 
                                                  alignment.endr1 - alignment.startr1, 
                                                  alignment.endr2 - alignment.startr2,
                                                  alignment.orientation)
                    if alignment.orientation == 'PAIREDEND' and alignment.iss >= minreasonable_insertsize \
                    and alignment.iss <= maxreasonable_insertsize:  # pe with reasonable insert size
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