from pysam.libcalignmentfile import AlignmentFile
from typing import Dict

from .containers import MappedRead, Mapping

from collections import namedtuple, defaultdict
from itertools import product
import logging

PAlignment = namedtuple('PAlignment',
                        'iss ref revr1 revr2 score startr1 endr1 startr2 endr2 orientation')


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


def mapped_reads_from_alignment(bam_file : AlignmentFile,
                                min_coverage : float = 0,
                                min_alength : int = 0,
                                extended_info=True) -> dict:
    """
    Get all the mapping for all the reads from the BAM file.

    Args:
        bam_file (str): Path to the BAM file.
        min_coverage (float): Minimum alignment coverage of the read
                              to be considered mapped.
        min_alength (int): Minimum alignment length of the read
                           to be considered mapped.
        extended_info (bool): If True, the Mapping object will contain
                              CIGAR and gapless blocks information.

    Returns:
        dict: A dictionary of MappedRead objects.
    """
    mapped: dict = {}
    alignments: AlignmentFile = AlignmentFile(str(bam_file), "rb")

    for alignment in alignments:

        # if alignment not eligible
        # dismiss straight away
        if not (alignment.get_tag('al') >= min_alength and alignment.get_tag('qc') >= min_coverage):
            continue

        data_tmp = alignment.qname.rsplit('/', 1)

        readname = data_tmp[0]
        orientation = data_tmp[1]

        if readname not in mapped:
            mapped[readname] = MappedRead(name=readname)
        mapped[readname][orientation] = Mapping.from_aligned_segment(alignment, extended=extended_info)

    alignments.close()

    return mapped


def _calculate_orientation(rev_r1: bool, pos_r1: int,
                           rev_r2: bool, pos_r2: int) -> str:

    """
    Based on orientation of both reads and their alignment start positions,
    estimate if a read is PE/SAME/OPR.

    Args:
        rev_r1 (bool): True if the read is mapped in reverse, False otherwise.
        pos_r1 (int): The start position of the read on the reference.
        rev_r2 (bool): True if the read is mapped in reverse, False otherwise.
        pos_r2 (int): The start position of the read on the reference.

    Returns:
        str: PAIREDEEND/SAME/OPR.
    """

    orientation = "PAIREDEND"
    if rev_r1 == rev_r2:
        orientation = "SAME"
    else:
        if rev_r1:
            rev_pos = pos_r1
            fwd_pos = pos_r2
        else:
            rev_pos = pos_r2
            fwd_pos = pos_r1

        if rev_pos < fwd_pos:
            orientation = "OPR"

    return orientation


def get_paired_alignments(read_mappings: Dict[str, MappedRead],
                          cutoff_bestscore: float = 0.95):

    cnt_singleend = 0
    cnt_pairedend = 0
    score_filtered = {}

    for query, alignments in read_mappings.items():
        if alignments.is_single_end():
            cnt_singleend += 1
            continue
        else:
            cnt_pairedend += 1

            samematches = []
            bestscore = 0
            matches = defaultdict(list)
            cutoff_bestscore = 0.95

            for alnr1, alnr2 in product(alignments['R1'],
                                        alignments['R2']):
                if alnr1.ref != alnr2.ref:  # if both reads align to 2 different references we continue
                    continue
                orientation: str = _calculate_orientation(alnr1.rev, alnr1.rstart, alnr2.rev,  alnr2.rstart)
                positions = [alnr1.rstart, alnr2.rstart, alnr1.rend, alnr2.rend]
                iss = abs(max(positions) - min(positions))
                score = alnr1.score + alnr2.score
                # add SAME oritentation to a separate container
                if orientation == 'SAME':
                    samematches.append(
                        PAlignment(iss=iss, revr1=alnr1.rev, revr2=alnr2.rev,
                                score=score, startr1=alnr1.rstart, endr1=alnr1.rend, startr2=alnr2.rstart, endr2=alnr2.rend,
                                orientation=orientation, ref=alnr1.ref))

                else:
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
                    break
                match = matches[score]
                tmpmatch = sorted(match, key=lambda x: x.iss)
                scoreinsertsizesortedmatches.extend(tmpmatch)

            for samematch in samematches:
                if samematch.score >= minscore:
                    scoreinsertsizesortedmatches.append(samematch)

            if not scoreinsertsizesortedmatches:
                continue

            score_filtered[query] = scoreinsertsizesortedmatches

    return score_filtered

def extract_insert_sizes(paired_alignments, n=20_000):
     # only use first N values
    insert_sizes = []
    for pe_alignments in paired_alignments.values():
        if len(pe_alignments) == 1 and pe_alignments[0].orientation == "PAIREDEND":
            insert_sizes.append(pe_alignments[0].iss)
            if len(insert_sizes) == n:
                break

    return insert_sizes

def _estimate_insert_size(insert_sizes, sd = 7):

    def get_mean(array: list):
        return int(sum(array) / len(array))

    def median(array: list):
        return int(sorted(array)[len(array) // 2])

    def pstdev(array: list):
        m = get_mean(array)
        numerator = sum([(x - m) ** 2 for x in array])
        denominator = len(array)
        return int((numerator / denominator) ** 0.5)

    logging.info(
        'Estimating insert size by mean/stdev convergence.'
    )
    min_val_tmp = -1
    max_val_tmp = -1
    iteration = 0
    while True:
        iteration += 1
        mean = get_mean(insert_sizes)
        stdev = pstdev(insert_sizes)

        # For mapping Illumina short-insert reads to the human genome,
        # x is about 6-7 sigma away from the mean.
        # http://bio-bwa.sourceforge.net/bwa.shtml section Estimating Insert Size Distribution
        min_val = mean - sd * stdev
        min_val = max(0, min_val)  # Ensure min_val is not negative
        max_val = mean + sd * stdev

        filtered_sizes = [val for val in insert_sizes if min_val <= val <= max_val]
        insert_sizes = filtered_sizes

        logging.info('Iteration %d. Min insert size: %d. Max insert size: %d. '
            'Mean insert size: %d. Stdev insert size: %d' % (iteration, min_val, max_val, mean, stdev))


        if min_val == min_val_tmp and max_val == max_val_tmp:
            break

        min_val_tmp = min_val
        max_val_tmp = max_val

    if min_val != 0:
        min_val = 0

    return min_val, max_val, median(insert_sizes)


# sizes = extract_insert_sizes(paired)
# _estimate_insert_size([1, 2, 10, 11, 11, 12, 14, 18], sd=2)


## TODO: write tests
## TODO: throw away found == False alignments
def _calc_primary_paired_alignment(insert2alignments, min_insert_size: int, max_insert_size: int):
    """
    Go through all paired alignments and find one that has a reasonable good score and is close to the insert size.
    1. if there is a PE alignment with reasonable insert size --> take that one (and all other equally good PE with reasonable IS)
    2. Otherwise take the best one
    3. If there are multiple best ones return all and state add the false flag to indicate that alignment could not be resolved to a single best alignment
    :param insert2alignments:
    :param insertsize:
    :return:
    """

    primary = []

    for qname, alignments in insert2alignments.items():

        if len(alignments) == 1:  # This is a unique alignment. No filtering required.
            primary.append((qname, alignments, True))
        else:
            # Check if there are any alignments with PE and reasonable insert size
            pe_candidates = [aln for aln in alignments if aln.orientation == 'PAIREDEND' and
                             min_insert_size <= aln.iss <= max_insert_size]
            if pe_candidates:
                 # if there are alignments with PE and reasonable insertsize
                 # then pick the ones with the best score
                bestscore = max(pe_candidates, key=lambda x: x.score).score
                pe_bestalns = [aln for aln in pe_candidates if aln.score == bestscore]
                if pe_bestalns:
                    primary.append((qname, pe_bestalns, True))

            else:
                primary.append((qname, alignments, False))

    return primary


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

    mapped = mapped_reads_from_alignment(out_bam_file, min_coverage=min_coverage,
                                         min_alength=min_alength, extended_info=False)
    paired = get_paired_alignments(mapped)

    insert_sizes = extract_insert_sizes(paired, n=20_000)

    logging.info('Start estimating insert size from paired end alignments.')
    min_insertsize, max_insertsize, median_insertsize = _estimate_insert_size(insert_sizes)
    logging.info('Finished estimating insert size from paired end alignments.')
    logging.info('Min reasonable insert size:\t{}'.format(min_insertsize))
    logging.info('Max reasonable insert size:\t{}'.format(max_insertsize))

    primary = _calc_primary_paired_alignment(paired, min_insert_size=min_insertsize,
                                             max_insert_size=max_insertsize)

    with open(opr_file, 'w') as handle:
        handle.write(
            '#MIN_REASONABLE_INSERTSIZE={}\n#MAX_REASONABLE_INSERTSIZE={}\n#ESTIMATED_INSERTSIZE={}\n#READNAME\tREFERENCE\tINSERT_SIZE\tR1_ORIENTATION\tR2_ORIENTATION\tBWA_SCORE\tR1_START\tR2_START\tR1_ALNLENGTH\tR2_ALNLENGTH\tREAD_ORIENTATION\n'.format(
                min_insertsize, max_insertsize, median_insertsize))



        template: str = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'

        logging.info('Start screening for OPRs and Paired-End inserts with unreasonable insert size.')
        alncnt = 0
        for alncnt, (query, alignments, found) in enumerate(primary):
            if found:
                for alignment in alignments:
                    printstring = template.format(query, alignment.ref, alignment.iss,
                                                  get_read_orientation(alignment.revr1),
                                                  get_read_orientation(alignment.revr2), alignment.score, alignment.startr1,
                                                  alignment.startr2, alignment.endr1 - alignment.startr1, alignment.endr2 - alignment.startr2,
                                                  alignment.orientation)
                    if alignment.orientation == 'PAIREDEND' and min_insertsize <= alignment.iss <= max_insertsize:  # pe with reasonable insert size
                        pass
                    else:  # pe with unreasonable insert size and opr/same
                        printstring = printstring.replace('PAIREDEND', 'IPR')
                        handle.write(printstring)
            else:  # undefined
                continue

            if alncnt % 500000 == 0:
                logging.info('Paired inserts processed:\t{}'.format(format(alncnt, ',d')))
        logging.info('Paired inserts processed:\t{}'.format(format(alncnt, ',d')))
        logging.info('Finished screening for OPRs and Paired-End inserts with unreasonable insert size.')