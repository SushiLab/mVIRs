import pysam
from pysam.libcalignmentfile import AlignmentFile
from typing import Dict

from .containers import MappedRead, Mapping

from collections import namedtuple, defaultdict
from itertools import product
import logging
from typing import Dict, List

PAlignment = namedtuple('PAlignment',
                        'iss ref revr1 revr2 score startr1 endr1 startr2 endr2 orientation')


def get_read_orientation(rev: bool) -> str:
    """
    Helper function to print if an alignments is forward or reverse
    """

    if rev:
        return 'reverse'
    else:
        return 'forward'


def mapped_reads_from_alignment(bam_file : AlignmentFile,
                                min_coverage : float = 0,
                                min_alength : int = 0,
                                max_lines : int = None,
                                extended_info = True,
                                threads: int = 1) -> Dict[str, MappedRead]:
    """
    Get all the mapping for all the reads from the BAM file.

    Args:
        bam_file: Path to the BAM file.
        min_coverage: Minimum alignment coverage of the read
          to be considered mapped.
        min_alength: Minimum alignment length of the read
          to be considered mapped.
        extended_info: If True, the Mapping object will contain
          CIGAR and gapless blocks information.
        threads: Number of threads to use.

    Returns:
        dict: A dictionary of MappedRead objects.
    """
    mapped = {}

    # https://github.com/pysam-developers/pysam/issues/939
    save = pysam.set_verbosity(0)
    alignments = AlignmentFile(str(bam_file), "rb", threads=threads)
    pysam.set_verbosity(save)


    for i, alignment in enumerate(alignments):

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
        if max_lines and i >= max_lines:
            break

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
                          cutoff_bestscore: float = 0.95
                          ) -> Dict[str, List[PAlignment]]:
    """
    Get the paired alignments from the read mappings.

    Args:
        read_mappings (dict): A dictionary of MappedRead objects.
        cutoff_bestscore (float): The cutoff to use for filtering alignments.

    Returns:
        dict: A dictionary of paired alignments.
    """


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

def extract_insert_sizes(paired_alignments: Dict[str, List[PAlignment]],
                         n: int = 20_000):
    """
    Extract the insert sizes from the N first paired alignments
    for insert size estimation.

    Args:
        paired_alignments (dict): A dictionary of PAlignment objects.
        n (int): The number of alignments to use for insert size estimation.
    """

     # only use first N values
    insert_sizes = []
    for pe_alignments in paired_alignments.values():
        if len(pe_alignments) == 1 and pe_alignments[0].orientation == "PAIREDEND":
            insert_sizes.append(pe_alignments[0].iss)
            if len(insert_sizes) == n:
                break

    return insert_sizes

def _estimate_insert_size(insert_sizes: list, sd : int = 7):
    """
    Estimate the insert size from the paired alignments.
    This is done by calculating the mean and standard deviation of the insert sizes until convergence.

    Args:
        insert_sizes (list): A list of insert sizes.
        sd (int): The standard deviation to use for filtering.

        Returns:
            tuple: A tuple of three elements. The first element is the minimum insert size. The second element is the
                maximum insert size. The third element is the median insert size."""


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
    1. Read the original bam file into memory and pair.
    2. Estimate the insert size from the unique paired mappers.
    We do that to get a reasonable insert size for filtering
    potentially good (paired-end) or strange (OPR) alignments.
    3. Calculate the best primary alignments.
    3.1  if there is a unique best alignment --> take that one
    3.2. elif there is paired-end alignment with reasonable insert size and within 5% of the best score --> take that one
    3.3. elif there is an OPR that is the best alignment for this insert --> take that one
    3.4  elif take every best mapper
    The output is a file that reports for each insert for which we can identify a best mapper and has either an:
    - an unreasonable insert size
    - an abnormal orientation (OPR <---- ---->) or (SAME ---> --->)

    Args:
        out_bam_file (str): Path to the BAM file.
        opr_file (str): Path to the output file.
        min_coverage (float): Minimum alignment coverage of the read.
            to be considered mapped.
        min_alength (int): Minimum alignment length of the read.

    Returns:
        None
    """

    logging.info('Start OPR finding step')
    logging.info('Input BAM File:\t%s', out_bam_file)
    logging.info('Output OPR File:\t%s', opr_file)

    mapped = mapped_reads_from_alignment(out_bam_file, min_coverage=min_coverage,
                                         min_alength=min_alength, extended_info=False)
    paired = get_paired_alignments(mapped)

    insert_sizes = extract_insert_sizes(paired, n=20_000)

    logging.info('Start estimating insert size from paired end alignments.')
    min_insertsize, max_insertsize, median_insertsize = _estimate_insert_size(insert_sizes)
    logging.info('Finished estimating insert size from paired end alignments.')
    logging.info('Min reasonable insert size:\t%s', min_insertsize)
    logging.info('Max reasonable insert size:\t%s', max_insertsize)

    primary = _calc_primary_paired_alignment(paired, min_insert_size=min_insertsize,
                                             max_insert_size=max_insertsize)

    with open(opr_file, 'w') as handle:
        handle.write(
            '#MIN_REASONABLE_INSERTSIZE={}\n#MAX_REASONABLE_INSERTSIZE={}\n#ESTIMATED_INSERTSIZE={}' \
            '\n#READNAME\tREFERENCE\tINSERT_SIZE\tR1_ORIENTATION\tR2_ORIENTATION\tBWA_SCORE\tR1_START'\
            '\tR2_START\tR1_ALNLENGTH\tR2_ALNLENGTH\tREAD_ORIENTATION\n'.format(
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
                logging.info('Paired inserts processed:\t%s', format(alncnt, ',d'))
        logging.info('Paired inserts processed:\t%s', format(alncnt, ',d'))
        logging.info('Finished screening for OPRs and Paired-End inserts with unreasonable insert size.')



# logging.info('Start clipped reads finding step')
# logging.info('Input BAM File:\t{}'.format(bam_file))
# logging.info('Output Clipped File:\t{}'.format(clipped_file))
def check_clipping(alignment: Mapping, clip_type: int):
    """
    Check if the alignment is clipped at the start and/or end.

    Args:
        alignment (Mapping): The alignment to check.
        clip_type (int): The type of clipping to check. 4 for softclipping, 5 for hardclipping
        (from pysam documentation).

    Returns:
        list: A list of two elements. The first element is 1 if the alignment is clipped at the start,
            0 otherwise. The second element is 1 if the alignment is clipped at the end, 0 otherwise.

    """

    clip_location = [0, 0]
    cigar_tuples = alignment.cigartuples

    if cigar_tuples[0][0] == clip_type:
        clip_location[0] = 1
    if cigar_tuples[-1][0] == clip_type:
        clip_location[1] = 1

    return clip_location


def find_clipped_reads(paired_alignments: Dict[str, MappedRead],
                       clipped_file: str) -> None:
    """
    Find clipped reads in a BAM file. The goal of the clipped reads is to find start/end positions of
    activated prophages. This can be done by looking at cigar string that contain a S (softclips).
    If the S is before a M then this is a START position. If the S is after an M then that is am=n END
    position.
    For convenience, CIGAR tuples instead of CIGAR strings are used.

    Args:
        paired_alignments (dict): A dictionary of MappedRead objects.
        clipped_file (str): Path to the output file.

    Returns:
        None
    """

    # pysam operations
    PYSAM_BAM_CSOFT_CLIP = 4
    PYSAM_BAM_CHARD_CLIP = 5
    # counters
    alignments_seen = 0
    alignments_softclipped = 0
    alignments_hardclipped = 0


    clipped_reads = defaultdict(list)
    for insert, ori_2_alignments in paired_alignments.items():
        for ori, alignments in ori_2_alignments.mapping.items():
            for alignment in alignments:
                alignments_seen += 1
                if alignments_seen % 500000 == 0:
                    logging.info('%s - %s - %s | Alignments - Softclips - Hardclips',
                                 format(alignments_seen, ",d"), format(alignments_softclipped, ",d"),
                                 format(alignments_hardclipped, ",d"))
                if len(alignment.cigartuples) > 1:
                    # 1st number in cigartuple denotes the operation
                    cigar_operations = set(map(lambda cigar: cigar[0], alignment.cigartuples))
                    # 4 is code for softclip
                    softclipped = PYSAM_BAM_CSOFT_CLIP in cigar_operations
                    # 5 is code for hardclip
                    hardclipped = PYSAM_BAM_CHARD_CLIP in cigar_operations

                    if softclipped:
                        clip_location = check_clipping(alignment, PYSAM_BAM_CSOFT_CLIP)
                        if clip_location == [1, 1]:
                            continue
                        alignments_softclipped += 1
                        clip_type = "S"

                    elif hardclipped and not softclipped:
                        clip_location = check_clipping(alignment, PYSAM_BAM_CHARD_CLIP)
                        if clip_location == [1, 1]:
                            continue
                        alignments_hardclipped += 1
                        clip_type = "H"

                    else:
                        continue

                    direction = '->' if clip_location[0] == 1 else '<-'
                    startpos = alignment.blocks[0][0] if clip_location[0] == 1 else alignment.blocks[-1][1]
                    clipped_reads[(insert, ori)].append((clip_type, direction, startpos, alignment.ref))

    logging.info('%s - %s - %s | Alignments - Softclips - Hardclips',
                 format(alignments_seen, ",d"), format(alignments_softclipped, ",d"),
                 format(alignments_hardclipped, ",d"))
    logging.info('Keeping only paired soft/hard-clips and writing unfiltered soft and filtered hard-clips to %s.', clipped_file)

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
        softclip_ori, softclip_coord, softclip_ref = softclipped[0][1:4]
        hardclipped_filtered = []
        for hardclip in hardclipped:
            hardclip_ori, hardclip_coordinate, hardclip_ref = hardclip[1:4]

            if hardclip_ref != softclip_ref:
                continue
            if softclip_ori == hardclip_ori:
                continue
            if softclip_ori == '->' and softclip_coord > hardclip_coordinate:
                continue
            else:
                hardclipped_filtered.append(hardclip)

        hardclips_written = hardclips_written + len(hardclipped_filtered)
        softclips_written = softclips_written + len(softclipped)
        for aln in (softclipped + hardclipped_filtered):
            out_file.write(f'{insert}\t{ori}\t{aln[0]}\t{aln[1]}\t{aln[2]}\t{aln[3]}\n')

    out_file.close()
    logging.info(f'Wrote {format(softclips_written, ",d")} soft and {format(hardclips_written, ",d")} hard-clips.')