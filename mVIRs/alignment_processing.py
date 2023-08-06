from pysam.libcalignmentfile import AlignmentFile
from typing import Dict

from .containers import MappedRead, Mapping

from collections import namedtuple, defaultdict
from itertools import product
PAlignment = namedtuple('PAlignment',
                        'iss ref revr1 revr2 score startr1 endr1 startr2 endr2 orientation')


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

    for query, alignments in read_mappings.items():
        if alignments.is_single_end():
            cnt_singleend += 1
            continue
        else:
            cnt_pairedend += 1
            score_filtered = {}

            matches = {}
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

            minscore: int = int(bestscore * cutoff_bestscore)
            scoreinsertsizesortedmatches = []
            for score in sorted(matches.keys(), reverse=True):
                if score < minscore:
                    break
                match = matches[score]
                tmpmatch = sorted(match, key=lambda x: x.iss)
                scoreinsertsizesortedmatches = scoreinsertsizesortedmatches + tmpmatch

            for samematch in samematches:
                print(samematch)
                if samematch.score >= minscore:
                    scoreinsertsizesortedmatches.append(samematch)

            if len(scoreinsertsizesortedmatches) == 0:
                continue

            score_filtered[query] = scoreinsertsizesortedmatches

            return score_filtered