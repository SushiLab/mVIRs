import collections
import logging
import sys
from typing import Dict, Generator, List, Tuple

from mVIRs import ( 
    get_read_orientation, 
    insertize_bamfile_by_name,
    _estimate_insert_size,
    _calc_primary_paired_alignment
)


debug = False

PYSAM_BAM_CSOFT_CLIP = 4
PYSAM_BAM_CHARD_CLIP = 5

PAlignment = collections.namedtuple('PAlignment',
                                    'iss, ref revr1 revr2 score startr1 endr1 startr2 endr2 orientation')
SAMLine = collections.namedtuple('SAMLine', 'rev ref rstart rend score cigartuples blocks')


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