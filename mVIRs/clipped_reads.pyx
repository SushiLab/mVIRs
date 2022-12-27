import logging
from collections import defaultdict
import sys 

from .oprs_c import insertize_bamfile_by_name


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