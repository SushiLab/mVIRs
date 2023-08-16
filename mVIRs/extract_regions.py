import logging

from mVIRs.utils import load_fasta
from operator import itemgetter


def read_clipped(clipped_file: str):
    """
    Reads clipped file.

    Args:
        clipped_file (str): Path to clipped file.
    """
    clipped_reads = {}
    soft_clipped_positions = {}

    with open(clipped_file, "r") as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            splits = line.strip().split()
            splits = line.strip().split()
            insert = splits[0]
            orientation = splits[1]
            clip_type = splits[2]
            direction = splits[3]
            position = int(splits[4])
            scaffold = splits[5]
            clipped_reads.setdefault((insert, orientation), [])
            clipped_reads[(insert, orientation)].append((clip_type, direction, position, scaffold))

            if clip_type == 'S':
                soft_clipped_positions.setdefault((position, scaffold), 0)
                soft_clipped_positions[(position, scaffold)] += 1

    return clipped_reads, soft_clipped_positions

def read_oprs_file(oprs_filepath: str):
    """
    Reads OPRs file.

    Args:
        oprs_filepath (str): Path to OPRs file.

    Returns:
        oprs_start_to_stop (dict): Dictionary with OPRs start to stop.
        max_reasonable_insert_size (int): Maximum reasonable insert size.
        estimated_insert_size (int): Estimated insert size.
    """



    oprs_start_to_stop = {}
    with open(oprs_filepath, "r") as handle:
        for line in handle:
            if line.startswith('#'):
                if line.startswith('#MAX_REASONABLE_INSERTSIZE'):
                    max_reasonable_insert_size = int(line.split('=')[1])
                elif line.startswith('#ESTIMATED_INSERTSIZE'):
                    estimated_insert_size = int(line.split('=')[1])
                continue
            splits = line.strip().split()
            if splits[-1] == 'OPR':
                pos1 = int(splits[6])
                pos2 = int(splits[7])
                scaffold = splits[1]

                if pos1 < pos2:
                    # all fine
                    pass
                else:
                    # order is mixed
                    pos1, pos2 = pos2, pos1

                oprs_start_to_stop.setdefault((pos1, pos2, scaffold), 0)
                oprs_start_to_stop[(pos1, pos2, scaffold)] += 1

    return oprs_start_to_stop, max_reasonable_insert_size, estimated_insert_size

def pair_hard_soft(clipped_reads):
    '''
    Pairing hard/soft clips
    '''
    soft_to_hardclip_pairs = {}

    for alignments in clipped_reads.values():
        softclipped = [aln for aln in alignments if aln[0] == 'S'][0]
        hardclipped = [aln for aln in alignments if aln[0] == 'H']
        if len(hardclipped) == 0:
            continue
        scaffold = softclipped[3]
        softclipped_pos = softclipped[2]
        for hardclip in hardclipped:
            hardclip_pos = hardclip[2]
            if hardclip_pos > softclipped_pos:
                # all fine, hardclip is in the beggining
                pass
            else:
                # reverse the positions, hardclip is in the end
                softclipped_pos, hardclip_pos = hardclip_pos, softclipped_pos

            soft_to_hardclip_pairs.setdefault((softclipped_pos, hardclip_pos, scaffold), 0)
            soft_to_hardclip_pairs[(softclipped_pos, hardclip_pos, scaffold)] += 1

    return soft_to_hardclip_pairs


def denoise_softclips(soft_clipped_positions, softclip_range: int):

    denoised_soft_clip_pos = {}

    for (softclip_position, softclip_scaffold), softclip_count in sorted(soft_clipped_positions.items(),
                                                                        key=itemgetter(1),
                                                                        reverse=True):
        candidates = {}

        for potential_pos in range(softclip_position - softclip_range,
                                   softclip_position + softclip_range):
            if (potential_pos, softclip_scaffold) in denoised_soft_clip_pos:
                candidates[(potential_pos, softclip_scaffold)] = denoised_soft_clip_pos[(potential_pos, softclip_scaffold)]

        if len(candidates) == 0:
            denoised_soft_clip_pos[(softclip_position, softclip_scaffold)] = softclip_count

        else:
            if len(candidates) == 1:
                winner = list(candidates.keys())[0]
            else:
                # most common candidate
                winner = max(candidates.items(), key=itemgetter(1))[0]

            denoised_soft_clip_pos.setdefault(winner, 0)
            denoised_soft_clip_pos[winner] += softclip_count

        return denoised_soft_clip_pos


def filter_singletons(soft_clip_positions):
    filtered_soft_clip = {}

    for (softclip_position, softclip_scaffold), softclip_count in sorted(soft_clip_positions.items(),
                                                                         key=itemgetter(1),
                                                                         reverse=True):
        if softclip_count > 1:
            filtered_soft_clip[(softclip_position, softclip_scaffold)] = softclip_count

    return filtered_soft_clip


def extract_regions(
    clipped_file,
    opr_file,
    reference_fasta_file,
    output_fasta_file,
    min_length=1000,
    max_length=1000000,
    allow_complete_scaffolds=1):

    # (start, stop, scaffold) --> Count
    # OPRS have an unprecise location of start and end position of virus but start and end are connected
    # Soft-Hard pairs are partly unprecise (But more precise then OPRs) but connect start and end
    # Soft-clips are precise and plenty but they don't connect start with end

    # 1. Create a start-end map with abundance from soft-hard pairs
    # 2. Add OPRs to the start-end map whereever they potentially fit +- 500 bp
    # 3. Get a more precise location with softclips

    logging.info('Finding potential viruses in the genome')

    reference_header_2_sequence = load_fasta(reference_fasta_file)
    clipped_reads, soft_clipped_positions = read_clipped(clipped_file)
    soft_cnt = len(soft_clipped_positions)
    reads_sum = sum(soft_clipped_positions.values())

    logging.info('Start finding start/end positions of non-continous alignment '
                 'regions using clipped and OPR alignments')
    logging.info('Creating initial start/end positions using hard-soft '
                 'alignment pairs')
    logging.info(f'Start denoising {soft_cnt} soft clipped '
                 f'positions from {reads_sum} reads.')

    updated_soft_clipped_positions = denoise_softclips(soft_clipped_positions,
                                                       softclip_range = 20)

    tmp_cnt = len(soft_clipped_positions)
    tmp_sum = sum(soft_clipped_positions.values())
    upd_cnt = len(updated_soft_clipped_positions)
    upd_sum = sum(updated_soft_clipped_positions.values())

    logging.info(f'Denoising soft clipped reads finished. '
                 f'{upd_cnt} ({int(upd_cnt * 100.0 / tmp_cnt)}%) '
                 f'positions from {tmp_sum} ({int(upd_sum * 100.0 / tmp_sum)}%) '
                 'reads were kept.')

    # Removing singletons

    filtered_soft_clipped_positions = filter_singletons(updated_soft_clipped_positions)

    soft_clip_percentage = int(len(filtered_soft_clipped_positions) * 100.0 / tmp_cnt)
    reads = sum(filtered_soft_clipped_positions.values())
    reads_percentage = int(reads * 100.0 / tmp_sum)

    logging.info(f'Removing singleton positions finished. {len(soft_clipped_positions)} '
                 f'({soft_clip_percentage}%) positions from {reads} '
                 f'({reads_percentage}%) reads were kept.')

    # Pair hard & soft clips

    soft_to_hardclip_pairs = pair_hard_soft(clipped_reads)

    logging.info(f'Found {len(soft_to_hardclip_pairs)} hardclip-softclip split '
                 f'alignment pairs from {sum(soft_to_hardclip_pairs.values())} '
                  'reads determining start/end positons.')

    # read opr file
    oprs_start_to_stop, max_reasonable_insert_size, estimated_insert_size = read_oprs_file(opr_file)

    logging.info(f'Addng OPR information from {len(oprs_start_to_stop)} positions '
                 f'and {sum(oprs_start_to_stop.values())} inserts.')


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

        if len(distances) > 0:
            max_hscnt = max([distance[3] for distance in distances.keys()])
            filtered_distances = {distance:insert_size for (distance,insert_size) in distances.items() if distance[3] == max_hscnt}
            if len(filtered_distances) == 1:
                winner = tuple(filtered_distances.keys())[0]
            else:
                min_dev_from_insert_size = min([abs(isize - estimated_insert_size) for isize in filtered_distances.values()])
                double_filtered_distances = {distance: isize for (distance, isize) in filtered_distances.items() if abs(isize - estimated_insert_size) == min_dev_from_insert_size}
                winner = tuple(double_filtered_distances.keys())[0]
        else:
            winner = (oprp1, oprp2, oprscaffold, 0)

        ## filtering rule 1: length
        length = winner[1] - winner[0] + 1

        if min_length >= length or length >= max_length:
            continue

        ## filtering rule 3: length relative to scaffold
        scaffold_length = len(reference_header_2_sequence[winner[2]])
        if not allow_complete_scaffolds:
            if scaffold_length * 0.99 <  length:
                continue

        entry = opr_supported_start_ends.get((winner[0], winner[1], winner[2]), (0,0))
        entry = (entry[0], entry[1] + oprcnt)
        opr_supported_start_ends[(winner[0], winner[1], winner[2])] = entry

    logging.info(f'Added OPR information. Now working with '
                 f'{len(opr_supported_start_ends)} start/end combinations.')

    minoprcount = 1
    minhscount = 1
    mincombcount = 5
    updated_opr_supported_start_ends = {}

    def get_sum(tpl):
        tuple_sum = sum(tpl[1])
        return tuple_sum


    for (start, end, scaffold), (hs_cnt, opr_cnt) in sorted(opr_supported_start_ends.items(),
                                                            key=get_sum, reverse=True):
        ## filtering rule 2: softclip and hardclip counts
        ## moved here, because total counts are required
        if opr_cnt < minoprcount or hs_cnt < minhscount or (opr_cnt + hs_cnt) < mincombcount:
            continue
        updated_opr_supported_start_ends[(start, end, scaffold)] = (hs_cnt, opr_cnt)
        # if total sum is lower than threshold the loop can be broken
        # as every next entry will satisfy this criteria because of sorting before
        if (opr_cnt + hs_cnt) < mincombcount:
            break

    opr_supported_start_ends = updated_opr_supported_start_ends
    logging.info(f'Finished filtering. Working with {len(opr_supported_start_ends)} start/end combinations.')

    ref_opr_supported_start_ends = {}

    softclip_range = 20

    for (start, end, scaffold), (hs_cnt, opr_cnt) in opr_supported_start_ends.items():
        found = False
        for (refstart, refend, refscaffold), (ref_hs_cnt, ref_opr_cnt) in ref_opr_supported_start_ends.items():
            if refscaffold == scaffold:
                if refstart - softclip_range <= start <= refstart + softclip_range:
                    if refend - softclip_range <= end <= refend + softclip_range:
                        ref_opr_supported_start_ends[(refstart, refend, refscaffold)] = (hs_cnt + ref_hs_cnt, opr_cnt + ref_opr_cnt)
                        found = True
                        break
        if not found:
            ref_opr_supported_start_ends[(start, end, scaffold)] = (hs_cnt, opr_cnt)

    with open(output_fasta_file, 'w') as outhandle:
         for (start, end, scaffold), (hs_cnt, opr_cnt) in ref_opr_supported_start_ends.items():
            scaffold_sequence = reference_header_2_sequence[scaffold]
            scaffold_length = len(scaffold_sequence)
            subsequence = scaffold_sequence[start : end + 1]
            scaffold_coverage = round(100.0 * len(subsequence) / scaffold_length, 6)
            outhandle.write(f'>{scaffold}:{start}-{end}\t'
                            f'OPRs={opr_cnt}-HSs={hs_cnt}-SF={scaffold_coverage}'
                            f'\n{subsequence}\n')