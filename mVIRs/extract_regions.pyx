import logging
from collections import defaultdict, Counter
import gzip


cdef dict load_fasta(str sequence_file):
    '''
    Read a fasta file and put it into a dictionary
    :param sequence_file:
    :return:
    '''
    cdef dict sequences = {}  
    cdef dict sequences2 = {} 
    cdef str current_header
    cdef str line, header, tmp_seq 
    cdef list seq, sequence

    if sequence_file.endswith('.gz'):
        handle = gzip.open(sequence_file, 'rt')
    else:
        handle = open(sequence_file)

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
    
    for header, sequence in sequences.items():
        tmp_seq = ''.join(sequence)
        sequences2[header] = tmp_seq
    
    return sequences2

cdef read_clipped_file(str clipped_file):

    cdef str line 
    cdef list splits
    clipped_reads = defaultdict(list)
    soft_clipped_positions = Counter()

    with open(clipped_file) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            splits = line.strip().split()
            clipped_reads[(splits[0], splits[1])].append((splits[2], splits[3], int(splits[4]), splits[5]))
            if splits[2] == 'S':
                soft_clipped_positions[(int(splits[4]), splits[5])] += 1
    
    return clipped_reads, soft_clipped_positions 

cdef read_oprs_file(str oprs_filepath):
    
    cdef int pos1, pos2 
    cdef int max_reasonable_insert_size = 0
    cdef int estimated_insert_size = 0
    cdef str line, scaffold
    cdef list splits
    oprs_start_to_stop = Counter() # (start, stop, scaffold) --> Count
    
    with open(oprs_filepath) as handle:
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

    return oprs_start_to_stop, max_reasonable_insert_size, estimated_insert_size

cdef denoise_softclips(soft_clipped_positions, 
                       int softclip_range):
    cdef int softclip_position, softclip_count
    cdef str softclip_scaffold
    cdef tuple winner 
    denoised_soft_clip_pos = Counter()

    for (softclip_position, softclip_scaffold), softclip_count in soft_clipped_positions.most_common():
        candidates = Counter()

        for potential_pos in range(softclip_position-softclip_range, softclip_position+softclip_range):
            if (potential_pos, softclip_scaffold) in denoised_soft_clip_pos:
                candidates[(potential_pos, softclip_scaffold)] = denoised_soft_clip_pos[(potential_pos, softclip_scaffold)]

        if len(candidates) == 0:
            denoised_soft_clip_pos[(softclip_position, softclip_scaffold)] = softclip_count
        
        elif len(candidates) == 1:
            winner = list(candidates.keys())[0]
            denoised_soft_clip_pos[winner] += softclip_count

        else:
            winner = candidates.most_common()[0]
            denoised_soft_clip_pos[winner] += softclip_count
        
        return denoised_soft_clip_pos


cdef pair_hard_soft(clipped_reads):
    '''
    Pairing hard/soft clips
    '''
    cdef int softclipped_pos, hardclip_pos
    cdef list hardclipped
    cdef tuple softclipped
    soft_to_hardclip_pairs = Counter()

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
                soft_to_hardclip_pairs[(softclipped_pos, hardclip_pos, scaffold)] += 1
            else:
                soft_to_hardclip_pairs[(hardclip_pos, softclipped_pos,  scaffold)] += 1

    return soft_to_hardclip_pairs


def extract_regions(
    str clipped_file, 
    str opr_file, 
    str reference_fasta_file, 
    str output_fasta_file, 
    int minmvirlength=1000, 
    int maxmvirlength=1000000, 
    bint allow_complete_scaffolds=1):
    
    # (start, stop, scaffold) --> Count
    cdef int softclip_range = 20
    # OPRS have an unprecise location of start and end position of virus but start and end are connected
    # Soft-Hard pairs are partly unprecise (But more precise then OPRs) but connect start and end
    # Soft-clips are precise and plenty but they don't connect start with end

    # 1. Create a start-end map with abundance from soft-hard pairs
    # 2. Add OPRs to the start-end map whereever they potentially fit +- 500 bp
    # 3. Get a more precise location with softclips

    logging.info('Finding potential viruses in the genome')
    logging.info('Reading reference fasta')

    cdef dict reference_header_2_sequence = load_fasta(reference_fasta_file)
    clipped_reads, soft_clipped_positions = read_clipped_file(clipped_file)
    cdef int soft_cnt = len(soft_clipped_positions)
    cdef int reads_sum = len(soft_clipped_positions.values())

    logging.info('Start finding start/end positions of non-continous alignment '
                 'regions using clipped and OPR alignments')
    logging.info('Creating initial start/end positions using hard-soft '
                 'alignment pairs')
    logging.info(f'Start denoising {soft_cnt} soft clipped '
                 f'positions from {reads_sum} reads.')

    updated_soft_clipped_positions = denoise_softclips(soft_clipped_positions,
                                                       softclip_range)

    cdef int cnt_tmp = len(soft_clipped_positions)
    cdef int sum_tmp = sum(soft_clipped_positions.values())
    cdef int upd_cnt = len(updated_soft_clipped_positions)
    cdef int upd_sum = sum(updated_soft_clipped_positions.values())

    logging.info(f'Denoising soft clipped reads finished. '
                 f'{upd_cnt} ({int(upd_cnt * 100.0 / cnt_tmp)}%) '
                 f'positions from {upd_sum} ({int(upd_sum * 100.0 / sum_tmp)}%) '
                 'reads were kept.')
    
    cdef int softclip_count 
    cdef int soft_clip_percentage = 0 
    cdef int reads = 0 
    cdef int reads_percentage = 0
    soft_clipped_positions = Counter()
    
    for (softclip_position, softclip_scaffold), softclip_count in updated_soft_clipped_positions.most_common():
        if softclip_count > 1:
            soft_clipped_positions[(softclip_position, softclip_scaffold)] = softclip_count
    
    soft_clip_percentage = int(len(soft_clipped_positions) * 100.0 / cnt_tmp)
    reads = sum(soft_clipped_positions.values())
    reads_percentage = int(reads * 100.0 / sum_tmp)

    logging.info(f'Removing singleton positions finished. {len(soft_clipped_positions)}'
                 f'({soft_clip_percentage}%) positions from {reads}' 
                 f'({reads_percentage}%) reads were kept.')

    soft_to_hardclip_pairs = pair_hard_soft(clipped_reads)

    logging.info(f'Found {len(soft_to_hardclip_pairs)} hardclip-softclip split' 
                 f'alignment pairs from {sum(soft_to_hardclip_pairs.values())}' 
                  'reads determining start/end positons.')

    oprs_start_to_stop, max_reasonable_insert_size, estimated_insert_size = read_oprs_file(opr_file)

    logging.info(f'Adding OPR information from {len(oprs_start_to_stop)} positions and {sum(oprs_start_to_stop.values())} inserts.')

    cdef dict opr_supported_start_ends = {}
    cdef int hsp1, hsp2
    cdef str hsscaffold
    for (hsp1, hsp2, hsscaffold), hscnt in soft_to_hardclip_pairs.items():
        opr_supported_start_ends[(hsp1, hsp2, hsscaffold)] = (hscnt, 0)

    cdef dict distances, filtered_distances, double_filtered_distances
    cdef int deltap1, deltap2, delta
    cdef int oprp1, oprp2
    cdef str oprscaffold
    cdef int min_dev_from_insert_size
    cdef int oprcnt

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
            max_hscnt = max([distance[3] for distance in distances])
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

    logging.info(f'Added OPR information. Now working with '
                 f'{len(opr_supported_start_ends)} start/end combinations.')
    #minmvirlength=1000, maxmvirlength=1000000, allow_complete_scaffolds=True
    cdef int minsize = minmvirlength
    cdef int maxsize = maxmvirlength
    cdef int minoprcount = 1
    cdef int minhscount = 1
    cdef int mincombcount = 5
    cdef int length, scaffold_length, opr_cnt, hs_cnt, start, end
    cdef dict updated_opr_supported_start_ends = {}

    logging.info(f'Filtering: \n\tby length {minsize} <= length <= {maxsize}. '
                 f'\n\t#OPRs >= {minoprcount}\n\t#HARD-SOFT >= {minhscount} '
                 f'\n\t#OPRS + #HARD-SOFT >= {mincombcount} '
                 f'\n\tAllow full scaffolds: {allow_complete_scaffolds}')
    
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

    cdef dict ref_opr_supported_start_ends = {}
    cdef bint found

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

    cdef str scaffold_sequence, subsequence, scaffold_coverage_str
    cdef double scaffold_coverage

    with open(output_fasta_file, 'wb') as outhandle:
        for (start, end, scaffold), (hs_cnt, opr_cnt) in sorted(ref_opr_supported_start_ends.items(), key=lambda i: sum(i[1]),reverse=True):
            scaffold_sequence = reference_header_2_sequence[scaffold]
            scaffold_length: int = len(scaffold_sequence)
            subsequence: str = scaffold_sequence[start:end+1]
            scaffold_coverage = 100.0 * float(len(subsequence))/float(scaffold_length)
            scaffold_coverage_str = "{:0.6f}".format(scaffold_coverage)
            outhandle.write(str.encode(f'>{scaffold}:{start}-{end}\t'
                                       f'OPRs={opr_cnt}-HSs={hs_cnt}-SF={scaffold_coverage_str}'
                                       f'\n{subsequence}\n'))