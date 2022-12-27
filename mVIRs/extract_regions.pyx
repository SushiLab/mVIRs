import logging
from collections import defaultdict
from pysam.libcfaidx cimport FastxFile, FastxRecord

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map as umap
from libcpp.map cimport map

from operator import itemgetter


cdef inline umap[string, string] load_fasta(string sequence_file):
    cdef umap[string, string] sequences 
    cdef FastxFile fasta = FastxFile(sequence_file)
    cdef FastxRecord entry

    for entry in fasta:
        sequences[entry.name.encode()] = entry.sequence.encode()

    return sequences


cdef read_clipped_file(str clipped_file):

    cdef string line 
    cdef vector[string] splits
    clipped_reads = defaultdict(list)
    cdef dict soft_clipped_positions = {}

    with open(clipped_file, "rb") as handle:
        for line in handle:
            if line.startswith(b'#'):
                continue
            splits = line.strip().split()
            clipped_reads[(splits[0], splits[1])].append((splits[2], splits[3], int(splits[4]), splits[5]))
            if splits[2] == b'S':

                soft_clipped_positions.setdefault((int(splits[4]), splits[5]), 0)
                soft_clipped_positions[(int(splits[4]), splits[5])] += 1
    
    return clipped_reads, soft_clipped_positions 


cdef tuple read_oprs_file(str oprs_filepath):
    
    cdef int pos1, pos2 
    cdef int max_reasonable_insert_size = 0
    cdef int estimated_insert_size = 0
    cdef string line, scaffold
    cdef list splits
    cdef dict oprs_start_to_stop = {}
    
    with open(oprs_filepath, "rb") as handle:
        for line in handle:
            if line.startswith(b'#'):
                if line.startswith(b'#MAX_REASONABLE_INSERTSIZE'):
                    max_reasonable_insert_size = int(line.split(b'=')[1])
                if line.startswith(b'#ESTIMATED_INSERTSIZE'):
                    estimated_insert_size = int(line.split(b'=')[1])
                continue
            splits = line.strip().split()
            if splits[-1] == b'OPR':
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


cdef denoise_softclips(soft_clipped_positions, int softclip_range):

    cdef int softclip_position, softclip_count
    cdef string softclip_scaffold
    cdef tuple winner 
    cdef dict denoised_soft_clip_pos = {}
    cdef dict candidates 

    for (softclip_position, softclip_scaffold), softclip_count in sorted(soft_clipped_positions.items(), 
                                                                        key=itemgetter(1),
                                                                        reverse=True):
        candidates = {}

        for potential_pos in range(softclip_position-softclip_range, softclip_position+softclip_range):
            if (potential_pos, softclip_scaffold) in denoised_soft_clip_pos:
                candidates[(potential_pos, softclip_scaffold)] = denoised_soft_clip_pos[(potential_pos, softclip_scaffold)]

        if len(candidates) == 0:
            denoised_soft_clip_pos[(softclip_position, softclip_scaffold)] = softclip_count
        
        else:
            if len(candidates) == 1:
                winner = list(candidates.keys())[0]
            else:
                # most common candidate
                winner = sorted(candidates.items(), key=itemgetter(1), reverse=True)[0]
                    
            denoised_soft_clip_pos.setdefault(winner, 0) 
            denoised_soft_clip_pos[winner] += softclip_count
        
        return denoised_soft_clip_pos


cdef pair_hard_soft(clipped_reads):
    '''
    Pairing hard/soft clips
    '''
    cdef int softclipped_pos, hardclip_pos
    cdef list hardclipped
    cdef tuple softclipped
    cdef dict soft_to_hardclip_pairs = {}

    for alignments in clipped_reads.values():
        softclipped = [aln for aln in alignments if aln[0] == b'S'][0]
        hardclipped = [aln for aln in alignments if aln[0] == b'H']
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

cdef int get_sum(tuple tpl):
    tuple_sum = sum(tpl[1])
    return tuple_sum

cpdef extract_regions(
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

    cdef umap[string, string] reference_header_2_sequence = load_fasta(reference_fasta_file.encode())
    clipped_reads, soft_clipped_positions = read_clipped_file(clipped_file)
    cdef int soft_cnt = len(soft_clipped_positions)
    cdef int reads_sum = sum(soft_clipped_positions.values())

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
    cdef int soft_clip_percentage, reads, reads_percentage 
    cdef dict filtered_soft_clipped_positions = {}
    
    for (softclip_position, softclip_scaffold), softclip_count in sorted(updated_soft_clipped_positions.items(), 
                                                                         key=itemgetter(1),
                                                                         reverse=True):
        if softclip_count > 1:
            filtered_soft_clipped_positions[(softclip_position, softclip_scaffold)] = softclip_count
    
    soft_clip_percentage = int(len(filtered_soft_clipped_positions) * 100.0 / cnt_tmp)
    reads = sum(filtered_soft_clipped_positions.values())
    reads_percentage = int(reads * 100.0 / sum_tmp)

    logging.info(f'Removing singleton positions finished. {len(soft_clipped_positions)} '
                 f'({soft_clip_percentage}%) positions from {reads} ' 
                 f'({reads_percentage}%) reads were kept.')

    soft_to_hardclip_pairs = pair_hard_soft(clipped_reads)

    logging.info(f'Found {len(soft_to_hardclip_pairs)} hardclip-softclip split ' 
                 f'alignment pairs from {sum(soft_to_hardclip_pairs.values())} ' 
                  'reads determining start/end positons.')

    oprs_start_to_stop, max_reasonable_insert_size, estimated_insert_size = read_oprs_file(opr_file)

    logging.info(f'Addng OPR information from {len(oprs_start_to_stop)} positions '
                 f'and {sum(oprs_start_to_stop.values())} inserts.')

    cdef dict opr_supported_start_ends = {}
    cdef int hsp1, hsp2
    cdef string hsscaffold
    for (hsp1, hsp2, hsscaffold), hscnt in soft_to_hardclip_pairs.items():
        opr_supported_start_ends[(hsp1, hsp2, hsscaffold)] = (hscnt, 0)

    cdef dict distances, filtered_distances, double_filtered_distances
    cdef int deltap1, deltap2, delta
    cdef int oprp1, oprp2
    cdef string oprscaffold
    cdef int min_dev_from_insert_size
    cdef int oprcnt, length, scaffold_length
    cdef int minsize = minmvirlength
    cdef int maxsize = maxmvirlength
    cdef (int, int) entry
    cdef tuple winner

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
            
        if minsize >= length or length >= maxsize:
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
    #minmvirlength=1000, maxmvirlength=1000000, allow_complete_scaffolds=True

    cdef int minoprcount = 1
    cdef int minhscount = 1
    cdef int mincombcount = 5
    cdef int opr_cnt, hs_cnt, start, end
    cdef string scaffold
    cdef dict updated_opr_supported_start_ends = {}
    

    logging.info(f'Filtering: \n\tby length {minsize} <= length <= {maxsize}. '
                 f'\n\t#OPRs >= {minoprcount}\n\t#HARD-SOFT >= {minhscount} '
                 f'\n\t#OPRS + #HARD-SOFT >= {mincombcount} '
                 f'\n\tAllow full scaffolds: {allow_complete_scaffolds}')
    
    for (start, end, scaffold), (hs_cnt, opr_cnt) in sorted(opr_supported_start_ends.items(), key=get_sum, reverse=True):

        # 2:
        if opr_cnt < minoprcount or hs_cnt < minhscount or (opr_cnt + hs_cnt) < mincombcount:
            continue
        updated_opr_supported_start_ends[(start, end, scaffold)] = (hs_cnt, opr_cnt)
        if (opr_cnt + hs_cnt) < mincombcount:
            break

        updated_opr_supported_start_ends[(start, end, scaffold)] = (hs_cnt, opr_cnt)
    opr_supported_start_ends = updated_opr_supported_start_ends
    logging.info(f'Finished filtering. Working with {len(opr_supported_start_ends)} start/end combinations.')

    cdef dict ref_opr_supported_start_ends = {}
    cdef int refstart, refend, ref_hs_cnt, ref_opr_cnt
    cdef string refscaffold
    cdef bint found

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

    cdef string scaffold_sequence
    cdef str subsequence
    cdef double scaffold_coverage

    with open(output_fasta_file, 'wb') as outhandle:
        for (start, end, scaffold), (hs_cnt, opr_cnt) in ref_opr_supported_start_ends.items():
            scaffold_sequence = reference_header_2_sequence[scaffold]
            scaffold_length: int = len(scaffold_sequence)
            subsequence: str = scaffold_sequence[start:end+1].decode()
            scaffold_coverage = round(100.0 * len(subsequence) / scaffold_length, 6)
            outhandle.write(str.encode(f'>{scaffold.decode()}:{start}-{end}\t'
                                       f'OPRs={opr_cnt}-HSs={hs_cnt}-SF={scaffold_coverage}'
                                       f'\n{subsequence}\n'))