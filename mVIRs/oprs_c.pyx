cimport cython
import gzip 
import collections
import logging


from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
import pathlib

cpdef get_read_orientation(bint rev):
    cdef str orientation
    if rev:
        orientation = "reverse"
    else:
        orientation = "forward" 
    
    return orientation

cdef dict load_fasta(str sequence_file):
    '''
    Read a fasta file and put it into a dictionary
    :param sequence_file:
    :return:
    '''
    cdef dict sequences = {}  
    cdef dict sequences2 = {} 
    cdef str current_header
    cdef str line, header 
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


cpdef str _calc_orientation(bint revr1, bint revr2, int posr1, int posr2):
    """
    Based on orientation of both reads and their alignment start positions, estimate if a read is PE/SAME/OPR
    :param revr1:
    :param revr2:
    :param posr1:
    :param posr2:
    :return: PAIREDEND, SAME, OPR
    """
    cdef str orientation = 'PAIREDEND'
    if revr1 == revr2:
        orientation = 'SAME'
        return orientation
    cdef int fwpos = 0
    cdef int revpos = 0
    if revr1:
        revpos = posr1
        fwpos = posr2
    else:
        revpos = posr2
        fwpos = posr1
    if revpos < fwpos:
        orientation = 'OPR'
    return orientation

SAMLine = collections.namedtuple('SAMLine', 
                                 'rev ref rstart rend score cigartuples blocks')

cdef inline extract_alignment_info(AlignedSegment alignment, bint need_extended):
    
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
        blocks: list = alignment.get_blocks()
        cigartuples: list = alignment.cigartuples
    samline = SAMLine(rev=reverse, ref=refname, rstart=refstart, 
                        rend=refend, score=ascore, cigartuples=cigartuples, blocks=blocks)
    
    return readname, orientation, samline


def insertize_bamfile_by_name(str bam_file,
                              double max_sam_lines = -1, 
                              double min_coverage = 0.0, 
                              int min_alength = 0, 
                              bint need_extended = 1):
    
    cdef AlignmentFile alignments = AlignmentFile(str(bam_file), "rb")
    
    cdef str current_name = None
    current_insert = collections.defaultdict(list)

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
            current_insert = collections.defaultdict(list)
            current_insert[orientation].append(samline)

    if len(current_insert) != 0:
        yield current_name, current_insert


def extract_regions(
    str clipped_file, 
    str opr_file, 
    str reference_fasta_file, 
    str output_fasta_file, 
    int minmvirlength=1000, 
    int maxmvirlength=1000000, 
    bint allow_complete_scaffolds=1):

    clipped_reads = collections.defaultdict(list)
    soft_clipped_positions = collections.Counter()
    soft_to_hardclip_pairs = collections.Counter() # (start, stop, scaffold) --> Count
    oprs_start_to_stop = collections.Counter() # (start, stop, scaffold) --> Count
    
    cdef int max_reasonable_insert_size = 0
    cdef int estimated_insert_size = 0

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

    with open(clipped_file) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            splits = line.strip().split()
            clipped_reads[(splits[0], splits[1])].append((splits[2], splits[3], int(splits[4]), splits[5]))
            if splits[2] == 'S':
                soft_clipped_positions[(int(splits[4]), splits[5])] += 1

    logging.info(f'Start finding start/end positions of non-continous alignment regions using clipped and OPR alignments')
    logging.info(f'Creating initial start/end positions using hard-soft alignment pairs')
    logging.info(f'Start denoising {len(soft_clipped_positions)} soft clipped positions from {sum(soft_clipped_positions.values())} reads.')
    updated_soft_clipped_positions = collections.Counter()

    for (softclip_position, softclip_scaffold), softclip_count in soft_clipped_positions.most_common():

        candidates = collections.Counter()
        for potential_pos in range(softclip_position-softclip_range, softclip_position+softclip_range):
            if (potential_pos, softclip_scaffold) in updated_soft_clipped_positions:
                candidates[(potential_pos, softclip_scaffold)] = updated_soft_clipped_positions[(potential_pos, softclip_scaffold)]

        if len(candidates) == 0:
            updated_soft_clipped_positions[(softclip_position, softclip_scaffold)] = softclip_count
        elif len(candidates) == 1:
            winner = list(candidates.keys())[0]
            updated_soft_clipped_positions[winner] += softclip_count
        else:
            winner = candidates.most_common()[0]
            updated_soft_clipped_positions[winner] += softclip_count

    cdef int cnt_tmp = len(soft_clipped_positions)
    cdef int sum_tmp = sum(soft_clipped_positions.values())
    logging.info(
        f'Denoising soft clipped reads finished. {len(updated_soft_clipped_positions)} ({int(len(updated_soft_clipped_positions) * 100.0 / cnt_tmp)}%) positions from {sum(updated_soft_clipped_positions.values())} ({int(sum(updated_soft_clipped_positions.values()) * 100.0 / sum_tmp)}%) reads were kept.')
    soft_clipped_positions = collections.Counter()
    for (softclip_position, softclip_scaffold), softclip_count in updated_soft_clipped_positions.most_common():
        if softclip_count > 1:
            soft_clipped_positions[(softclip_position, softclip_scaffold)] = softclip_count
            soft_clip_percentage = int(len(soft_clipped_positions)*100.0 / cnt_tmp)
            reads = sum(soft_clipped_positions.values())
            reads_percentage = int(reads*100.0 / sum_tmp)
    logging.info(f'Removing singleton positions finished. {len(soft_clipped_positions)}'
                 f'({soft_clip_percentage}%) positions from {reads}' 
                 f'({reads_percentage}%) reads were kept.')

    '''
    Pairing hard/soft clips
    '''
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
    logging.info(f'Found {len(soft_to_hardclip_pairs)} hardclip-softclip split' 
                 f'alignment pairs from {sum(soft_to_hardclip_pairs.values())}' 
                  'reads determining start/end positons.')


    with open(opr_file) as handle:
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

    logging.info(f'Adding OPR information from {len(oprs_start_to_stop)} positions and {sum(oprs_start_to_stop.values())} inserts.')

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
    minsize = minmvirlength
    maxsize = maxmvirlength
    minoprcount = 1
    minhscount = 1
    mincombcount = 5
    
    logging.info(f'Filtering: \n\tby length {minsize} <= length <= {maxsize}. \n\t#OPRs >= {minoprcount}\n\t#HARD-SOFT >= {minhscount}\n\t#OPRS + #HARD-SOFT >= {mincombcount}\n\tAllow full scaffolds: {allow_complete_scaffolds}')
    updated_opr_supported_start_ends = {}
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


    ref_opr_supported_start_ends = {}
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

    with open(output_fasta_file, 'w') as outhandle:
        for (start, end, scaffold), (hs_cnt, opr_cnt) in sorted(ref_opr_supported_start_ends.items(), key=lambda i: sum(i[1]),reverse=True):
            scaffold_sequence = reference_header_2_sequence[scaffold]
            scaffold_length = len(scaffold_sequence)
            subsequence = scaffold_sequence[start:end+1]
            scaffold_coverage = 100.0 * float(len(subsequence))/float(scaffold_length)
            scaffold_coverage = "{:0.6f}".format(scaffold_coverage)
            outhandle.write(f'>{scaffold}:{start}-{end}\tOPRs={opr_cnt}-HSs={hs_cnt}-SF={scaffold_coverage}\n{subsequence}\n')