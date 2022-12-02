cimport cython
import gzip 
from collections import Counter, defaultdict, namedtuple
import logging
import statistics
from typing import List, Dict, Tuple

from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
import pathlib

PAlignment = namedtuple('PAlignment',
                        'iss, ref revr1 revr2 score startr1 endr1 startr2 endr2 orientation')
SAMLine = namedtuple('SAMLine', 'rev ref rstart rend score cigartuples blocks')

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


cdef inline str _calc_orientation(bint revr1, bint revr2, int posr1, int posr2):
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

SAMLine = namedtuple('SAMLine', 
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
    current_insert: defaultdict = defaultdict(list)

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


def extract_regions(
    str clipped_file, 
    str opr_file, 
    str reference_fasta_file, 
    str output_fasta_file, 
    int minmvirlength=1000, 
    int maxmvirlength=1000000, 
    bint allow_complete_scaffolds=1):

    clipped_reads: defaultdict = defaultdict(list)
    soft_clipped_positions: Counter = Counter()
    soft_to_hardclip_pairs: Counter = Counter() # (start, stop, scaffold) --> Count
    oprs_start_to_stop: Counter = Counter() # (start, stop, scaffold) --> Count
    
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
    updated_soft_clipped_positions: Counter = Counter()

    for (softclip_position, softclip_scaffold), softclip_count in soft_clipped_positions.most_common():

        candidates: Counter = Counter()

        for potential_pos in range(softclip_position-softclip_range, softclip_position+softclip_range):
            if (potential_pos, softclip_scaffold) in updated_soft_clipped_positions:
                candidates[(potential_pos, softclip_scaffold)] = updated_soft_clipped_positions[(potential_pos, softclip_scaffold)]

        if len(candidates) == 0:
            updated_soft_clipped_positions[(softclip_position, softclip_scaffold)] = softclip_count
        elif len(candidates) == 1:
            winner: tuple = list(candidates.keys())[0]
            updated_soft_clipped_positions[winner] += softclip_count
        else:
            winner = candidates.most_common()[0]
            updated_soft_clipped_positions[winner] += softclip_count

    cdef int cnt_tmp = len(soft_clipped_positions)
    cdef int sum_tmp = sum(soft_clipped_positions.values())
    logging.info(
        f'Denoising soft clipped reads finished. {len(updated_soft_clipped_positions)} ({int(len(updated_soft_clipped_positions) * 100.0 / cnt_tmp)}%) positions from {sum(updated_soft_clipped_positions.values())} ({int(sum(updated_soft_clipped_positions.values()) * 100.0 / sum_tmp)}%) reads were kept.')
    soft_clipped_positions: Counter = Counter()
    for (softclip_position, softclip_scaffold), softclip_count in updated_soft_clipped_positions.most_common():
        if softclip_count > 1:
            soft_clipped_positions[(softclip_position, softclip_scaffold)] = softclip_count
            soft_clip_percentage: int = int(len(soft_clipped_positions)*100.0 / cnt_tmp)
            reads: int = sum(soft_clipped_positions.values())
            reads_percentage: int = int(reads*100.0 / sum_tmp)
    logging.info(f'Removing singleton positions finished. {len(soft_clipped_positions)}'
                 f'({soft_clip_percentage}%) positions from {reads}' 
                 f'({reads_percentage}%) reads were kept.')

    '''
    Pairing hard/soft clips
    '''
    for alignments in clipped_reads.values():
        softclipped: tuple = [aln for aln in alignments if aln[0] == 'S'][0]
        hardclipped: list = [aln for aln in alignments if aln[0] == 'H']
        if len(hardclipped) == 0:
            continue
        scaffold: str = softclipped[3]
        softclipped_pos: int = softclipped[2]
        for hardclip in hardclipped:
            hardclip_pos: int = hardclip[2]
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
            splits: list = line.strip().split()
            if splits[-1] == 'OPR':
                pos1: int = int(splits[6])
                pos2: int = int(splits[7])
                scaffold = splits[1]
                if pos1 < pos2:
                    oprs_start_to_stop[(pos1, pos2, scaffold)] += 1
                else:
                    oprs_start_to_stop[(pos2, pos1, scaffold)] += 1

    logging.info(f'Adding OPR information from {len(oprs_start_to_stop)} positions and {sum(oprs_start_to_stop.values())} inserts.')

    opr_supported_start_ends: dict = {}
    for (hsp1, hsp2, hsscaffold), hscnt in soft_to_hardclip_pairs.items():
        opr_supported_start_ends[(hsp1, hsp2, hsscaffold)] = (hscnt, 0)

    for (oprp1, oprp2, oprscaffold), oprcnt in oprs_start_to_stop.items():
        distances: dict = {}
        for (hsp1, hsp2, hsscaffold), hscnt in soft_to_hardclip_pairs.items():
            if oprscaffold == hsscaffold:
                deltap1: int = abs(hsp1 - oprp1)
                deltap2: int = abs(hsp2 - oprp2)
                delta: int = deltap1 + deltap2
                if delta <= max_reasonable_insert_size:
                    distances[(hsp1, hsp2, hsscaffold, hscnt)] = delta
        winner: tuple = (oprp1, oprp2, oprscaffold, 0)
        if len(distances) > 0:
            max_hscnt: int = max([distance[3] for distance in distances])
            filtered_distances: dict = {distance:insert_size for (distance,insert_size) in distances.items() if distance[3] == max_hscnt}
            if len(filtered_distances) == 1:
                winner = list(filtered_distances.items())[0][0]
            else:
                min_dev_from_insert_size: int = min([abs(isize - estimated_insert_size) for isize in filtered_distances.values()])
                double_filtered_distances: dict = {distance: isize for (distance, isize) in filtered_distances.items() if abs(isize - estimated_insert_size) == min_dev_from_insert_size}
                winner = list(double_filtered_distances.items())[0][0]
        else:
            winner: tuple = (oprp1, oprp2, oprscaffold, 0)

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
    cdef int length, scaffold_length
    
    logging.info(f'Filtering: \n\tby length {minsize} <= length <= {maxsize}. \n\t#OPRs >= {minoprcount}\n\t#HARD-SOFT >= {minhscount}\n\t#OPRS + #HARD-SOFT >= {mincombcount}\n\tAllow full scaffolds: {allow_complete_scaffolds}')
    cdef dict updated_opr_supported_start_ends = {}
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

    cdef str scaffold_sequence, subsequence

    with open(output_fasta_file, 'w') as outhandle:
        for (start, end, scaffold), (hs_cnt, opr_cnt) in sorted(ref_opr_supported_start_ends.items(), key=lambda i: sum(i[1]),reverse=True):
            scaffold_sequence = reference_header_2_sequence[scaffold]
            scaffold_length: int = len(scaffold_sequence)
            subsequence: str = scaffold_sequence[start:end+1]
            scaffold_coverage = 100.0 * float(len(subsequence))/float(scaffold_length)
            scaffold_coverage = "{:0.6f}".format(scaffold_coverage)
            outhandle.write(f'>{scaffold}:{start}-{end}\tOPRs={opr_cnt}-HSs={hs_cnt}-SF={scaffold_coverage}\n{subsequence}\n')


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

    for query, alignments in insert2alignments:
        if len(alignments) < 2:
            cnt_singleend += 1
            continue
        else:
            cnt_pairedend += 1
            matches: defaultdict = defaultdict(list)
            samematches: list = []
            bestscore: int = 0
            alnr1: List[SAMLine]
            alnr2: List[SAMLine]

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
                                       score=alnr1.score + alnr2.score, startr1=alnr1.rstart, endr1=alnr1.rend, startr2=alnr2.rstart, endr2=alnr2.rend,
                                       orientation=orientation, ref=alnr1.ref))
                        continue
                    score: int = alnr1.score + alnr2.score
                    if bestscore < score:
                        bestscore = score
                    matches[alnr1.score + alnr2.score].append(
                        PAlignment(iss=iss, revr1=alnr1.rev, revr2=alnr2.rev, score=score,
                                   startr1=alnr1.rstart, endr1=alnr1.rend, startr2=alnr2.rstart,
                                   endr2=alnr2.rend, orientation=orientation, ref=alnr1.ref))


            minscore: int = int(bestscore * cutoff_bestscore)
            scoreinsertsizesortedmatches: DefaultDict[str, List[PAlignment]] = []
            score: int
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


def _estimate_insert_size(insert2alignments: Dict[str, Dict[str, List[SAMLine]]]):
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
    
    while True:
        iteration += 1
        logging.info(
            'Estimating insert size by mean/stdev convergence. Iteration:\t{}. Min insert size:\t{}. Max insert size:\t{}. Mean insert size:\t{}. Stdev insert size:\t{}'.format(
                iteration, format(min(so), ',d'), format(max(so), ',d'), format(int(statistics.mean(so)), ',d'), format(int(statistics.stdev(so)), ',d')))
        tmp: List[int] = []
        mean: int = int(statistics.mean(so))
        stdev: int = int(statistics.pstdev(so))
        #print(stdev)
        minval: int = mean - 7 * stdev
        if minval < 0:
            minval = 0
        maxval: int = mean + 7 * stdev
        val: int
        for val in so:
            if minval > val:
                continue
            if maxval < val:
                continue
            tmp.append(val)
        so: List[int] = tmp

        if minvaltmp == minval and maxval == maxvaltmp:
            break

        minvaltmp: int = minval
        maxvaltmp: int = maxval


    # For mapping Illumina short-insert reads to the human genome, x is about 6-7 sigma away from the mean
    # http://bio-bwa.sourceforge.net/bwa.shtml section Estimating Insert Size Distribution

    minvaltmp: int = mean - 7 * int(statistics.pstdev(so))

    if minvaltmp != 0:
        minvaltmp = 0

    maxvaltmp: int = mean + 7 * int(statistics.pstdev(so))
    logging.info(
        'Estimating insert size by mean/stdev convergence. Iteration:\t{}. Min insert size:\t{}. Max insert size:\t{}. Mean insert size:\t{}. Stdev insert size:\t{}'.format(
            iteration, format(minvaltmp, ',d'), format(maxvaltmp, ',d'), format(int(statistics.mean(so)), ',d'),
            format(int(statistics.stdev(so)), ',d')))

    return minvaltmp, maxvaltmp, int(statistics.median(so))


def _calc_primary_paired_alignment(insert2alignments: Dict[str, Dict[str, List[SAMLine]]], 
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

    for (qname, alignments) in _generate_paired_alignments(insert2alignments):
        if len(alignments) == 1:  # This is a unique alignment. No filtering required.
            yield (qname, alignments, True)
        else:
            # Check if there are any alignments with PE and reasonable insert size
            pe_candidates: List[PAlignment] = []
            alignment: PAlignment
            for alignment in alignments:
                if alignment.orientation == 'PAIREDEND' and alignment.iss >= minreasonable_insertsize and alignment.iss <= maxreasonable_insertsize:
                    pe_candidates.append(alignment)
            if len(pe_candidates) != 0:
                # if there are alignments with PE and reasonable insertsize then pick the one(s) with the best score
                bestscore: int = -1
                pe_bestalns: List[PAlignment] = []
                alignment: PAlignment
                for alignment in sorted(pe_candidates, key=lambda x: x.score, reverse=True):
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
