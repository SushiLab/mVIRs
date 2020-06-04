import pysam
import pathlib
import collections
import statistics
from typing import List, Dict, Tuple, Counter, DefaultDict, Iterable, IO, Generator
import logging
import argparse
import sys





PAlignment = collections.namedtuple('PAlignment',
                                    'iss, ref revr1 revr2 score startr1 endr1 startr2 endr2 orientation')
SAMLine = collections.namedtuple('SAMLine', 'rev ref rstart rend score')


def _calc_orientation(revr1: bool, revr2: bool, posr1: int, posr2: int) -> str:
    """
    Based on orientation of both reads and their alignment start positions, estimate if a read is PE/SAME/OPR
    :param revr1:
    :param revr2:
    :param posr1:
    :param posr2:
    :return: PAIREDEND, SAME, OPR
    """
    orientation: str = 'PAIREDEND'
    if revr1 == revr2:
        orientation = 'SAME'
        return orientation
    fwpos: int = 0
    revpos: int = 0
    if revr1:
        revpos = posr1
        fwpos = posr2
    else:
        revpos = posr2
        fwpos = posr1
    if revpos < fwpos:
        orientation = 'OPR'
    return orientation


def _generate_paired_alignments(insert2alignments: Dict[str, Dict[str, List[SAMLine]]],
                                cutoff_bestscore: float = 0.95) -> Generator[
    Tuple[str, DefaultDict[str, List[PAlignment]]], None, None]:
    """
    Iterates over all inserts and generates pairs with a cross product. Pairs with best score are returned.
    Careful --> throws away single end mappers. Considers only paired reads
    :param insert2alignments:
    :return:
    """

    cnt_singleend: int = 0
    cnt_pairedend: int = 0
    query: str
    alignments: Dict[str, List[SAMLine]]
    for query, alignments in insert2alignments.items():
        if len(alignments) < 2:
            cnt_singleend += 1
            continue
        else:
            cnt_pairedend += 1
            matches: DefaultDict[str, List[PAlignment]] = collections.defaultdict(list)
            samematches: List[PAlignment] = []
            bestscore: int = 0
            alnr1: List[SAMLine]
            alnr2: List[SAMLine]

            for alnr1 in alignments['R1']:
                for alnr2 in alignments['R2']:
                    if alnr1.ref != alnr2.ref:  # if both reads align to 2 different references we continue
                        continue
                    orientation: str = _calc_orientation(alnr1.rev, alnr2.rev, alnr1.rstart, alnr2.rstart)
                    isp = abs(alnr1.rstart - alnr2.rstart)
                    positions = [alnr1.rstart, alnr2.rstart, alnr1.rend, alnr2.rend]
                    iss = abs(max(positions) - min(positions))
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
                match: List[PAlignment] = matches[score]
                tmpmatch: Iterable[PAlignment] = sorted(match, key=lambda x: x.iss)
                scoreinsertsizesortedmatches = scoreinsertsizesortedmatches + tmpmatch
            samematch: PAlignment
            for samematch in samematches:
                if samematch.score >= minscore:
                    scoreinsertsizesortedmatches.append(samematch)
            if len(scoreinsertsizesortedmatches) == 0:
                continue
            yield (query, scoreinsertsizesortedmatches)


def _read_bam_file(bam_file: pathlib.Path, max_sam_lines: int = -1) -> Dict[str, Dict[str, List[SAMLine]]]:
    """
    Reads a bam file into memory and splits by R1/R2/S. Each insert is either singleton or paired.
    Singletons are removed
    :param bamfile:
    :return: Dictionary with inserts and the location where reads align
    """
    logging.info('Start reading alignments from file:\t{}'.format(bam_file))
    bf: pysam.libcalignmentfile.AlignmentFile = pysam.AlignmentFile(str(bam_file), "rb")

    insert2alignments: Dict[str, Dict[str, List[SAMLine]]] = collections.defaultdict(
        lambda: collections.defaultdict(list))
    cnt: int
    alignment: pysam.AlignedSegment
    for cnt, alignment in enumerate(bf):
        if cnt % 1000000 == 0:
            logging.info('Alignments read:\t{}\tInserts found:\t{}'.format(format(cnt, ',d'), format(len(insert2alignments), ',d')))
        qname: str = alignment.qname
        splits: List[str] = qname.rsplit('/', 1)
        qname: str = splits[0]

        if splits[1] == 'S':  # Remove singletons
            continue
        ascore: int = alignment.get_tag('AS')

        orientation: str = splits[1]
        reverse: bool = True if alignment.is_reverse else False
        refname: str = alignment.reference_name
        refstart: int = alignment.reference_start
        refend: int = alignment.reference_end


        insert2alignments[qname][orientation].append(
            SAMLine(rev=reverse, ref=refname, rstart=refstart, rend=refend, score=ascore))
        if max_sam_lines <= cnt:
            break
    logging.info('Alignments read:\t{}\tInserts found:\t{}'.format(format(cnt, ',d'), format(len(insert2alignments), ',d')))

    bf.close()
    return insert2alignments


def _estimate_insert_size(insert2alignments: Dict[str, Dict[str, List[SAMLine]]]) -> Tuple[int, int, int]:
    """
    Goes through UNIQUE PAIRED END alignments and estimates the median insert size.
    :param insert2alignments:
    :return:
    """

    cnter: Counter = collections.Counter()
    insertsizes: List[int] = []

    alignments: List[PAlignment]
    for (_, alignments) in _generate_paired_alignments(insert2alignments):
        if len(alignments) == 1 and alignments[0].orientation == 'PAIREDEND':
            insertsizes.append(alignments[0].iss)

    so: List[int] = insertsizes
    minvaltmp: int = -1
    maxvaltmp: int = -1
    iteration = 0
    while True:
        iteration += 1
        logging.info(
            'Estimating insert size by mean/stdev convergence. Iteration:\t{}. Min insert size:\t{}. Max insert size:\t{}'.format(
                iteration, format(min(so), ',d'), format(max(so), ',d')))
        tmp: List[int] = []
        mean: int = int(statistics.mean(so))
        stdev: int = int(statistics.pstdev(so))
        minval: int = mean - 2 * stdev
        if minval < 0:
            minval = 0
        maxval: int = mean + 2 * stdev
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

    minvaltmp: int = minvaltmp - 7 * int(statistics.pstdev(so))
    if minvaltmp != 0:
        minvaltmp = 0
    maxvaltmp: int = maxvaltmp + 7 * int(statistics.pstdev(so))
    return minvaltmp, maxvaltmp, int(statistics.median(so))


def _calc_primary_paired_alignment(insert2alignments: Dict[str, Dict[str, List[SAMLine]]], insertsize: int,
                                   minreasonable_insertsize: int, maxreasonable_insertsize: int) -> Generator[
    Tuple[str, List[PAlignment], bool], None, None]:
    """
    Go through all paired alignments and find one that has a reasonable good score and is close to the insert size.
    1. if there is a PE alignment with reasonable insert size --> take that one (and all other equally good PE with reasonable IS)
    2. Otherwise take the best one
    3. If there are multiple best ones return all and state add the false flag to indicate that alignment could not be resolved to a single best alignment

    :param insert2alignments:
    :param insertsize:
    :return:
    """

    qname: str
    alignments: List[PAlignment]

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













def startup():
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=logging.INFO)
    logging.info('Starting OPR Finder')


def shutdown(status=0):
    logging.info('Finishing OPR Finder with status:\t{}'.format(status))
    sys.exit(status)


def find() -> None:
    '''
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
    '''

    parser = argparse.ArgumentParser(
        description='Takes a BAM alignment file and detects OPRs.')

    parser.add_argument('-i', action='store', dest='bam', help='Input BAM/SAM file',
                        default='', required=True)
    parser.add_argument('-o', action='store', dest='oprfile',
                        help='OPR output file',
                        default=None, required=True)

    try:
        results = parser.parse_args(sys.argv[2:])
    except:
        parser.print_help()
        shutdown(1)



    bam_file: pathlib.Path = pathlib.Path(results.bam)
    opr_file = pathlib.Path = pathlib.Path(results.oprfile)

    logging.info('Input BAM File:\t{}'.format(bam_file))
    logging.info('Output OPR File:\t{}'.format(opr_file))


    # START DEBUG PARAMETERS

    max_sam_lines = 1000000
    # END DEBUG PARAMETERS

    insert2alignments: Dict[str, Dict[str, List[SAMLine]]] = _read_bam_file(bam_file, max_sam_lines)
    logging.info('Start estimating insert size from paired end alignments.')
    minreasonable_insertsize, maxreasonable_insertsize, estimated_insertsize = _estimate_insert_size(insert2alignments)
    logging.info('Finished estimating insert size from paired end alignments.')
    logging.info('Min reasonable insert size:\t{}'.format(minreasonable_insertsize))
    logging.info('Max reasonable insert size:\t{}'.format(maxreasonable_insertsize))





    with open(opr_file, 'w') as handle:
        handle.write(
            '#MIN_REASONABLE_INSERTSIZE={}\n#MAX_REASONABLE_INSERTSIZE={}\n#READNAME\tREFERENCE\tINSERT_SIZE\tORIENTATION_R1\tORIENTATION_R2\tSCORE\tPOSR1\tPOSR2\tQLENGTHR1\tQLENGTHR2\tPAIR_LAYOUT\n'.format(
                minreasonable_insertsize, maxreasonable_insertsize))



        template: str = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'


        logging.info('Start screening for OPRs and Paired-End inserts with unreasonable insert size.')
        alncnt = 0
        for alncnt, (query, alignments, found) in enumerate(
                _calc_primary_paired_alignment(insert2alignments, estimated_insertsize, minreasonable_insertsize,
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
                        handle.write(printstring)
            else:  # undefined
                continue

            if alncnt % 100000 == 0:
                logging.info('Paired inserts processed:\t{} / {}'.format(format(alncnt, ',d'), format(len(insert2alignments), ',d')))
        logging.info('Paired inserts processed:\t{} / {}'.format(format(alncnt, ',d'), format(len(insert2alignments), ',d')))
        logging.info('Singleton inserts tossed:\t{} / {}'.format(format(len(insert2alignments) - alncnt, ',d'), format(len(insert2alignments), ',d')))
        logging.info('Finished screening for OPRs and Paired-End inserts with unreasonable insert size.')


def align():
    x = 0

def main():
    startup()

    parser = argparse.ArgumentParser(
        description='A toolkit to align reads and find OPRs', usage='''oprfinder <command> [<args>]

    Command options
        align     align read reads using bwa
        find    find OPRs in alignment files
    ''')
    parser.add_argument('command', help='Subcommand to run: align|find')

    args = parser.parse_args(sys.argv[1:2])

    if args.command == 'align':
        align()
    elif args.command == 'find':
        find()
    else:
        print('Unrecognized command')
        parser.print_help()
        shutdown(1)
    shutdown(0)



if __name__ == '__main__':
    main()