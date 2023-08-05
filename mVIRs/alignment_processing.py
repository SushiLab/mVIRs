from pysam.libcalignmentfile import AlignmentFile

from .containers import MappedRead, Mapping


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
        mapped[readname][orientation] = Mapping(alignment, extended=extended_info)

    alignments.close()

    return mapped



