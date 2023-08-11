from pysam import AlignedSegment

class Mapping:
    """
    Mapping object is a container for extracting and storing information about a single
    mapping of a read to a reference genome from SAM/BAM file.

    Attributes:
        rev: True if the read is mapped in reverse, False otherwise.
        ref: The reference name.
        rstart: The start position of the read on the reference.
        rend: The end position of the read on the reference.
        score: The alignment score.
        blocks: A list of start and end positions of aligned gapless blocks.
        cigartuples: The CIGAR alignment. The alignment is returned as a list of tuples of (operation, length).


    Methods:
        __init__: initialize the Mapping object.
        __repr__: represent the Mapping object.
        __str__: represent the Mapping object.


    Examples:
        >>> from pysam import AlignmentFile
        >>> from mVIRs.containers import Mapping
        >>> bam_file = AlignmentFile("tests/data/ERR4552622_100k_mVIRs.bam", "rb")
        >>> for read in bam_file:
        >>>     mapping = Mapping(read)
        >>>     print(mapping)
        Mapping(reverse=False, reference_name=NC_045512.2, map_start=1, map_end=29903, map_score=60, cirgar_tuples=[(0, 29903)], blocks=[(0, 29903)])
    """


    def __init__(self, rev, ref, rstart, rend, score,
                 blocks=None, cigartuples=None):
        """
        Initialize the Mapping object from pysam.AlignedSegment.

        Args:
            rev (bool): True if the read is mapped in reverse, False otherwise.
            ref (str): The reference name.
            rstart (int): The start position of the read on the reference.
            rend (int): The end position of the read on the reference.
            score (int): The alignment score.
            blocks (list): A list of start and end positions of aligned gapless blocks.
            cigartuples (list): The CIGAR alignment. The alignment is returned as a list of tuples of (operation, length).
        """

        self.rev: bool = rev
        self.ref: str = ref
        self.rstart: int = rstart
        self.rend: int = rend
        self.score: int = score
        self.blocks: list = blocks
        self.cigartuples: list = cigartuples

    @classmethod
    def from_aligned_segment(cls, alignment: AlignedSegment, extended: bool = True):
        rev: bool = True if alignment.is_reverse else False
        ref: str = alignment.reference_name
        rstart: int = alignment.reference_start
        rend: int = alignment.reference_end
        score: int = alignment.get_tag('AS')
        blocks = None
        cigartuples = None
        if extended:
            blocks = alignment.get_blocks()
            cigartuples = alignment.cigartuples

        return cls(rev, ref, rstart, rend, score,
                   blocks, cigartuples)

    def __repr__(self):
        return f'Mapping(reverse={self.rev}, reference_name={self.ref}, ' \
               f'map_start={self.rstart}, map_end={self.rend}, map_score={self.score}, '\
               f'cirgar_tuples={self.cigartuples}, blocks={self.blocks})'

    def __str__(self):
        return f'Mapping(reverse={self.rev}, reference_name={self.ref}, ' \
               f'map_start={self.rstart}, map_end={self.rend}, map_score={self.score}, '\
               f'cirgar_tuples={self.cigartuples}, blocks={self.blocks})'

class MappedRead:
    """
    MappedRead object is a container storing information about a single
    read from SAM/BAM file and all it's alignments to reference.

    Attributes:
        name: The name of the read.
        mapping: A dictionary of mappings for each read orientation.
    """

    def __init__(self, name):
        self.name = name
        self.mapping = {}

    def __setitem__(self, orientation: str, mapping: Mapping):
        # initiate the list if absent
        if orientation not in self.mapping:
            self.mapping[orientation] = []
        self.mapping[orientation].append(mapping)

    def __getitem__(self, orientation: str):
        return self.mapping[orientation]

    def __repr__(self) -> str:
        return f'MappedRead(name={self.name})'

    def __str__(self) -> str:
        return f'MappedRead(name={self.name})'

    def is_single_end(self) -> bool:
        """
        Check if the read is single-end.
        """
        return True if len(self.mapping) < 2 else False

    def total_mappings(self) -> int:
        """
        Return the total number of mappings for the read.
        """
        return sum([len(self.mapping[orientation]) \
                    for orientation in self.mapping])