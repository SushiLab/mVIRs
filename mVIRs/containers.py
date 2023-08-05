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


    def __init__(self, alignment: AlignedSegment, extended: bool = False):
        """
        Initialize the Mapping object from pysam.AlignedSegment.
        """

        self.rev: bool = True if alignment.is_reverse else False
        self.ref: str = alignment.reference_name
        self.rstart: int = alignment.reference_start
        self.rend: int = alignment.reference_end
        self.score: int = alignment.get_tag('AS')
        self.blocks = None
        self.cigartuples = None
        if extended:
            self.blocks: list = alignment.get_blocks()
            self.cigartuples: list = alignment.cigartuples

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