
cimport cython
from libcpp cimport bool

try:
    import typing
    import dataclasses
except ImportError:
    pass

@cython.dataclasses.dataclass
cdef class AlignedSegmentInfo:
    cdef str readname
    cdef str orientation
    cdef int ascore
    cdef bool reverse
    cdef str refname 
    cdef int refstart
    cdef int refend
    cdef (int, int) blocks
    cdef (int, int) cigartuples


cpdef AlignedSegmentInfo extract_info_aligned_segment(
    pysam.AlignedSegment aligned_segment,
    bool need_extended = True):
    """
    Extracts information from a pysam.AlignedSegment object
    :param aligned_segment:
    :return: (reference_name, reverse, alignment_start, alignment_end, alignment_score)
    """
    data_tmp: List[str] = aligned_segment.qname.rsplit('/', 1)
    readname: str = data_tmp[0]
    orientation: str = data_tmp[1]
    ascore: int = aligned_segment.get_tag('AS')
    reverse: bool = True if aligned_segment.is_reverse else False
    refname: str = aligned_segment.reference_name
    refstart: int = aligned_segment.reference_start
    refend: int = aligned_segment.reference_end
    blocks = None
    cigartuples = None
    if need_extended:
        blocks: Tuple[int, int] = aligned_segment.get_blocks()[0]
        cigartuples: Tuple[int, int] = aligned_segment.cigartuples[0]
    
    return readname, orientation, ascore, reverse, refname, refstart, refend, blocks, cigartuples