from libc.string cimport strtok, strcmp, strncmp, strlen, strcat, strcspn, strcpy, memset
from libc.stdio cimport FILE, fopen, fclose, getline
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.string cimport string
from libcpp.pair cimport pair

ctypedef pair[string, string] idx

cdef struct ClippedRead:
    string name
    bint read_orientation
    bint clip_type
    bint clip_orientation
    string clip_start
    string scaffold

cdef class ClippedReads:
     cdef public vector[ClippedRead] _reads

cdef ClippedRead read_from_line(char * line) nogil:
    cdef char * delim = b"\t"
    cdef char * token
    cdef ClippedRead read

    # allocate memory
    token = strtok(line, delim)
    read.name = token
    token = strtok(NULL, delim)
    comp = strcmp(token, b"R1")

    if comp == 0:
        read.read_orientation = 0
    else:
        read.read_orientation = 1

    token = strtok(NULL, delim)
    comp = strcmp(token, b"S")

    if comp == 0:
        read.clip_type = 0
    else:
        read.clip_type = 1

    token = strtok(NULL, delim)
    comp = strcmp(token, b"->")

    if comp == 0:
        read.clip_orientation = 0
    else:
        read.clip_orientation = 1

    read.clip_start = strtok(NULL, delim)
    read.scaffold = strtok(NULL, delim)
    # remove trailng character in scaffold
    read.scaffold.erase(read.scaffold.length() - 1)

    return read

cpdef tuple read_clipped_file(str filename):
    cdef int header
    cdef idx map_idx
    cdef ClippedRead read_obj
    cdef vector[ClippedRead] reads
    cdef map[pair[string, string], int] softclips

    filename_byte_string = filename.encode("UTF-8")
    cdef char * fname = filename_byte_string
    cdef FILE * cfile
    cfile = fopen(fname, "rb")
    if cfile == NULL:
        raise FileNotFoundError(2, "No such file or directory: '%s'" % filename)

    cdef char * line = NULL
    cdef size_t l = 0
    cdef ssize_t read

    # skip header
    header = getline(&line, &l, cfile)

    with nogil:
        while True:
            read = getline(&line, &l, cfile)

            if read == -1:
                break

            read_obj = read_from_line(line)
            reads.push_back(read_obj)

            if read_obj.clip_type == 0:
                map_idx = idx(read_obj.scaffold, read_obj.clip_start)

                if softclips.find(map_idx) == softclips.end():
                    softclips[map_idx] = 1
                else:
                    softclips[map_idx] += 1

    fclose(cfile)

    cdef ClippedReads reads_list = ClippedReads.__new__(ClippedReads)
    reads_list._reads = reads

    return reads_list, softclips

# cdef umap denoise_softclips(umap softclipped_positions, int softclip_range):
#     cdef umap[string, int] denoised_softclips
#     cdef umap[string, int].iterator it


#     # iterate over map
#     with
