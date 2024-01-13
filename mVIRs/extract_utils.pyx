from libc.stdlib cimport atoi
from libc.string cimport strtok, strcmp, strncmp, strlen, strcat, strcspn, strcpy, memset
from libc.stdio cimport FILE, fopen, fclose, getline
from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map as umap
from libcpp.string cimport string


cdef struct ClippedRead:
    char* name
    bint read_orientation
    bint clip_type
    bint clip_orientation
    char* clip_start
    char* scaffold


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
    read.scaffold[strlen(read.scaffold) - 1] = 0

    return read

cpdef tuple read_clipped_file(filename):
    cdef int header
    cdef char idx[400]
    memset(idx, 0, 400)

    cdef vector[ClippedRead] reads
    cdef umap[string, int] softclips
    cdef dict softclipped_positions = {}


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

    while True:
        read = getline(&line, &l, cfile)

        if read == -1:
            break

        read_obj = read_from_line(line)
        reads.push_back(read_obj)

        if read_obj.clip_type == 0:
            strcpy(idx, read_obj.clip_start)
            strcat(idx, b"_")
            strcat(idx, read_obj.scaffold)

            if softclips.find(idx) == softclips.end():
                softclips[idx] = 1
            else:
                softclips[idx] += 1

    fclose(cfile)

    softclipped_positions = dict(softclips)
    softclipped_positions = {tuple(k.split(b"_")) :v for k, v in softclipped_positions.items()}

    return reads, softclipped_positions