from libc.stdlib cimport atoi
from libc.stdio cimport printf
from libc.string cimport strtok
from libc.stdio cimport FILE, fopen, fclose, getline

ctypedef enum ReadOrientation:
    R1 = 0,
    R2 = 1

ctypedef enum ClipType:
    SOFT = 0,
    HARD = 1

ctypedef enum ClipOrientation:
    FORWARD = 0,
    REVERSE = 1


cdef class ClippedRead:
    cdef public char name
    cdef public ReadOrientation read_orientation
    cdef public ClipType clip_type
    cdef public ClipOrientation clip_orientation
    cdef public int clip_start
    cdef public char scaffold

    def __init__(self, name, read_orientation, clip_type, clip_orientation, clip_start, scaffold):
        self.name = name
        self.read_orientation = read_orientation
        self.clip_type = clip_type
        self.clip_orientation = clip_orientation
        self.clip_start = clip_start
        self.scaffold = scaffold

    def __str__(self):
        return "ClippedRead(%s, %s, %s, %s, %s, %s)" % (self.name, self.read_orientation, self.clip_type, self.clip_orientation, self.clip_start, self.scaffold)

    @property
    def read_id(self):
        # full id is name + read orientation
        return self.name + ("R1" if self.read_orientation == ReadOrientation.R1 else "R2")

# cdef class ClippedReadsList:
#     cdef public unordered_map[string, ClippedRead] reads


cpdef read_clipped_file(filename):
    cdef char* name
    cdef ReadOrientation read_orientation
    cdef ClipType clip_type
    cdef ClipOrientation clip_orientation
    cdef int clip_start
    cdef char* scaffold

    filename_byte_string = filename.encode("UTF-8")
    cdef char* fname = filename_byte_string
    cdef FILE* cfile
    cfile = fopen(fname, "rb")
    if cfile == NULL:
        raise FileNotFoundError(2, "No such file or directory: '%s'" % filename)

    cdef char * line = NULL
    cdef size_t l = 0
    cdef ssize_t read
    cdef char * delim = b"\t"
    cdef char * token


    while True:
        read = getline(&line, &l, cfile)
        if read == -1:
            break
        # skip header
        if line.startswith(b"#"):
            continue

        name = strtok(line, delim)
        clip_start       = atoi(strtok(NULL, delim))
        scaffold         = strtok(NULL, delim)
        read_orientation = int(0 if strtok(NULL, delim) == b"R2" else 1)
        clip_type        = int(0 if strtok(NULL, delim) == b"S" else 1)
        clip_orientation = int(0 if strtok(NULL, delim) == b"->" else 1)

        read = ClippedRead(name, read_orientation, clip_type, clip_orientation, clip_start, scaffold)
        print(read)

    fclose(cfile)

    # return clipped_reads_list