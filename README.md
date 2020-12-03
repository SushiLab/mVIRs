# mVIRs - A bioinformatic approach to locate prophages by orientation and coverage in NGS data





## Overview


The tool was designed by Mirjam Zuend, Hans-Joachim Ruscheweyh and Shinichi Sunagawa and distributed under the GPLv3 license. 

Questions/Comments? Write a github issue.






## Installation



The tool is written in `python (>=3.6)` and requires `bwa` and `samtools` to be installed upfront. Installation of the `mVIRs` tool can be done with pip:

```bash
$git clone https://github.com/SushiLab/mVIRs

$cd mVIRs

#Installs the package locally

$pip install -r requirements.txt -e .

#Test

$ mvirs -h

2020-12-03 15:00:47,977 INFO: Starting mVIRs
usage: mvirs <command> [<args>]

    Command options
        oprs    align reads and find OPRs


Bioinformatic toolkit for finding prophages in sequencing data

positional arguments:
  command     Subcommand to run: oprs

optional arguments:
  -h, --help  show this help message and exit

```

## Usage


The `mVIRs` toolkit currently aligns and detects prophages with a single command: `oprs`


### OPRS


This step takes 2 FastQ files and a reference genome as bwa index as input, performs alignment, alignment filtering and detect OPRs:


```

usage: mvirs oprs [-h] [-t THREADS] i1 i2 r b o

Align reads against a reference and find OPRs and IPRs.

positional arguments:
  i1          Forward reads file
  i2          Reverse reads file
  r           BWA reference
  b           Output bam file
  o           Output OPR file

optional arguments:
  -h, --help  show this help message and exit
  -t THREADS  Number of threads to use. (Default = 1)

  
# The oprs subcommand can then be executed with the following command:

$ mvirs oprs r1.fq.gz r2.fq.gz reference.fasta output.bam output.oprs

```




## IPRs, OPRs and SAME


A paired-end read can align in the following orientations:

```
# IPRs --> When the insert orientation matches the reference sequence

REFERENCE ---------------------------------------
R1              -------->
R2                        <--------

# OPRs --> When the insert orientation doesnt match the reference genome, e.g. circularized virus particles

REFERENCE ---------------------------------------
R1              <--------
R2                        -------->


# SAME -> When the insert orientation doesnt match the reference genome, e.g. rearrangements

REFERENCE ---------------------------------------
R1              -------->
R2                        -------->
 
```


This tool reports IPRs with unreasonable insert sizes and OPRs.

#### Algorithm


The algorithm works the following:

1. Reads were grouped to inserts by name. High and low boundaries for insert sizes of properly paired inserts are estimated using the mean insert sizes and +/- 7 StDev of uniquely mapping inward-oriented paired-end reads (IPRs).
2. OPRs were detected the following: For each insert 
  - Find the best scoring alignment pairs within 3% of the best alignment score.
  - Report the insert as OPR if there is no IPR with reasonable insert size within the 3% cutoff and if the OPR is the best scoring alignment.




## The Output File


The output file from the `oprfinder find` script has the following columns:

- `MIN_REASONABLE_INSERTSIZE`: The low boundary for regular insert sizes
- `MAX_REASONABLE_INSERTSIZE`: The high boundary for regular insert sizes
- `INSERTNAME`: The name of the insert
- `REFERENCE`: The name of the reference
- `INSERT_SIZE`: The insert size of this insert
- `R1_ORIENTATION`: The orientation of the r1 read on the reference
- `R2_ORIENTATION`: The orientation of the r2 read on the reference
- `BWA_SCORE`: The additive score reported by bwa for this insert
- `R1_START`: The lowest coordinate on the reference where the r1 read aligns
- `R2_START`: The lowest coordinate on the reference where the r2 read aligns
- `R1_ALNLENGTH`: The length of the alignment of the r1 read on the reference
- `R2_ALNLENGTH`: The length of the alignment of the r2 read on the reference
- `READ_ORIENTATION`: The orientation of the both reads to each other. Can either be IPR, OPR or SAME


An example output is below:

```
#MIN_REASONABLE_INSERTSIZE=0
#MAX_REASONABLE_INSERTSIZE=867
```
**INSERTNAME**|**REFERENCE**|**INSERT\_SIZE**|**R1\_ORIENTATION**|**R2\_ORIENTATION**|**BWA\_SCORE**|**R1\_START**|**R2\_START**|**R1\_ALNLENGTH**|**R2\_ALNLENGTH**|**READ\_ORIENTATION**
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
A00917:31:HNNVNDSXX:2:1101:2031:3944|KB1|680277|forward|forward|300|549739|1229866|150|150|SAME
A00917:31:HNNVNDSXX:2:1101:2230:18067|KB1|2158899|reverse|reverse|300|132294|2291043|150|150|SAME
A00917:31:HNNVNDSXX:2:1101:2871:4773|KB1|936577|forward|forward|299|383083|1319511|150|149|SAME
A00917:31:HNNVNDSXX:2:1101:3884:4053|KB1|571933|reverse|forward|300|50395|622178|150|150|OPR
A00917:31:HNNVNDSXX:2:1101:5132:13197|KB1|3025889|forward|reverse|300|3025830|91|150|150|OPR
A00917:31:HNNVNDSXX:2:1101:5358:12242|KB1|3025889|forward|reverse|295|3025830|91|150|150|OPR
A00917:31:HNNVNDSXX:2:1101:5710:28604|KB1|406369|reverse|reverse|295|708899|302675|145|150|SAME
A00917:31:HNNVNDSXX:2:1101:8395:16939|KB1|361736|reverse|forward|297|684580|322991|147|150|IPR
A00917:31:HNNVNDSXX:2:1101:11071:18505|KB1|201078|reverse|forward|300|447787|246859|150|150|IPR
A00917:31:HNNVNDSXX:2:1101:12346:27571|KB1|2071244|forward|reverse|300|309857|2380951|150|150|IPR
A00917:31:HNNVNDSXX:2:1101:13593:10848|KB1|391337|forward|forward|300|1063656|1454843|150|150|SAME
A00917:31:HNNVNDSXX:2:1101:13829:12727|KB1|1414532|reverse|reverse|300|1886023|471641|150|150|SAME
A00917:31:HNNVNDSXX:2:1101:15619:8312|KB1|86962|forward|forward|295|1744410|1657598|150|150|SAME
A00917:31:HNNVNDSXX:2:1101:15845:8296|KB1|86962|forward|forward|295|1744410|1657598|150|150|SAME
A00917:31:HNNVNDSXX:2:1101:15845:31939|KB1|1279181|forward|reverse|295|1269110|2548141|150|150|IPR



 
