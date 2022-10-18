




  
<p align="center">
<img src="https://raw.githubusercontent.com/SushiLab/mVIRs/master/pics/mVIRs.png" width="500"  />
</p>


# mVIRs: Localisation of inducible prophages using NGS data
 
mVIRs is a command-line tool that localizes and extracts genome sequences of inducible prophages in bacterial host genomes using paired-end DNA sequencing data as input. The approach relies on identifying DNA segments that are predicted to exist in a circularized or concatenated form upon induction. To achieve this, mVIRs uses information on the orientation of short, paired-end sequencing (e.g., Illumina) reads that are aligned to the genome of a lysogenic host as a reference. The identified segments can be length-filtered and classified by prediction tools (e.g., [VirSorter2](https://doi.org/10.1186/s40168-020-00990-y), [VirFinder](https://doi.org/10.1186/s40168-017-0283-5), [VIBRANT](https://doi.org/10.1186/s40168-020-00867-0) or [Prophage Hunter](https://doi.org/10.1093/nar/gkz380)), to identify putative prophage candidates.

The tool was developed by Hans-Joachim Ruscheweyh, Mirjam Zünd and Shinichi Sunagawa. It is distributed under [![License GPL v3](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html). 

If you use **mVIRs**, please cite:
> [Zünd M, Ruscheweyh HJ, Field CM, Meyer N, Cuenca M, Hoces D, Hardt WD, Sunagawa S. **High throughput sequencing provides exact genomic locations of inducible prophages and accurate phage-to-host ratios in gut microbial strains.** *Microbiome*, 9: 77, 2021.](https://doi.org/10.1186/s40168-021-01033-w)

Analyses in the publication were executed using version 1.0.0.

Questions/Comments? Write a github issue.


## Installation

The tools is written in `Python` and has the following dependencies:

- [Python](https://www.python.org/) >= 3.7
- [BWA](https://github.com/lh3/bwa) (tested: v0.7.17-r1188)
- [samtools](https://github.com/samtools/samtools) (tested: v1.9)


### Installation using conda

The easiest way to install **mVIRs** is to use the conda package manager using the [bioconda channel](https://bioconda.github.io/user/install.html#set-up-channels), which will automatically create an environment with the correct versions of the dependencies and then install **mVIRs** using `pip`.

```bash
# Install dependencies

$ conda create -n mvirs python==3.7 pip bwa samtools pysam -c bioconda 
$ conda activate mvirs
$ python -m pip install mvirs

# Test installation

$ mvirs -h

Program: mVIRs - Localisation of inducible prophages using NGS data
Version: 1.1.1
Reference: Zünd, Ruscheweyh, et al.
High throughput sequencing provides exact genomic locations of inducible
prophages and accurate phage-to-host ratios in gut microbial strains.
Microbiome (2021). doi:10.1186/s40168-021-01033-w

Usage: mvirs <command> [options]
Command:

    index   create index files for reference used in the
            mvirs oprs routine

    oprs    align reads against reference and used clipped
            alignment positions and OPRs to extract potential
            prophages

    test    run mVIRs for a public dataset
  
```


### Manual installation

Although installation using conda is recommended, manual installation of dependencies is also possible. `pip` is then used to install **mVIRs**:

```bash
# Manually install dependencies
...

# Install mVIRs
$ python -m pip install mvirs

# Test installation
$ mvirs -h

...

```

## Running mVIRs


The `mVIRs` toolkit includes three commands, `index`, `oprs` and `test`. The `index` command takes a reference sequence file as input and builds the reference database files that are needed for the execution of the `oprs` command. The `oprs` command aligns paired-end reads against the reference database to detect so called outward-oriented paired-end reads (OPRs) and uses soft-clipped alignments (clipped reads) to identify the location and extract the sequence of potential prophages. The test function executes `mVIRs` on a test dataset that is downloaded as part of this function.  

### `mvirs index`

This step takes a FASTA-formatted reference sequence file as input and builds an index using the `bwa index` command. This command needs to be executed before running the `oprs` command.

```bash
$ mvirs index

Program: mVIRs - Localisation of inducible prophages using NGS data
Version: 1.1.1
Reference: Zünd, Ruscheweyh, et al.
High throughput sequencing provides exact genomic locations of inducible
prophages and accurate phage-to-host ratios in gut microbial strains.
Microbiome, 9: 77, 2021. doi:10.1186/s40168-021-01033-w

Usage: mvirs index [options]

    Input:
        -f  FILE   Reference FASTA file. Can be gzipped. [Required]
        -o  DIR    Output directory. [Optional]

```

**Example**

```bash
mvirs index reference.fasta
```



### `mvirs oprs`

This step takes paired-end read files as input (one for each forward and reverse reads) and the name of the reference database produced by `mvirs index`. It aligns the reads against the reference database and uses the alignment information to identify potential prophage sequences within the reference genome using coverage information from OPRS and clipped reads.


```bash
$ mvirs oprs

Program: mVIRs - Localisation of inducible prophages using NGS data
Version: 1.1.1
Reference: Zünd, Ruscheweyh, et al.
High throughput sequencing provides exact genomic locations of inducible
prophages and accurate phage-to-host ratios in gut microbial strains.
Microbiome, 9: 77, 2021. doi:10.1186/s40168-021-01033-w

Usage: mvirs oprs [options]

    Input:
        -f  FILE   Forward reads file. FastA/Q. Can be gzipped. [Required]
        -r  FILE   Reverse reads file. FastA/Q. Can be gzipped. [Required]
        -db FILE   Reference database file (prefix) created by mvirs index. [Required]

    Output:
        -o  PATH   Prefix for output files. [Required]

    Options:
        -t  INT    Number of threads. [1] 
        -ml INT    Minimum sequence length for extraction. [4000]
        -ML INT    Maximum sequence length for extraction. [800000]
        -m         Allow full contigs/scaffolds/chromosomes to be reported 
	           (When OPRs and clipped reads are found at the start and end of contigs/scaffolds/chromosomes)

```

**Example**

Run `mvirs oprs` on the read files (`reads.1.fq.gz` and `reads.2.fq.gz`) using the same reference sequence file (`reference.fasta`) that was used as input for `mvirs index`.

```bash
$ mvirs oprs -f reads.1.fq.gz -r reads.2.fq.gz -db reference.fasta -o mvirs.output

# Will produce the following files (see below for explanation of the files)
# mvirs.output.bam --> Alignments
# mvirs.output.oprs --> The OPR positions in a tab separated file
# mvirs.output.clipped --> The clipped alignment positions in a tab separated file
# mvirs.output.fasta --> The potential prophage regions as a fasta file
```

### `mvirs test`

The `mvirs test` command downloads example read and reference files and launches the default `mvirs oprs` using the downloaded files as input.



```bash
$ mvirs.py test

2021-09-20 08:44:31,283 INFO: Starting mVIRs
Program: mVIRs - Localisation of inducible prophages using NGS data
Version: 1.1.1
Reference: Zünd, Ruscheweyh, et al.
High throughput sequencing provides exact genomic locations of inducible
prophages and accurate phage-to-host ratios in gut microbial strains.
Microbiome, 9: 77, 2021. doi:10.1186/s40168-021-01033-w
Usage: mvirs test [options]

    Input:
        -o  PATH   Output folder. [Required]
```

**Example**

```bash
$ mvirs test ~/mVIRs_test/
# Will produce the following files (see below for explanation of the files)
# mvirs.output.bam --> Alignments
# mvirs.output.oprs --> The OPR positions in a tab separated file
# mvirs.output.clipped --> The clipped alignment positions in a tab separated file
# mvirs.output.fasta --> The potential prophage regions as a fasta file

```


## Output files

`mvirs oprs` produces four output files (`mvirs.output.bam`. `mvirs.output.oprs`, `mvirs.output.clipped` and `mvirs.output.fasta`) of which we will explain the latter three.



### mvirs.output.fasta

The `mvirs.output.fasta` is a FASTA-formatted file with the potential prophage sequences that were extracted from the reference genome. The FASTA headers include information on the 

- source scaffold
- start and end coordinates of the extracted sequence
- number of supporting OPRs
- number of supporting clipped alignments
- fraction of the scaffold length that is covered by the extracted region

**Example**

```
>SalmonellaLT2:1213986-1255756	ORPs=3868-HSs=1473-SF=0.852597
ATTCGTAATGCGAAGGTCGTAGGTTCGACTCCTATTATCGGCACCAGTTAAATCAAATACTTAC...

# SalmonellaLT2:1213986-1255756 --> Scaffold:START-STOP
# ORPs=3868 --> Number of OPRs = 3868
# HSs=1473 --> Number of hard- and soft-clipped alignments
# SF=0.852597 --> 0.85% of the scaffold length is covered by the extracted region
```


### mvirs.output.oprs 

The `mvirs.output.oprs` file lists the inserts of paired-end reads that align either with an unusual orientation (e.g. OPR or SAME) or have an unexpected large insert size (IPR) when compared to the estimated insert size (See section [Concepts](#concepts) for the definition of OPR, SAME and IPR)

The columns of the file are:
   
- `INSERTNAME`: The name of the insert
- `REFERENCE`: Name of the scaffold/genomic region of the reference sequence
- `INSERT_SIZE`: The size of the insert
- `R1_ORIENTATION`: The orientation of the R1 read on the reference
- `R2_ORIENTATION`: The orientation of the R2 read on the reference
- `BWA_SCORE`: The sum of scores reported by bwa for the insert
- `R1_START`: The leftmost coordinate on the reference where the R1 read aligns
- `R2_START`: The leftmost coordinate on the reference where the R2 read aligns
- `R1_ALNLENGTH`: The length of the alignment of the R1 read on the reference
- `R2_ALNLENGTH`: The length of the alignment of the R2 read on the reference
- `INSERT_ORIENTATION`: The orientation of the both reads to each other. Can either be IPR, OPR or SAME


An example output is below:

```
#MIN_REASONABLE_INSERTSIZE=0
#MAX_REASONABLE_INSERTSIZE=1628
```
| #INSERTNAME                               	| REFERENCE             	| INSERT_SIZE 	| R1_ORIENTATION 	| R2_ORIENTATION 	| BWA_SCORE 	| R1_START 	| R2_START 	| R1_ALNLENGTH 	| R2_ALNLENGTH 	| INSERT_ORIENTATION 	|
|-----------------------------------------	|-----------------------	|-------------	|----------------	|----------------	|-----------	|----------	|----------	|--------------	|--------------	|--------------------	|
| K00206:180:H2CJWBBXY:8:1107:6644:49230  	| SalmonellaLT2         	| 41477       	| forward        	| reverse        	| 297       	| 1255437  	| 1214111  	| 151          	| 147          	| OPR                	|
| K00206:180:H2CJWBBXY:8:1107:7182:12181  	| SalmonellaLT2         	| 41392       	| forward        	| reverse        	| 288       	| 1255606  	| 1214365  	| 151          	| 143          	| OPR                	|
| K00206:180:H2CJWBBXY:8:1107:7436:46873  	| SalmonellaLT2         	| 41449       	| reverse        	| forward        	| 302       	| 1214126  	| 1255424  	| 151          	| 151          	| OPR                	|
| K00206:180:H2CJWBBXY:8:1107:8582:43304  	| SalmonellaLT2         	| 1351429     	| reverse        	| reverse        	| 225       	| 4216570  	| 2865291  	| 150          	| 80           	| SAME               	|
| K00206:180:H2CJWBBXY:8:1107:9404:2176   	| SalmonellaLT2         	| 41222       	| forward        	| reverse        	| 291       	| 1255124  	| 1214053  	| 151          	| 145          	| OPR                	|
| K00206:180:H2CJWBBXY:8:1107:10470:14959 	| SalmonellaLT2         	| 41453       	| reverse        	| forward        	| 302       	| 1214268  	| 1255570  	| 151          	| 151          	| OPR                	|
| K00206:180:H2CJWBBXY:8:1107:10724:29958 	| SalmonellaLT2         	| 140         	| reverse        	| forward        	| 201       	| 475072   	| 475135   	| 140          	| 76           	| OPR                	|

### mvirs.output.clipped


The `mvirs.output.clipped` file contains the name, orientation and position of aligned reads that were clipped (i.e., the read could only be partially aligned).

The columns of the file are:

- `INSERTNAME`: Name of the insert
- `READ ORIENTATION`: R1 or R2 
- `HARD/SOFTCLIP`: Reported clip type by BWA (Soft --> longer part of the alignment. Hard --> shorter part of the alignment)
- `DIRECTION`: Direction of the alignment with respect to the reference sequence
- `POSITION`: Leftmost coordinate of the alignment on the reference sequence
- `REFERENCE`: Name of the scaffold/genomic region of the reference sequence

An example output is below:

| #INSERT                                | READORIENTATION | HARD/SOFTCLIP | DIRECTION | POSITION | REFERENCE      |
|----------------------------------------|-----------------|---------------|-----------|----------|---------------|
| K00206:180:H2CJWBBXY:8:1101:1336:1982  | R2              | S             | ->        | 252072   | SalmonellaLT2 |
| K00206:180:H2CJWBBXY:8:1101:1418:9895  | R2              | S             | <-        | 1946641  | SalmonellaLT2 |
| K00206:180:H2CJWBBXY:8:1101:1468:42231 | R2              | S             | ->        | 1492581  | SalmonellaLT2 |
| K00206:180:H2CJWBBXY:8:1101:1864:3670  | R2              | S             | <-        | 2886283  | SalmonellaLT2 |
| K00206:180:H2CJWBBXY:8:1101:1915:46346 | R1              | S             | <-        | 1255756  | SalmonellaLT2 |
| K00206:180:H2CJWBBXY:8:1101:1915:46346 | R1              | H             | ->        | 1213986  | SalmonellaLT2 |
| K00206:180:H2CJWBBXY:8:1101:1996:25738 | R2              | S             | <-        | 4147832  | SalmonellaLT2 |
| K00206:180:H2CJWBBXY:8:1101:2037:2457  | R2              | S             | ->        | 465110   | SalmonellaLT2 |
| K00206:180:H2CJWBBXY:8:1101:2087:10422 | R1              | S             | <-        | 4769629  | SalmonellaLT2 |



  

## Concepts

### IPRs, OPRs and SAME


A paired-end read can align in the following orientations:

```
# IPRs --> If the insert orientation matches the reference sequence

REFERENCE ---------------------------------------
R1              -------->
R2                        <--------

# OPRs --> If the insert orientation does not match the reference genome (e.g., paired-end reads from circularized phage genomes spanning across attP sites)

REFERENCE ---------------------------------------
R1              <--------
R2                        -------->


# SAME -> If the insert orientation does not match the reference genome (e.g., due to inversions)

REFERENCE ---------------------------------------
R1              -------->
R2                        -------->
 
```


This tool reports IPRs with unreasonable insert sizes and OPRs.

## Algorithm

The algorithm for potential prophage genome detection consists of three steps.The first two steps scan the alignment file (`mvirs.output.bam`) and report OPRs and clipped alignments. The last step uses the information generated in the first two steps to detect potential prophages and report their sequences as a FASTA formatted output file. 


### 1. Read Alignment Orientation

1. Reads are conceptually paired as inserts according to their naming.
2. Upper and lower maxima for reasonable insert sizes are estimated using the mean insert size +/- seven standard deviations for uniquely mapping inward-oriented paired-end reads (IPRs).
3. For each insert:
  - Find the best-scoring alignment pairs within 3% of the best alignment score.
  - Report the insert as OPR if there is no IPR with a reasonable insert size within the 3% cutoff and if the OPR is the best scoring alignment.

### 2. Clipped Alignments

The alignment of a read is clipped if it can not be fully aligned against the reference genome. Two reasons for clipped alignments are:

1. The representation of circular bacterial chromosome sequence in linear form. Reads that align at the beginning of the linearised chromosome (here named reference) will also align at the end. A single full length alignment is not possible. In the example below, the first three bases of a given read align at the end of the reference, the last 3 bases at the beginning:

	```
	REFERENCE ---------------------------------------
	READ      ---                                 ---
	          456                                 123  
	```
		
	Similarly, if a read originates from a circularized phage genome that is also encoded as a prophage in the reference, the reads will align at the end and the beginning of the integrated prophage genome:  	
 
	```
	# Phage genome integrated in reference chromosome (denoted as P)
	REFERENCE ----------------PPPPPPPPPPP------------
	READ                      ---     ---
	                          456     123
	```
	

2. The reference genome contains an element (e.g., an integrated phage), but the read originates from a genome of a naive host (i.e., without the prophage). In the example below, the read originates from a naive host that does not contain the prophage genome (denoted as P), and thus flanks the phage integration site. 

	```
	REFERENCE ----------------PPPPPPPP---------------
	READ                   ---        ---
	                       123        456
	```



### 3. Identification of potential prophages

Regions where an accumulation of clipped alignments and OPRs are detected are reported as potential phages in the output fasta file.

The start and end positions of OPRs are distributed around the phage insertion sites. As such, they are indicative for potential prophages; however, they cannot be used to determine their exact positions. The locations of clipped alignments are precise; however, they often miss part of the alignment due to ambiguous alignments or not passing  the criteria for minimal alignment lengths. Therefore,  a genomic region requires the support from at least 1 OPR and at least 1 clipped alignment to be considered as a potential inducible  prophage. Furthermore, potential sites must have in total a sum of 5 or more OPRs and clipped alignments.

Potential prophage regions are also filtered by length, with a minimum requirement of 4kb and a maximum of 800kb. These limits can be modified with the `-ml` (minimum) and `-ML` (maximum) parameters.

If the `-m` flag is set, potential prophage regions that cover entire contigs/scaffolds will be reported, otherwise they will be discarded.


