




  
<p align="center">
<img src="pics/mVIRs.png" width="500"  />
</p>


# mVIRs: Localisation of inducible prophages using NGS data

- [mVIRs: Localisation of inducible prophages using NGS data](#mvirs--localisation-of-inducible-prophages-using-ngs-data)
- [Installation](#installation)
  * [Installation using conda](#installation-using-conda)
  * [Manual installation](#manual-installation)
- [Usage](#usage)
  * [INDEX](#index)
  * [OPRS](#oprs)
- [Output files](#output-files)
    + [mvirs.output.oprs](#mvirsoutputoprs)
    + [mvirs.output.clipped](#mvirsoutputclipped)
    + [mvirs.output.fasta](#mvirsoutputfasta)
- [Concepts](#concepts)
  * [IPRs, OPRs and SAME](#iprs--oprs-and-same)
    + [Algorithm](#algorithm)
  * [Clipped Alignments](#clipped-alignments)
 
mVIRs is a tool that locates the integration sites of inducible prophages in bacterial genomes. It extracts information on the alignment orientation of paired-end Illumina reads that are mapped to lysogenic host reference genome. The aim is to identify DNA segments that are predicted to exist in circularized or concatenated form after induction. These segments can be length-filtered and classified by prediction tools, such as VirSorter2, VirFinder, VIBRANT or Prophage Hunter, to identify putative prophage candidates.

The tool was developed by Mirjam Zuend, Hans-Joachim Ruscheweyh and Shinichi Sunagawa. It is distributed under [![License GPL v3](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html). 

If you use mVIRs, please cite:
> [Zünd M, Ruscheweyh HJ, Field CM, Meyer N, Cuenca M, Hoces D, Hardt WD, Sunagawa S. **High throughput sequencing provides exact genomic locations of inducible prophages and accurate phage-to-host ratios in gut microbial strains.** *Microbiome*, 2021.](https://doi.org/10.1186/s40168-021-01033-w)

Analyses in the publication were executed using version 1.0.0.

Questions/Comments? Write a github issue.








# Installation

The tools is written in `python` and has the following dependencies:

- `Python>=3.7`
	- `pysam>=0.15.2`
- [BWA](https://github.com/lh3/bwa) (tested: v0.7.17-r1188)
- [samtools](https://github.com/samtools/samtools) (tested: v1.9)


## Installation using conda

The easiest way to install mVIRs is to use the conda package manager, which will automatically creat an environment with the correct versions of the dependencies.

```bash
# Install dependencies
$ git clone https://github.com/SushiLab/mVIRs
$ cd mVIRs   

$ conda env create -f conda_env_mvirs.yaml
$ conda activate mvirs   

#Installs the mVIRs package 
$ pip install -r requirements.txt -e .
# Add --user for installation with user permissions only


#Test
$ mvirs -h

2021-06-08 11:00:57,878 INFO: Starting mVIRs

Program: mVIRs - Localisation of inducible prophages using NGS data
Version: 1.1.0
Reference: Zuend, Ruscheweyh, et al.
High throughput sequencing provides exact genomic locations of inducible
prophages and accurate phage-to-host ratios in gut microbial strains.
Microbiome (2021). doi:10.1186/s40168-021-01033-w

Usage: mvirs <command> [options]
Command:

    oprs    align reads against reference and used clipped
            alignment positions and OPRs to extract potential
            prophages

    index   create index files for reference used in the
            mvirs oprs routine
            
    test    run mVIRs for a public dataset


mvirs.py: error: the following arguments are required: command
  
```


## Manual installation

Manual installation is possible but not recommended. Install via pip after installation of dependencies:

```bash
$ git clone https://github.com/SushiLab/mVIRs
$ cd mVIRs

#Installs the package locally
$ pip install -r requirements.txt -e .
# Add --user for installation with user permissions only

#Test
$ mvirs -h

Program: mVIRs - Localisation of inducible prophages using NGS data
Version: 1.1.0
Reference: Zuend, Ruscheweyh, et al.
High throughput sequencing provides exact genomic locations of inducible
prophages and accurate phage-to-host ratios in gut microbial strains.
Microbiome (2021). doi:10.1186/s40168-021-01033-w

Usage: mvirs <command> [options]
Command:

    oprs    align reads against reference and used clipped
            alignment positions and OPRs to extract potential
            prophages

    index   create index files for reference used in the
            mvirs oprs routine


mvirs.py: error: the following arguments are required: command
```

# Usage


The `mVIRs` toolkit includes 3 functions, `oprs`, `index` and `test`. The `index` commands takes a reference file as input and builds the index that is needed for the execution of the `oprs` command. The `oprs` command aligns paired-end reads against the reference database to detect so called outward orientated paired-end reads (OPRs) and uses soft-clipped alignments (clipped reads) to extract regions from the database that are potential prophages. The test function executes `mVIRs` on a test dataset that is downloaded as part of this function.  

## INDEX

This step takes the reference sequence file as input and builds an `bwa` index using the `bwa index` command. This command has to be executed before running the `oprs` command

```bash

Program: mVIRs - Localisation of inducible prophages using NGS data
Version: 1.1.0
Reference: Zuend, Ruscheweyh, et al.
High throughput sequencing provides exact genomic locations of inducible
prophages and accurate phage-to-host ratios in gut microbial strains.
Microbiome (2021). doi:10.1186/s40168-021-01033-w

Usage: mvirs index [options]

    Input:
        -f  FILE   Reference FastA file. Can be gzipped. [Required]

# The index subcommand can then be executed with the following command:

$ mvirs index -f reference.fasta
```


## OPRS


This step needs two sequence read files and a reference database (produced from `mvirs index`). It aligns the reads against the database and uses the alignments information to extract potential prophage sequences from the reference genome using coverage information from OPRS and Clipped reads.


```bash
$ mvirs oprs

Program: mVIRs - Localisation of inducible prophages using NGS data
Version: 1.1.0
Reference: Zuend, Ruscheweyh, et al.
High throughput sequencing provides exact genomic locations of inducible
prophages and accurate phage-to-host ratios in gut microbial strains.
Microbiome (2021). doi:10.1186/s40168-021-01033-w

Usage: mvirs oprs [options]

    Input:
        -f  FILE   Forward reads file. Can be gzipped. [Required]
        -r  FILE   Reverse reads file. Can be gzipped. [Required]
        -db FILE   BWA reference. Has to be created upfront. [Required]

    Output:
        -o  PATH   Prefix for output file. [Required]

    Options:
        -t  INT    Number of threads. [1] 
	-ml INT    Minimum length to extract. [4000]
        -ML INT    Maximum length to extract. [800000]
        -m         Allow full scaffolds to be reported 
	           (When OPRs and Clipped reads are found at start and end of scaffolds)

# Example
$ mvirs oprs -f reads.1.fq.gz -r reads.2.fq.gz -db reference.fasta -o mvirs.output


# Will produce the following files (see below for explanation of the files)
# mvirs.output.bam --> Alignments
# mvirs.output.oprs --> The raw OPRS positions
# mvirs.output.clipped --> The raw clipped alignment positions
# mvirs.output.fasta --> The potential prophage regions

```

## TEST

```bash
$ mvirs.py test

2021-09-20 08:44:31,283 INFO: Starting mVIRs
Program: mVIRs - Localisation of inducible prophages using NGS data
Version: 1.1.0
Reference: Zuend, Ruscheweyh, et al.
High throughput sequencing provides exact genomic locations of inducible
prophages and accurate phage-to-host ratios in gut microbial strains.
Microbiome (2021). doi:10.1186/s40168-021-01033-w
Usage: mvirs test [options]

    Input:
        -o  PATH   Output folder [Required]

# Example
$ mvirs test ~/mVIRs_test/
# Will produce the following files (see below for explanation of the files)
# ERR4552622_100k_mVIRs.bam --> Alignments
# ERR4552622_100k_mVIRs.oprs --> The raw OPRS positions
# ERR4552622_100k_mVIRs.clipped --> The raw clipped alignment positions
# ERR4552622_100k_mVIRs.fasta --> The potential prophage regions

```



# Output files

`mvirs oprs` produces 4 output files (`mvirs.output.bam`. `mvirs.output.oprs`, `mvirs.output.clipped` and `mvirs.output.fasta`) of which we will explain the latter 3.

### mvirs.output.oprs 

The `mvirs.output.oprs` file lists the inserts of paired-end reads that align either with an unusual orientation (e.g. OPR or SAME) or have an unexpected large insert size (IPR) when compared to the estimated insert size.

The columns of the file are:   

- `MIN_REASONABLE_INSERTSIZE`: The low boundary for regular insert sizes
- `MAX_REASONABLE_INSERTSIZE`: The high boundary for regular insert sizes
- `ESTIMATED_INSERTSIZE`: ??
- `READNAME`: The name of the insert
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
#MAX_REASONABLE_INSERTSIZE=1628
#ESTIMATED_INSERTSIZE=
```
| #READNAME                               	| REFERENCE             	| INSERT_SIZE 	| R1_ORIENTATION 	| R2_ORIENTATION 	| BWA_SCORE 	| R1_START 	| R2_START 	| R1_ALNLENGTH 	| R2_ALNLENGTH 	| INSERT_ORIENTATION 	|
|-----------------------------------------	|-----------------------	|-------------	|----------------	|----------------	|-----------	|----------	|----------	|--------------	|--------------	|--------------------	|
| K00206:180:H2CJWBBXY:8:1107:6644:49230  	| SalmonellaLT2         	| 41477       	| forward        	| reverse        	| 297       	| 1255437  	| 1214111  	| 151          	| 147          	| OPR                	|
| K00206:180:H2CJWBBXY:8:1107:7182:12181  	| SalmonellaLT2         	| 41392       	| forward        	| reverse        	| 288       	| 1255606  	| 1214365  	| 151          	| 143          	| OPR                	|
| K00206:180:H2CJWBBXY:8:1107:7436:46873  	| SalmonellaLT2         	| 41449       	| reverse        	| forward        	| 302       	| 1214126  	| 1255424  	| 151          	| 151          	| OPR                	|
| K00206:180:H2CJWBBXY:8:1107:8582:43304  	| SalmonellaLT2         	| 1351429     	| reverse        	| reverse        	| 225       	| 4216570  	| 2865291  	| 150          	| 80           	| SAME               	|
| K00206:180:H2CJWBBXY:8:1107:9404:2176   	| SalmonellaLT2         	| 41222       	| forward        	| reverse        	| 291       	| 1255124  	| 1214053  	| 151          	| 145          	| OPR                	|
| K00206:180:H2CJWBBXY:8:1107:10470:14959 	| SalmonellaLT2         	| 41453       	| reverse        	| forward        	| 302       	| 1214268  	| 1255570  	| 151          	| 151          	| OPR                	|
| K00206:180:H2CJWBBXY:8:1107:10724:29958 	| SalmonellaLT2         	| 140         	| reverse        	| forward        	| 201       	| 475072   	| 475135   	| 140          	| 76           	| OPR                	|

### mvirs.output.clipped


The `mvirs.output.clipped` file contains name, orientation and position of aligned reads that were clipped (An alignment is clipped when not the entire read can be mapped consecutively).

The columns of the file are:

- `Insert`: Name of the insert
- `READORIENTATION`: Orientation of thr read (R1 or R2)
- `HARD/SOFTCLIP`: Reported clip type by BWA (Soft --> longer part of the alignment. Hard --> shorter part of the alignment)
- `DIRECTION`: Direction of the alignment in respect to the reference sequence
- `POSITON`: Leftmost aligned based on the reference sequence
- `SCAFFOLD`: Name of the scaffold/genomic region of the reference sequence

An example output is below:

| #INSERT                                | READORIENTATION | HARD/SOFTCLIP | DIRECTION | POSITION | SCAFFOLD      |
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


### mvirs.output.fasta

The `mvirs.output.fasta` is a fasta file with the potential prophage regions that were extracted from the reference. The header of the fasta includes information on source scaffold, start and end coordinates of the potential prophage, number of supporting OPRS, the number of supporting clipped alignments and the fraction of the scaffold that is covered by the extracted region.

```
>SalmonellaLT2:1213986-1255756	ORPs=3868-HSs=1473-SF=0.852597
ATTCGTAATGCGAAGGTCGTAGGTTCGACTCCTATTATCGGCACCAGTTAAATCAAATACTTAC...
>SalmonellaLT2:1849457-1892188	ORPs=743-HSs=338-SF=0.872212
TCCTTTCAGTGATTGCATAACCACTTAACATCTTGTTTTATCTAAATAAAATTAAGCATGTTAT...
>SalmonellaLT2_plasmid:0-93927  OPRs=15-HSs=10-SF=100.000000
ACAGCGGCCGCCGAATTTCATACGCAGTCTTTCCCGCTGTTTTCTGGTCAGCGCCATAACATCT...

# SalmonellaLT2:1213986-1255756 --> Scaffold:START-STOP
# ORPs=3868 --> Number of OPRs = 3868
# HSs=1473 --> Number of clipped alignments at the location of the OPR = 1473
# SF=0.852597 --> 0.85% of the scaffold is covered by the extracted region
```
  

# Concepts

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

### Algorithm


The algorithm works the following:

1. Reads were grouped to inserts by name. Upper and lower maxima for insert sizes of paired inserts are estimated using the mean insert sizes and +/- 7 StDev's of uniquely mapping inward-oriented paired-end reads (IPRs).
2. OPRs are detected the following: For each insert 
  - Find the best scoring alignment pairs within 3% of the best alignment score.
  - Report the insert as OPR if there is no IPR with reasonable insert size within the 3% cutoff and if the OPR is the best scoring alignment.

## Clipped Alignments

The alignment of a read is clipped if it can't align against the reference in full length. Two reasons why alignments can be clipped are:

1. A circular chromosome is linearized. Reads that align to the beginning of the linearized chromosome (here named reference) will also align at the end. A single full length alignment is not possible. E.g. The first three bases of the read align at the end of the genome. The last 3 bases at the beginning

	```
	REFERENCE ---------------------------------------
	READ      ---                                 ---
	          456                                 123  
	```
		
	Or assuming a phage is integrated in the chromosome but the read comes from a the same, but circularized phage.  
 
	```
	# Phage integrated in reference chromosome denoted as P
	REFERENCE ----------------PPPPPPPPPPP------------
	READ                      ---     ---
	                          456     123
	```
	

2. The reference contains an element (e.g an integrated phage) that was not part of the genome from which the read was sequenced. This can happen when genomes from the same bacterial species contain or not contain an integrated phage. E.g here the bases denoted as P are from a phage that is absent in the read. 

	```
	REFERENCE ----------------PPPPPPPP---------------
	READ                   ---        ---
	                       123        456
	```

Regions where an accumulation of both cases and OPRs are detected are reported as potential phages in the output fasta file.

### Algorithm

OPRs have start and end positions but are not precise enough to detect exact excision positions. Clipped alignments are precise but often miss part of the alignment due to multimappers or the restriction on minimal alignment length. Therefore a region must have at least 1 OPR and at least 1 clipped alignment to be considered as a potential site for an integrated phage. Further, potential sites must have a sum of 5 or more OPRs and clipped alignments total.

Sites are also filtered by length, with a minimum requirement of 4kb and a maximum of 800kb. These limits can be modified with the -ml (minimum) and -ML (maximum) parameters on the command line.

If the -m flag is set, sites that cover entire scaffolds will be reported, otherwise they will be discounted.

2. Regions that cover an entire
