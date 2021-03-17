## Localisation of inducible prophages using NGS data

<p align="center">
<img src="pics/mVIRs.png" width="500"  />
</p>

mVIRs is a tool that locates integration sites of inducible prophages in bacterial genomes. It extracts information (i) on the alignment orientation of paired-end Illumina reads and (ii) locates partially aligned (soft/hard clipped) reads that are mapped to lysogenic host genome sequences to identify DNA segments that are predicted to exist in circularized form. These segments can be length-filtered and classified by prediction tools, such as VirSorter2, VirFinder, VIBRANT or Prophage Hunter, to identify putative prophage candidates.

The tool was developed by Mirjam Zuend, Hans-Joachim Ruscheweyh and Shinichi Sunagawa. It is distributed under [![License GPL v3](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html). 

If you use mVIRs, please cite:
> [Zünd M, Ruscheweyh HJ, Field CM, Meyer N, Cuenca M, Hoces D, Hardt WD, Sunagawa S. **High throughput sequencing provides exact genomic locations of inducible prophages and accurate phage-to-host ratios in gut microbial strains.** *Microbiome*, 2021.](https://doi.org/10.1186/s40168-021-01033-w)

Analyses in the publication were executed using version 1.0.

Questions/Comments? Write a github issue.

## Installation



The tool is written in `python (>=3.6.1)` and requires `bwa` (tested with v0.7.17-r1188) and `samtools` (tested with v1.9) to be installed upfront. Installation of the `mVIRs` tool can be done with pip:

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

usage: mvirs oprs [-h] [-t THREADS] i1 i2 r b o c f

Align reads against a reference and find OPRs and IPRs.

positional arguments:
  i1          Forward reads file
  i2          Reverse reads file
  r           BWA reference
  b           Output bam file
  o           Output OPR file
  c           Output Clipped file   
  f           Output Fasta file   

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


The output file from the `mvirs oprs` script has the following columns:

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
#MAX_REASONABLE_INSERTSIZE=1628
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
| K00206:180:H2CJWBBXY:8:1107:11464:41159 	| SalmonellaLT2         	| 42330       	| reverse        	| forward        	| 296       	| 1849504  	| 1891684  	| 151          	| 150          	| OPR                	|
| K00206:180:H2CJWBBXY:8:1107:12124:11354 	| SalmonellaLT2         	| 41504       	| reverse        	| forward        	| 293       	| 1214125  	| 1255487  	| 151          	| 142          	| OPR                	|
| K00206:180:H2CJWBBXY:8:1107:12195:32754 	| SalmonellaLT2         	| 41626       	| reverse        	| forward        	| 302       	| 1214118  	| 1255593  	| 151          	| 151          	| OPR                	|
| K00206:180:H2CJWBBXY:8:1107:14042:47630 	| SalmonellaLT2         	| 41508       	| forward        	| reverse        	| 292       	| 1255546  	| 1214189  	| 151          	| 146          	| OPR                	|
| K00206:180:H2CJWBBXY:8:1107:15077:21043 	| SalmonellaLT2_plasmid 	| 93601       	| reverse        	| forward        	| 270       	| 317      	| 93783    	| 145          	| 135          	| OPR                	|
| K00206:180:H2CJWBBXY:8:1107:15625:19953 	| SalmonellaLT2_plasmid 	| 93573       	| reverse        	| forward        	| 287       	| 59       	| 93481    	| 151          	| 151          	| OPR                	|
| K00206:180:H2CJWBBXY:8:1107:15777:9842  	| SalmonellaLT2_plasmid 	| 93693       	| reverse        	| forward        	| 275       	| 234      	| 93803    	| 151          	| 124          	| OPR                	|
| K00206:180:H2CJWBBXY:8:1112:7710:5675   	| SalmonellaLT2         	| 764617      	| forward        	| reverse        	| 258       	| 1243497  	| 2007992  	| 151          	| 122          	| IPR                	|
| K00206:180:H2CJWBBXY:8:1112:8308:5270   	| SalmonellaLT2         	| 729216      	| reverse        	| reverse        	| 276       	| 1667393  	| 938328   	| 151          	| 125          	| SAME               	|
| K00206:180:H2CJWBBXY:8:1112:8674:4813   	| SalmonellaLT2         	| 729220      	| reverse        	| reverse        	| 275       	| 1667393  	| 938324   	| 151          	| 129          	| SAME               	|
| K00206:180:H2CJWBBXY:8:1112:8937:34248  	| SalmonellaLT2         	| 41470       	| reverse        	| forward        	| 301       	| 1214052  	| 1255372  	| 151          	| 150          	| OPR                	|
| K00206:180:H2CJWBBXY:8:1112:7710:5675   	| SalmonellaLT2         	| 764617      	| forward        	| reverse        	| 258       	| 1243497  	| 2007992  	| 151          	| 122          	| IPR                	|
| K00206:180:H2CJWBBXY:8:1112:8308:5270   	| SalmonellaLT2         	| 729216      	| reverse        	| reverse        	| 276       	| 1667393  	| 938328   	| 151          	| 125          	| SAME               	|
| K00206:180:H2CJWBBXY:8:1112:8674:4813   	| SalmonellaLT2         	| 729220      	| reverse        	| reverse        	| 275       	| 1667393  	| 938324   	| 151          	| 129          	| SAME               	|
| K00206:180:H2CJWBBXY:8:1112:8937:34248  	| SalmonellaLT2         	| 41470       	| reverse        	| forward        	| 301       	| 1214052  	| 1255372  	| 151          	| 150          	| OPR                	|
| K00206:180:H2CJWBBXY:8:1112:10379:15715 	| SalmonellaLT2         	| 139         	| reverse        	| forward        	| 259       	| 3518002  	| 3518011  	| 139          	| 130          	| OPR                	|
| K00206:180:H2CJWBBXY:8:1112:10774:11020 	| SalmonellaLT2         	| 1336253     	| reverse        	| forward        	| 297       	| 2055194  	| 719092   	| 151          	| 151          	| IPR                	|



 
