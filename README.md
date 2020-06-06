# Impact of mass drug administration on the genomic diversity of Schistosoma mansoni populations

## Table of contents
0. [Project overview](#overview)
1. [Project setup](#setup)
2. [Mapping](#mapping)
3. [Variant calling](#variantcalling)
4. [Quality control](#qc)
5. [Analysis](#analysis)

## 00 - Project overview

The aim of this work is analyse the impact of mass drug adminstration of populations of *S. mansoni*.

## 01 - Project setup <a name="setup"></a>
### Setup a working environment for the analysis.
``` shell
mkdir ${HOME}/PZQ_POPGEN
cd ${HOME}/PZQ_POPGEN
WORKING_DIR=${PWD}

# make working directories
mkdir 00_METADATA 01_REFERENCES 02_RAWDATA 03_MAPPING 04_VCALLING 05_QC 06_ANALYSIS
```

### Reference genomes
Download the reference genomes
```
cd 01_REFERENCES

# Download the S. mansoni reference genome
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa.gz

# Unzip
gunzip schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa.gz

# Exclude haplotype scaffolds and trim scaffold names
seqtk subseq schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa <(grep "Retained" ../00_METADATA/supplementary_table_3.txt | cut -f1 | cat) | cut -f1 -d " " > Sm_v7_nohap.fa
```

### Sample metadata

## 02 - Mapping <a name="setup"></a>
### Download and map raw reads to the reference genome
```
# Make a list of FASTQ files to download
cut -f13 ${WORKING_DIR}/00_METADATA/supplementary_table_2.txt | grep gz | tr ';' '\t' > fastq.list

# Download FASTQ files
parallel -j4 --colsep '\t' "wget {1} {2}" :::: fastq.list

```
