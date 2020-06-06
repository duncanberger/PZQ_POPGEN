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
mkdir 00_SCRIPTS 01_REFERENCES 02_RAWDATA 03_MAPPING 04_VCALLING 05_QC 06_ANALYSIS
```

### Reference genomes
Get the reference genomes for mapping
```
cd 01_REFERENCES

# *S. mansoni reference genome*
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa.gz

# unzip
gunzip schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa.gz

# Select 
