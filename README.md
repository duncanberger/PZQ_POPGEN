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
seqtk subseq schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa <(grep "Retained" ${WORKING_DIR}/00_METADATA/supplementary_table_3.txt | cut -f1 | cat) | cut -f1 -d " " > Sm_v7_nohap.fa

# Create indexes and a sequence dictionary for the reference genome
bwa index Sm_v7_nohap.fa

samtools faidx Sm_v7_nohap.fa

gatk CreateSequenceDictionary --REFERENCE Sm_v7_nohap.fa
```

### Sample metadata

## 02 - Mapping <a name="setup"></a>
### Download sequence reads
```
# Download FASTQ files in parallel
parallel -j4 --colsep '\t' "wget {1} {2}" :::: <(cat ${WORKING_DIR}/00_METADATA/supplementary_table_2.tx| cut -f13 | grep 'gz' | tr ';' '\t')

```

### Map sequence reads to reference genome
```
cd ${WORKING_DIR}/03_MAPPING

parallel --dry-run --colsep '\t' "bwa mem -t 6 Sm_v7_nohap.fa {8}_1.fastq.gz {8}_2.fastq.gz | samtools sort -@6  {1}.bam -" :::: <(cat ../00_METADATA/supplementary_table_2.txt | grep "gz")
```

The parallel command will write each mapping command to screen, which can be run individually or in batches. It will name the output BAM file with the sample name. 
```
# For example:
bwa mem -t 6 Sm_v7_nohap.fa ERR3173238_1.fastq.gz ERR3173238_2.fastq.gz | samtools sort -@6 -o PZQ_popgen6472766.bam -
```

### Mark PCR duplicates
```
parallel --dry-run --colsep '\t' "gatk MarkDuplicates --INPUT {1}.bam --OUTPUT {1}.markdup.bam --METRICS_FILE {1}.metrics.txt" :::: <(cat ../00_METADATA/supplementary_table_2.txt | grep "gz")
```
The parallel command will write each markduplicate command to screen, which can be run individually or in batches. It will name the output BAM file with the sample name. 
```
# For example:
gatk MarkDuplicates --INPUT PZQ_popgen6472766.bam --OUTPUT PZQ_popgen6472766.markdup.bam --METRICS_FILE PZQ_popgen6472766.metrics.txt
```






