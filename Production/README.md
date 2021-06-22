# Impact of mass drug administration on the genomic diversity of Schistosoma mansoni populations

## Table of contents
0. [Project overview](#overview)
1. [Project setup](#setup)
2. [Raw data](#raw)
3. [Mapping](#mapping)
4. [Variant calling](#variantcalling)
5. [Quality control](#qc)

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
## 02 - Raw data <a name="raw"></a>
### Reference genome
Download the reference genome
```
cd 01_REFERENCES

# Download the S. mansoni reference genome
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa.gz

# Unzip
gunzip schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa.gz

# Exclude haplotype scaffolds and trim scaffold names
seqtk subseq schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa <(grep "Retained" ${WORKING_DIR}/00_METADATA/supplementary_data_11.txt | cut -f1 | cat) | cut -f1 -d " " > Sm_v7_nohap.fa

# Create indexes and a sequence dictionary for the reference genome
bwa index Sm_v7_nohap.fa

samtools faidx Sm_v7_nohap.fafsu

gatk CreateSequenceDictionary --REFERENCE Sm_v7_nohap.fa
```
## 03 - Mapping <a name="mapping"></a>
### Download sequence reads
```
# Download FASTQ files in parallel
parallel -j4 --colsep '\t' "wget {1} {2}" :::: <(cat ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | cut -f18 | grep 'gz' | tr ';' '\t')
```
### Map sequence reads to reference genome
```
cd ${WORKING_DIR}/03_MAPPING

parallel --dry-run --colsep '\t' "bwa mem -t 6 Sm_v7_nohap.fa {12}_1.fastq.gz {12}_2.fastq.gz | samtools sort -@6 {1}.bam -" :::: <(cat ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep "gz")

#The parallel command will write each mapping command to screen, which can be run individually or in batches. It will name the output BAM file with the sample name. For example:

bwa mem -t 6 ${WORKING_DIR}/01_REFERENCES/Sm_v7_nohap.fa ERR3173238_1.fastq.gz ERR3173238_2.fastq.gz | samtools sort -@6 -o PZQ_popgen6472766.bam -
```
### Mark PCR duplicates
```
parallel --dry-run --colsep '\t' "gatk MarkDuplicates --INPUT {1}.bam --OUTPUT {1}.markdup.bam --METRICS_FILE {1}.metrics.txt" :::: <(cat ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep "gz" | cut -f1)
```
The parallel command will write each markduplicate command to screen, which can be run individually or in batches. It will name the output BAM file with the sample name. For example:
```
gatk MarkDuplicates --INPUT PZQ_popgen6472766.bam --OUTPUT PZQ_popgen6472766.markdup.bam --METRICS_FILE PZQ_popgen6472766.metrics.txt

# Index all BAM files
parallel -j1 --colsep '\t' "samtools index {1}.markdup.bam" <(cat ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep "gz")
```
### Calculate coverage
```
# Create makewindows input
cut -f1,2 ${WORKING_DIR}/01_REFERENCES/Sm_v7_nohap.fa.fai > Sm_v7_nohap.txt

# Create 25 kb windows
bedtools makewindows -g Sm_v7_nohap.txt -w 25000 > Sm_v7_nohap.25kb.bed

# Calculate per-sample coverage
parallel "coverageBed -a Sm_v7_nohap.25kb.bed -b {1}.markdup.bam > {1}.cov" :::: <(cat ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep "gz")

# Name and merge output
parallel "awk '{print \$1,\$2,\$3,\$4,FILENAME}' {1}.cov | sed 's/.cov//g' | sed 's/SM_V7_//g' | sed 's/ZW/Z/g' > {1}.renamed.cov" :::: <(cat ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep "gz")
cat *.renamed.cov | datamash -g1,2,3 median 4 stdev 4 > median.coverage.txt
```
## 04 - Variant calling <a name="variantcalling"></a>

### Per-sampling variant calling
```
cd ${WORKING_DIR}/04_VCALLING

# Variant call each sample
parallel --dry-run "gatk HaplotypeCaller --emit-ref-confidence GVCF -I ${WORKING_DIR}/03_MAPPING/{1}.remarkdup.bam -R ${WORKING_DIR}/01_REFERENCES/Sm_v7_nohap.fa -O {1}.g.vcf" :::: <(cat ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep "gz")

#For example:
samtools index PZQ_popgen6472766
gatk HaplotypeCaller --emit-ref-confidence GVCF -I PZQ_popgen6472766.remarkdup.bam -R ${WORKING_DIR}/01_REFERENCES/Sm_v7_nohap.fa -O PZQ_popgen6472766.g.vcf
```
### Rename samples in each gVCF
```
parallel --dry-run --colsep '\t' "gatk RenameSampleInVcf --INPUT {1}.g.vcf --OUTPUT {1}.renamed.g.vcf --NEW_SAMPLE_NAME {1}" :::: <(cat ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep "gz")
```
### Combine all samples into a single gVCF and genotype
```
# Create and list of arguments for input into CombineGVCFs
parallel -j1 --colsep '\t' "echo '--INPUT {1}.renamed.g.vcf'" :::: <(cat ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep "gz" | cut -f1) > argument.list

# Combine gVCFs
gatk CombineGVCFs --arguments_file argument.list --reference ${WORKING_DIR}/01_REFERENCES/Sm_v7_nohap.fa --output merged_all_samples.g.vcf

# Genotype
gatk GenotypeGVCFs --reference ${WORKING_DIR}/01_REFERENCES/Sm_v7_nohap.fa --variant merged_all_samples.g.vcf --output merged_all_samples.vcf
```

## 05 - Quality control <a name="qc"></a>
### Calculate quality scores for all variant sites
```
# Produce a table of quality scores for each variant site
gatk VariantsToTable --variant merged_all_samples.vcf -F CHROM -F POS -F TYPE -F QD -F FS -F MQ -F MQRankSum -F ReadPosRankSum -F SOR -F InbreedingCoeff -R ${WORKING_DIR}/01_REFERENCES/Sm_v7_nohap.fa --output cohort.genotyped.txt

```
### Separate and filter SNPs
```
# Select SNPs
gatk SelectVariants -R ${WORKING_DIR}/01_REFERENCES/Sm_v7_nohap.fa --variant merged_all_samples.vcf --select-type-to-include SNP --output merged_all_samples.SNPs.vcf

# Tag low-quality SNPs
gatk VariantFiltration \
--filter-expression "ReadPosRankSum < -8.0" --filter-name "RPRS8" \
--filter-expression "QD < 2.0" --filter-name "QD2" \
--filter-expression "FS > 60.0" --filter-name "FS60" \
--filter-expression "MQ < 40.0" --filter-name "MQ40" \
--filter-expression "MQRankSum < -12.5" --filter-name "MQ12.5" \
--filter-expression "SOR > 3.0" --filter-name "SOR3" \
--variant merged_all_samples.SNPs.vcf \
-R ${WORKING_DIR}/01_REFERENCES/Sm_v7_nohap.fa  \
--output merged_all_samples.SNPs.tagged.vcf

# Remove low-quality sites
gatk SelectVariants -R ${WORKING_DIR}/01_REFERENCES/Sm_v7_nohap.fa --variant merged_all_samples.SNPs.tagged.vcf --exclude-filtered --output merged_all_samples.SNPs.filtered.vcf
```
### Separate and filter indels and mixed sites
```
# Select indels and mixed sites
gatk SelectVariants -R ${WORKING_DIR}/01_REFERENCES/Sm_v7_nohap.fa --variant merged_all_samples.vcf --select-type-to-include SNP --output merged_all_samples.indels_mixed.vcf

# Tag low-quality indels and mixed sites
/lustre/scratch118/infgen/team133/db22/software/gatk-4.1.0.0/gatk VariantFiltration \
--filter-expression "QD < 2.0" --filter-name "QD2" \
--filter-expression "FS > 200.0" --filter-name "FS200" \
--filter-expression "ReadPosRankSum < -20.0" --filter-name "RPRS20" \
--filter-expression "SOR > 10.0" --filter-name "SOR10" \
--variant merged_all_samples.indels_mixed.vcf \
-R ${WORKING_DIR}/01_REFERENCES/Sm_v7_nohap.fa  \
--output merged_all_samples.indels_mixed.tagged.vcf

# Remove low-quality sites
gatk SelectVariants -R ${WORKING_DIR}/01_REFERENCES/Sm_v7_nohap.fa --variant merged_all_samples.indels_mixed.tagged.vcf --exclude-filtered --output merged_all_samples.indels_mixed.filtered.vcf
```
### Recombine filtered variants
```
gatk MergeVcfs --INPUT merged_all_samples.SNPs.filtered.vcf --INPUT merged_all_samples.indels_mixed.filtered.vcf --OUTPUT merged_all_samples.filtered.vcf
```
### Remove low-quality samples and variants 
```
# Calculate per-individual missingness rate 
vcftools --vcf merged_all_samples.filtered.vcf.FL1.vcf --missing-indv --out missing_indv

# Filter out individuals with high rates of missing variant calls
awk '$6<0.55' missing_site.imiss | grep -v "MISS" | cut -f1  > retain.samples.list
vcftools --vcf merged_all_samples.filtered.vcf.FL1.vcf --keep retain.samples.list --recode-INFO-all --recode --out merged_all_samples.filtered.vcf.FL1.vcf  

# Calculate per-site missingness rate 
vcftools --vcf merged_all_samples.filtered.vcf --missing-site --out missing_site

# Filter out sites with high rates of missing variant calls
awk '$6<0.1' missing_site.lmiss | grep -v "MISS" | cut -f1,2 > retain.variants.list
vcftools --vcf merged_all_samples.filtered.vcf.FL1.vcf --postions retain.variants.list --recode-INFO-all --recode --out merged_all_samples.filtered.vcf.FL2.vcf 

# Calculate inbreeding coefficient
vcftools --vcf merged_all_samples.filtered.vcf.FL2.vcf --het --out ib

# Remove samples with excessively low inbreeding coefficient
awk '$5<-0.3' ib.het | grep -v "INDV" | cut -f1 > retain.IB.samples.list
vcftools --vcf merged_all_samples.filtered.vcf.FL2.vcf --keep retain.IB.samples.list --recode-INFO-all --recode --out merged_all_samples.filtered.vcf.FL3.vcf

# Remove sites with low minor allele frequency and exclude all variants not on chromosomes 1-7 or Z. 
vcftools --vcf merged_all_samples.filtered.vcf.FL3.vcf --recode --recode-INFO-all --maf 0.01 --out merged_all_samples.filtered.vcf.FL4.vcf
```
### Functionally annotate variant calls
```
# Using an existing SnpEFF database 
java -jar snpEff.jar Sm_v7.2 merged_all_samples.filtered.vcf.FL4.vcf > merged_all_samples.filtered.vcf.FL4.SNPEFF.vcf
```
### Produce an allsites VCF (for analyses with PIXY)
```
# Genotype the gVCF again include invariant sites (can be run for each chromosome with '-L' option). 
gatk GenotypeGVCFs --reference REF --variant merged_all_samples.g.vcf --include-non-variant-sites --output merged_all_samples.allsites.vcf

# Remove samples that were filtered out 
cat retain.samples.list retain.IB.samples.list | sort | uniq > all.retained.list
vcftools --vcf merged_all_samples.allsites.vcf --keep all.retained.list --recode-INFO-all --recode --out merged_all_samples.allsites.198.vcf
```
### Move final versions of VCFs to the analysis folder
```
mkdir ${WORKING_DIR}/06_ANALYSIS/FREEZE
mv merged_all_samples.filtered.vcf.FL4.SNPEFF.vcf ${WORKING_DIR}/06_ANALYSIS/FREEZE/PZQ_POPGEN.snpeff.vcf
mv merged_all_samples.filtered.vcf.FL4.vcf ${WORKING_DIR}/06_ANALYSIS/FREEZE/PZQ_POPGEN.vcf
mv merged_all_samples.allsites.198.vcf ${WORKING_DIR}/06_ANALYSIS/FREEZE/PZQ_POPGEN.allsites.vcf
```
