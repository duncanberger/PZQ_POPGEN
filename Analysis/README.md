# Analysis 

## Table of contents
0. [Project setup](#setup)
1. [Population Structure](#structure)
2. [Selection](#Selection)
3. [Association](#association)

## 00 - Project setup <a name="setup"></a>
### Setup a working environment for the analysis.
In my work environment, it is: ${WORKING_DIR}
```
cd ${HOME}/HELMINTH_EXTRACTION_WGS/
WORKING_DIR=${PWD}
cd ${WORKING_DIR}/06_ANALYSIS

# make working directories
mkdir 00_SCRIPTS 01_STRUCTURE 02_SELECTION 03_ASSOCIATION

#Kinship
#Recombination
#Structure - PCA, tree, admixture, pixy
#Selection - beagle, sweed, ihs, xpehh, fst
#Assocation
```
## 02 - Population structure <a name="setup"></a>
### Create plink files
```
# Convert vcf to plink .bed format, select autosomes. 
plink2 --vcf ${WORKING_DIR}/06_ANALYSIS/PZQ_POPGEN.vcf --chr SM_V7_1, SM_V7_2, SM_V7_3, SM_V7_4, SM_V7_5, SM_V7_6, SM_V7_7 --make-bed --allow-extra-chr --set-all-var-ids @_# --out autosomes_unfiltered

# Remove variants in strong linkage disequilibrium
plink2 --bfile autosomes_unfiltered --allow-extra-chr --set-all-var-ids @_# --indep-pairwise 50 10 0.15
plink2 --bfile autosomes_unfiltered --allow-extra-chr --set-all-var-ids @_# --extract autosomes_unfiltered.prune.in --out prunedData --make-bed
```
### Principal component analysis
```
plink2 --bfile prunedData --allow-extra-chr --set-all-var-ids @_# --pca

# Output can be passed to figure_2.R
```
### Neighbour-joining phylogeny
```
# Produce a distance matrix
plink2 --bfile prunedData --allow-extra-chr --set-all-var-ids @_# --distance square 1-ibs
paste <( cut -f2 prunedData_tree.mdist.id) prunedData_tree.mdist | cat <(cut -f2 prunedData_tree.mdist.id | tr '\n' '\t' | sed -e '1s/^/\t/g') - > autosomes.mdist

# Output can be passed to figure_2.R
```
### ADMIXTURE
```
#Fix scaffold names in bim file (ADMIXTURE accepts numerical scaffold names only)
sed -i 's/SM_V7_//g' prunedData.bim

# Produce random list of seeds
shuf -i 0-10000 | head -10 > seed.list

# Run ADMIXTURE of values of K:1-20, using 10 randomly generated seeds.
# There is no way of renaming admixture output based on seed value, so to avoid overwriting output files for each seed replicate, run in their own directory or batch run each seed one at a time. 
parallel --dry-run "admixture -j2 --seed={1} -B1000 prunedData.bed {2} --cv=10" :::: seed.list ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 > run_admix.sh

# For example:
mkdir 6127
cd 6127
admixture -j2 --seed=6127 -B1000 ../prunedData.bed 1 --cv=10
admixture -j2 --seed=6127 -B1000 ../prunedData.bed 2 --cv=10
# ... up to 20

# Add K values to ADMIXTURE output files
parallel --dry-run "sed -e 's/^/{1} /g' autosomes.{1}.Q" ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
cat *.Q > admixture_all.txt

# Output can be passed to figure_2.R
```
## 02 - Selection <a name="setup"></a>
### Create input VCFs
```
cd ${WORKING_DIR}/06_ANALYSIS/02_SELECTION
# Select only biallelic variants
vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/PZQ_POPGEN.vcf --recode --recode-INFO-all --out PZQ_POPGEN.biallelic.vcf --min-alleles 2 --max-alleles 2

# Split VCF into per-chromosome VCFs
parallel --dry-run "vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/PZQ_POPGEN.vcf --recode --recode-INFO-all --out PZQ_POPGEN.biallelic.{}.vcf --chr {}" ::: SM_V7_1 SM_V7_2 SM_V7_3 SM_V7_4 SM_V7_5 SM_V7_6 SM_V7_7 SM_V7_ZW

# E.g.
vcftools --vcf PZQ_POPGEN.biallelic.vcf --recode --recode-INFO-all --out PZQ_POPGEN.biallelic.SM_V7_1.vcf --chr SM_V7_1
vcftools --vcf PZQ_POPGEN.biallelic.vcf --recode --recode-INFO-all --out PZQ_POPGEN.biallelic.SM_V7_2.vcf --chr SM_V7_2
```

### Statistically phase variants using BEAGLE
```
# Phase variants for each vcf file using BEAGLE
parallel --dry-run "java -jar beagle.28Sep18.793.jar gt=PZQ_POPGEN.biallelic.{}.vcf out=PZQ_POPGEN.biallelic.{}.beagle map={}.gmap nthreads=4 iterations=100 burnin=10 ne=65000" ::: SM_V7_1 SM_V7_2 SM_V7_3 SM_V7_4 SM_V7_5 SM_V7_6 SM_V7_7 SM_V7_ZW
```
### Create genetic maps
```
INSERT METHOD HERE
```
### Produce subsets for each population
```
cat ${WORKING_DIR}/00_METADATA/supplementary_table_2.txt | grep "gz" | grep 'Mayuge' | cut -f4 > mayuge.list
cat ${WORKING_DIR}/00_METADATA/supplementary_table_2.txt | grep "gz" | grep 'Tororo' | cut -f4 > tororo.list
#EDIT KEEP MISSING VARIANTS
parallel --dry-run "vcftools --vcf PZQ_POPGEN.biallelic.{}.beagle.vcf.gz --recode --recode-INFO-all --out PZQ_POPGEN.biallelic.{1}.mayuge.vcf --chr {1} --keep {2}" ::: SM_V7_1 SM_V7_2 SM_V7_3 SM_V7_4 SM_V7_5 SM_V7_6 SM_V7_7 SM_V7_ZW :::: mayuge.list
parallel --dry-run "vcftools --vcf PZQ_POPGEN.biallelic.{}.beagle.vcf.gz --recode --recode-INFO-all --out PZQ_POPGEN.biallelic.{1}.tororo.vcf --chr {1} --keep {2}" ::: SM_V7_1 SM_V7_2 SM_V7_3 SM_V7_4 SM_V7_5 SM_V7_6 SM_V7_7 SM_V7_ZW ::: tororo.list
```
### Calculate IHS for each district (per chromosome)
```
# Calculate IHS
parallel --dry-run "selscan --ihs --vcf PZQ_POPGEN.biallelic.{1}.{2}.vcf --map {1}.gmap --threads 2 --out {1}.{2}.MAYUGE" ::: SM_V7_1 SM_V7_2 SM_V7_3 SM_V7_4 SM_V7_5 SM_V7_6 SM_V7_7 SM_V7_ZW ::: mayuge tororo

# Normalize
norm --ihs --files SM_V7_1.mayuge.ihs.out SM_V7_2.mayuge.ihs.out SM_V7_3.mayuge.ihs.out SM_V7_4.mayuge.ihs.out SM_V7_5.mayuge.ihs.out SM_V7_6.mayuge.ihs.out SM_V7_7.mayuge.ihs.out SM_V7_ZW.mayuge.ihs.out
norm --ihs --files SM_V7_1.tororo.ihs.out SM_V7_2.tororo.ihs.out SM_V7_3.tororo.ihs.out SM_V7_4.tororo.ihs.out SM_V7_5.tororo.ihs.out SM_V7_6.tororo.ihs.out SM_V7_7.tororo.ihs.out SM_V7_ZW.tororo.ihs.out

# Add chromosome names to normalized IHS results
parallel --dry-run "sed -e 's/^/{1}\t/g' {1}.{2}.ihs.out.100bins.norm > {1}.{2}.ihs.out.100bins.norm.temp" ::: SM_V7_1 SM_V7_2 SM_V7_3 SM_V7_4 SM_V7_5 SM_V7_6 SM_V7_7 SM_V7_ZW ::: mayuge tororo

cat *.mayuge.ihs.out.100bins.norm.temp | sed 's/SM_V7_//g' > mayuge.ihs.out.100bins.norm.all
cat *.tororo.ihs.out.100bins.norm.temp | sed 's/SM_V7_//g' > tororo.ihs.out.100bins.norm.all

# Output can be passed to figure_3.R
```
### Calculate XP-EHH between districts (per chromosome)
```
parallel --dry-run "selscan --xpehh --vcf PZQ_POPGEN.biallelic.{}.mayuge.vcf --vcf-ref PZQ_POPGEN.biallelic.{}.tororo.vcf --map {}.gmap --threads 2 --out {}" ::: SM_V7_1 SM_V7_2 SM_V7_3 SM_V7_4 SM_V7_5 SM_V7_6 SM_V7_7 SM_V7_ZW

#Normalize
norm --xpehh --files SM_V7_1.xpehh.out SM_V7_2.xpehh.out SM_V7_3.xpehh.out SM_V7_4.xpehh.out SM_V7_5.xpehh.out SM_V7_6.xpehh.out SM_V7_7.xpehh.out SM_V7_ZW.xpehh.out


# Create header for XP-EHH results
head -1 SM_V7_7.xpehh.out.norm | sed -e '1s/id/chromosome\tid/g' 
head -1 SM_V7_7.xpehh.out.norm | sed -e '1s/id/chromosome\tid/g'  > header.xpehh
# Add chromosome names to normalised XP-EHH results
parallel --dry-run "sed -e 's/^/{1}\t/g' {}.xpehh.out.norm > {}.xpehh.out.norm.temp" ::: SM_V7_1 SM_V7_2 SM_V7_3 SM_V7_4 SM_V7_5 SM_V7_6 SM_V7_7 SM_V7_ZW 

cat *.xpehh.out.norm.temp > xpehh.out.norm.all
# Output can be passed to figure_3.R
```
### Calculate CLR with districts (per chromosome)
```
```
### Calculate FST between districts (per chromosome)
```
```
