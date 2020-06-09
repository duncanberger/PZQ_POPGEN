# Analysis 

## Table of contents
0. [Project setup](#setup)
1. [Population Structure](#structure)
2. [Selection](#Selection)
3. [Association](#association)

## 00 - Project setup <a name="setup"></a>
### Setup a working environment for the analysis.
```
WORKING_DIR=${PWD}
cd ${WORKING_DIR}/06_ANALYSIS

# make working directories
mkdir 00_SCRIPTS 01_STRUCTURE 02_SELECTION 03_ASSOCIATION
```
___
## 01 - Population structure <a name="setup"></a>
### Create plink files
```
cd 01_STRUCTURE

# Convert vcf to plink .bed format, select autosomes. 
plink2 --vcf ${WORKING_DIR}/06_ANALYSIS/FREEZE/PZQ_POPGEN.vcf --chr SM_V7_1, SM_V7_2, SM_V7_3, SM_V7_4, SM_V7_5, SM_V7_6, SM_V7_7 --make-bed --allow-extra-chr --set-all-var-ids @_# --out autosomes_unfiltered

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
### Admixture
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

# Get a table of CV scores, found in the stdout files (in our case *.o files)
cat *.o | grep CV | cut -f2 -d "=" | sed 's/)://g' | tr ' ' '\t' > cv_scores.txt

# Output can be passed to supplementary_figure_2.R
```
### Nucleotide diversity, F<sub>ST</sub> and d<sub>XY</sub>
```
# Create list for each level of population
cut -f1,2 ${WORKING_DIR}/00_METADATA/supplementary_table_2.txt | grep PZQ > host.list
cut -f1,3 ${WORKING_DIR}/00_METADATA/supplementary_table_2.txt | grep PZQ > school.list
cut -f1,4 ${WORKING_DIR}/00_METADATA/supplementary_table_2.txt | grep PZQ > district.list
cut -f1,3,6 ${WORKING_DIR}/00_METADATA/supplementary_table_2.txt | grep -v Kocoge | grep PZQ > treatment.list

# Subset the allsites VCF for each chromosome
parallel --dry-run "vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FREEZE/PZQ_POPGEN.allsites.vcf --chr {} --recode-INFO-all --recode --out PZQ_POPGEN.allsites.{}.vcf" ::: SM_V7_1 SM_V7_2 SM_V7_3 SM_V7_4 SM_V7_5 SM_V7_6 SM_V7_7 SM_V7_ZW

# Run PIXY (this is an example which will write command for calculate the 3 statistics for each population and chromosome)
parallel --dry-run "pixy --stats fst dxy pi 
--vcf ${WORKING_DIR}/06_ANALYSIS/FREEZE/PZQ_POPGEN.allsites.{1}.vcf 
--zarr_path zarr/ 
--window_size {2}
#--reuse_zarr yes # Add after the first run has been completed
--populations {3}.list # Repeat for each of the 3 other lists
--variant_filter_expression 'DP>=10,GQ>=20,RGQ>=20' 
--invariant_filter_expression 'DP>=10,RGQ>=20' 
--outfile_prefix output/pixy.{1}.{2}.{3}" ::: SM_V7_1 SM_V7_2 SM_V7_3 SM_V7_4 SM_V7_5 SM_V7_6 SM_V7_7 SM_V7_ZW ::: 5000 25000 ::: host school district treatment

# Create headers for each statistic (doesn't matter which file you use, as long as it's specific to each statistic)
head -1 pixy.SM_V7_1.5000.school_fst.txt > fst.header
head -1 pixy.SM_V7_1.5000.school_dxy.txt > dxy.header
head -1 pixy.SM_V7_1.5000.school_pi.txt > pi.header

# Combine files for each chromosome
cat pixy.SM_V7_*.5000.school_pi.txt | grep -v pop | cat pi.header - > pi.school.txt # Output can be passed to figure_2.R
cat pixy.SM_V7_*.5000.host_pi.txt | grep -v pop | cat pi.header - > pi.host.txt # Output can be passed to supplementary_figure_5.R
cat pixy.SM_V7_*.5000.treatment_pi.txt | grep -v pop | cat pi.header - > pi.treatment.txt # Output can be passed to figure_4.R
cat pixy.SM_V7_*.5000.school_fst.txt | grep -v pop | cat fst.header - > fst.school.txt # Output can be passed to figure_2.R
cat pixy.SM_V7_*.5000.school_dxy.txt | grep -v pop | cat dxy.header - > .school.txt # Output can be passed to figure_2.R
cat pixy.SM_V7_*.5000.treatment_fst.txt | grep -v pop | cat fst.header - > fst.treatment.txt # Output can be passed to figure_4.R
cat pixy.SM_V7_*.5000.treatment_dxy.txt | grep -v pop | cat dxy.header - > dxy.treatment.txt # Output can be passed to figure_4.R
```
### Recombination
```
# Create input files for each population
cut -f1,3 ${WORKING_DIR}/00_METADATA/supplementary_table_2.txt | grep PZQ | grep "Kocoge" > kocoge.list
cut -f1,3 ${WORKING_DIR}/00_METADATA/supplementary_table_2.txt | grep PZQ | grep -v "Kocoge" | shuf | head -17 > mayuge_1.list
cut -f1,3 ${WORKING_DIR}/00_METADATA/supplementary_table_2.txt | grep PZQ | grep -v "Kocoge" | grep -v -f mayuge_1.list | shuf | head -17 > mayuge_2.list
cut -f1,3 ${WORKING_DIR}/00_METADATA/supplementary_table_2.txt | grep PZQ | grep -v "Kocoge" | grep -v -f <(cat mayuge_1.list mayuge_2.list) | shuf | head -17 > mayuge_3.list
cut -f1,3 ${WORKING_DIR}/00_METADATA/supplementary_table_2.txt | grep PZQ | grep -v "Kocoge" | grep -v -f <(cat mayuge_1.list mayuge_2.list mayuge_3.list) | shuf | head -17 > mayuge_4.list
cut -f1,3 ${WORKING_DIR}/00_METADATA/supplementary_table_2.txt | grep PZQ | grep -v "Kocoge" | grep -v -f <(cat mayuge_1.list mayuge_2.list mayuge_3.list mayuge_4.list) | shuf | head -17 > mayuge_5.list

parallel --dry-run "vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FREEZE/PZQ_POPGEN.vcf --geno-r2 --out {.} --keep {} --min-r2 0.1 --ld-window-bp 50000 --maf 0.05" ::: kocoge.list mayuge_1.list mayuge_2.list mayuge_3.list mayuge_4.list mayuge_5.list

# Calculate median r2 for each distance
parallel "awk '{print $1,$3-$2,$5}' {}.geno.ld | grep -v CHR | tr ' ' '\t' | sort -k1,1 -k2,2 | datamash -g 1,2 median 3 > {}.median.txt" ::: mayuge_1 mayuge_2 mayuge_3 mayuge_4 mauge_5 kocoge

# Combine samples from Mayuge district
cat mayuge_*.median.txt > all.median.mayuge.geno.ld
mv kocoge.median.txt all.median.kocoge.geno.ld

# Output can be passed to supplementary_figure_6.R
```
### Kinship
```
# Fix plink bim files
sed -i 's/SM_V7_//g' autosomes_unfiltered.bim

# Run King
king -bfile ${WORKING_DIR}/06_ANALYSIS/01_STRUCTURE/autosomes_unfiltered.bed --related --prefix autosomes_unfiltered
```
___
## 02 - Selection <a name="setup"></a>
### Create input VCFs
```
cd ${WORKING_DIR}/06_ANALYSIS/02_SELECTION
# Select only biallelic variants
vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FREEZE/PZQ_POPGEN.vcf --recode --recode-INFO-all --out PZQ_POPGEN.biallelic.vcf --min-alleles 2 --max-alleles 2

# Split VCF into per-chromosome VCFs
parallel --dry-run "vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FREEZE/PZQ_POPGEN.vcf --recode --recode-INFO-all --out PZQ_POPGEN.biallelic.{}.vcf --chr {}" ::: SM_V7_1 SM_V7_2 SM_V7_3 SM_V7_4 SM_V7_5 SM_V7_6 SM_V7_7 SM_V7_ZW

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
# Calculate CLR per district, per chromosome
parallel --dry-run "SweeD-P -name {1}.{2} -input PZQ_POPGEN.biallelic.{1}.{2}.vcf -grid 20000" ::: SM_V7_1 SM_V7_2 SM_V7_3 SM_V7_4 SM_V7_5 SM_V7_6 SM_V7_7 SM_V7_ZW ::: mayuge tororo

# Add chromosome name to output
parallel --dry-run "grep -v '/' SweeD_Report.{1}.{2} | grep -v 'al' -i | grep -v '^$' | sed -e 's/^/{1}\t/g' > {1}.{2}.swd.temp" ::: SM_V7_1 ::: mayuge tororo

# Combine output for each district
cat *.mayuge.swd.temp > all.mayuge.swd
cat *.tororo.swd.temp > all.tororo.swd

# Output can be passed to figure_3.R
```
### Calculate FST between districts (per chromosome)
```
# Run PIXY (described above too, the specific command is)
parallel --dry-run "pixy --stats fst 
--vcf PZQ_POPGEN.allsites.{1}.vcf 
--zarr_path zarr/ 
--window_size 25000
--reuse_zarr yes
--populations district.list 
--variant_filter_expression 'DP>=10,GQ>=20,RGQ>=20' 
--invariant_filter_expression 'DP>=10,RGQ>=20' 
--outfile_prefix output/pixy.{1}.25000.district" ::: SM_V7_1 SM_V7_2 SM_V7_3 SM_V7_4 SM_V7_5 SM_V7_6 SM_V7_7 SM_V7_ZW

cat pixy.SM_V7_*.25000.district_fst.txt | grep -v pop | cat fst.header - > fst.district.txt 

# Output can be passed to figure_3.R
```
___
## 03 - Association <a name="setup"></a>
```
# Create phenotype and covariate files
cut -f1,4 /00_METADATA/supplementary_table_2.txt | grep PZQ | sed -e 's/^/0\t/g' | sed 's/Good_clearers/1/g' | sed 's/Post-treatment/2/g' | grep -v 'Pre' > bin_treatment.pheno
cut -f1,7 ${WORKING_DIR}/00_METADATA/supplementary_table_2.txt | grep PZQ | sed -e 's/^/0\t/g' > quant_ERR.pheno
cut -f2,3,4,5,6 ${WORKING_DIR}/06_ANALYSIS/01_STRUCTURE/prunedData.eigenvec | grep -v 'PC' | sed -e 's/^/0\t/g' > pca_covar_4.txt
```
### Binary trait association (Mayuge good vs post-treatment)
```
# Run the association test (use plink1.9)
plink --logistic --bfile ${WORKING_DIR}/06_ANALYSIS/01_STRUCTURE/prunedData --allow-extra-chr --allow-no-sex --out BIN_assoc_covar4_mayuge_maf --adjust --covar pca_covar_4.txt --pheno bin_treatment.pheno --set-missing-var-ids @[# --keep mayuge.list --chr 1,2,3,4,5,6,7 --maf 0.05 --keep treatment.list

# Fix output
cat BIN_assoc_covar4_mayuge_maf.logistic.adjusted | tr -s ' ' | sed 's/^ //g' | sed -e '2,$s/_/ /' | sed 's/SNP/CHR SNP/g' | tr ' ' '\t' > BIN_assoc_covar4_mayuge_maf.logistic.adjusted.tbl

# Output can be passed to figure_4.R
```
### Quantitative trait association (Using host ERR values)
```
# Run the association test (use plink1.9)
plink --linear --bfile ${WORKING_DIR}/06_ANALYSIS/01_STRUCTURE/prunedData --allow-extra-chr --allow-no-sex --out ERR_linear_covar4_mayuge_maf --adjust --covar pca_covar_4.txt --pheno quant_ERR.pheno --all-pheno --set-missing-var-ids @[# --keep treatment.list --maf 0.05 --chr 1,2,3,4,5,6,7

# Fix output
cat ERR_linear_covar4_mayuge_maf.linear.adjusted | tr -s ' ' | sed 's/^ //g' | sed -e '2,$s/_/ /' | sed 's/SNP/CHR SNP/g' | tr ' ' '\t' > ERR_linear_covar4_mayuge_maf.linear.adjusted.tbl

# Output can be passed to figure_4.R
```
