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
...
```



