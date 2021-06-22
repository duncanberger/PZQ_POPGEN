# Analysis 

## Table of contents
0. [Project setup](#setup)
1. [Population Structure](#structure)
2. [Selection](#selection)
3. [Association](#association)
4. [Population bottleneck](#bottleneck)

## 00 - Project setup <a name="setup"></a>
### Setup a working environment for the analysis.
```
WORKING_DIR=${PWD}
cd ${WORKING_DIR}/06_ANALYSIS

# make working directories
mkdir 00_SCRIPTS 01_STRUCTURE 02_SELECTION 03_ASSOCIATION
```
___
## 01 - Population structure <a name="structure"></a>
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
```
### Neighbour-joining phylogeny
```
# Produce a distance matrix
plink2 --bfile prunedData --allow-extra-chr --set-all-var-ids @_# --distance square 1-ibs
paste <( cut -f2 prunedData_tree.mdist.id) prunedData_tree.mdist | cat <(cut -f2 prunedData_tree.mdist.id | tr '\n' '\t' | sed -e '1s/^/\t/g') - > autosomes.mdist
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

# Get a table of CV scores, found in the stdout files (in our case *.o files)
cat *.o | grep CV | cut -f2 -d "=" | sed 's/)://g' | tr ' ' '\t' > cv_scores.txt
```
### Nucleotide diversity, F<sub>ST</sub> and d<sub>XY</sub>
```
# Create list for each level of population
cut -f1,2 ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep PZQ > host.list
cut -f1,3 ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep PZQ > school.list
cut -f1,4 ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep PZQ > district.list
cut -f1,3,6 ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep -v Kocoge | grep PZQ > treatment.list

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
cat pixy.SM_V7_*.5000.school_pi.txt | grep -v pop | cat pi.header - | grep -v ZW  > all.pi.pixy.schools.txt
cat pixy.SM_V7_*.5000.host_pi.txt | grep -v pop | cat pi.header - | grep -v ZW > pi.host.txt 
cat pixy.SM_V7_*.5000.treatment_pi.txt | grep -v pop | cat pi.header - | grep -v ZW > all.pi.treat.fix.txt
cat pixy.SM_V7_*.5000.school_fst.txt | grep -v pop | cat fst.header - | grep -v ZW > autosomes.fst.5kb.schools.txt
cat pixy.SM_V7_*.5000.school_dxy.txt | grep -v pop | cat dxy.header - | grep -v ZW  > autosomes.dxy.5kb.schools.txt 
cat pixy.SM_V7_*.5000.treatment_fst.txt | grep -v pop | cat fst.header - | grep -v ZW > fst.treatment.txt 
cat pixy.SM_V7_*.5000.treatment_dxy.txt | grep -v pop | cat dxy.header - | grep -v ZW > autosomes.dxy.5kb.treatment.txt 
cat pixy.SM_V7_*.5000.treatment_fst.txt | grep -v pop | cat fst.header - | grep -v ZW > autosomes.fst.5kb.treatment.txt
cat pixy.SM_V7_*.2000.treatment_fst.txt | grep -v pop | cat fst.header - | grep -v ZW > fst.treatment.2kb.txt 
cat pixy.SM_V7_*.5000.host_fst.txt | grep -v pop | cat fst.header | grep -v ZW > fst.host.txt 
```
### Recombination
```
# Create input files for each population
cut -f1,3 ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep PZQ | grep "Kocoge" > kocoge.list
cut -f1,3 ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep PZQ | grep -v "Kocoge" | shuf | head -17 > mayuge_1.list
cut -f1,3 ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep PZQ | grep -v "Kocoge" | grep -v -f mayuge_1.list | shuf | head -17 > mayuge_2.list
cut -f1,3 ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep PZQ | grep -v "Kocoge" | grep -v -f <(cat mayuge_1.list mayuge_2.list) | shuf | head -17 > mayuge_3.list
cut -f1,3 ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep PZQ | grep -v "Kocoge" | grep -v -f <(cat mayuge_1.list mayuge_2.list mayuge_3.list) | shuf | head -17 > mayuge_4.list
cut -f1,3 ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep PZQ | grep -v "Kocoge" | grep -v -f <(cat mayuge_1.list mayuge_2.list mayuge_3.list mayuge_4.list) | shuf | head -17 > mayuge_5.list

parallel --dry-run "vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FREEZE/PZQ_POPGEN.vcf --geno-r2 --out {.} --keep {} --min-r2 0.1 --ld-window-bp 50000 --maf 0.05" ::: kocoge.list mayuge_1.list mayuge_2.list mayuge_3.list mayuge_4.list mayuge_5.list

# Calculate median r2 for each distance
parallel "awk '{print $1,$3-$2,$5}' {}.geno.ld | grep -v CHR | tr ' ' '\t' | sort -k1,1 -k2,2 | datamash -g 1,2 median 3 > {}.median.txt" ::: mayuge_1 mayuge_2 mayuge_3 mayuge_4 mauge_5 kocoge

# Combine samples from Mayuge district
cat mayuge_*.median.txt > mayuge_median.ld
mv kocoge.median.txt kocoge_median.ld.txt
```
### Kinship
```
# Fix plink bim files
sed -i 's/SM_V7_//g' autosomes_unfiltered.bim

# Run King
king -bfile ${WORKING_DIR}/06_ANALYSIS/01_STRUCTURE/autosomes_unfiltered.bed --related --prefix autosomes_unfiltered
```
___
## 02 - Selection <a name="selection"></a>
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
### Create genetic maps
```
# Produce a list of variant sites
bcftools query -f '%CHROM\t%POS\n' PZQ_POPGEN.biallelic.vcf > PZQ_POPGEN.biallelic.txt

# Produce an approximate genetic map for each chromosome (using a per-chromosome recombination rate estimates). Commands vary depending on location of first variant ('$2-2389' etc.)
awk '{$3=((($2-2389)*3.509)/1000000)}{print $1,".",$3,$2}' <(grep "SM_V7_1" PZQ_POPGEN.biallelic.txt) > SM_V7_1.gmap
awk '{$3=((($2-7397)*4.975)/1000000)}{print $1,".",$3,$2}' <(grep "SM_V7_2" PZQ_POPGEN.biallelic.txt) > SM_V7_2.gmap
awk '{$3=((($2-19901)*3.236)/1000000)}{print $1,".",$3,$2}' <(grep "SM_V7_3" PZQ_POPGEN.biallelic.txt) > SM_V7_3.gmap
awk '{$3=((($2-1216)*4.103)/1000000)}{print $1,".",$3,$2}' <(grep "SM_V7_4" PZQ_POPGEN.biallelic.txt) > SM_V7_4.gmap
awk '{$3=((($2-14768)*3.910)/1000000)}{print $1,".",$3,$2}' <(grep "SM_V7_5" PZQ_POPGEN.biallelic.txt) > SM_V7_5.gmap
awk '{$3=((($2-51031)*4.623)/1000000)}{print $1,".",$3,$2}' <(grep "SM_V7_6" PZQ_POPGEN.biallelic.txt) > SM_V7_6.gmap
awk '{$3=((($2-19246)*5.979)/1000000)}{print $1,".",$3,$2}' <(grep "SM_V7_7" PZQ_POPGEN.biallelic.txt) > SM_V7_7.gmap
awk '{$3=((($2-1972)*4.031)/1000000)}{print $1,".",$3,$2}' <(grep "SM_V7_ZW" PZQ_POPGEN.biallelic.txt) > SM_V7_ZW.gmap
```
### Statistically phase variants using BEAGLE
```
# Phase variants for each vcf file using BEAGLE
parallel --dry-run "java -jar beagle.28Sep18.793.jar gt=PZQ_POPGEN.biallelic.{}.vcf out=PZQ_POPGEN.biallelic.{}.beagle map={}.gmap nthreads=4 iterations=100 burnin=10 ne=65000" ::: SM_V7_1 SM_V7_2 SM_V7_3 SM_V7_4 SM_V7_5 SM_V7_6 SM_V7_7 SM_V7_ZW
```
### Produce subsets for each population
```
cat ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep "gz" | grep 'Mayuge' | cut -f4 > mayuge.list
cat ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep "gz" | grep 'Tororo' | cut -f4 > tororo.list
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

cat *.mayuge.ihs.out.100bins.norm.temp | sed 's/SM_V7_//g' > ALL.MAYUGE.IHS.ihs.out.100bins.norm.txt
cat *.tororo.ihs.out.100bins.norm.temp | sed 's/SM_V7_//g' > ALL.TORORO.IHS.ihs.out.100bins.norm.txt
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

cat *.xpehh.out.norm.temp > ALL.MAYUGEvsTORORO.xpehh.xpehh.out.norm.txt
```
### Calculate F<sub>ST</sub> between districts (per chromosome)
```
# Run PIXY (described above too, the specific command is)
parallel --dry-run "pixy --stats fst 
--vcf PZQ_POPGEN.allsites.{1}.vcf 
--zarr_path zarr/ 
--window_size 2000 # Or 5kb where needed
--reuse_zarr yes
--populations district.list 
--variant_filter_expression 'DP>=10,GQ>=20,RGQ>=20' 
--invariant_filter_expression 'DP>=10,RGQ>=20' 
--outfile_prefix output/pixy.{1}.25000.district" ::: SM_V7_1 SM_V7_2 SM_V7_3 SM_V7_4 SM_V7_5 SM_V7_6 SM_V7_7 SM_V7_ZW

cat pixy.SM_V7_*.2000.district_fst.txt | grep -v pop | cat fst.header - > fst.district.txt 

# RUN VCFtools (alternative to PIXY for FST calculations) for between district comparisons

vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FREEZE/PZQ_POPGEN.vcf --weir-fst-pop MAYUGE.list --weir-fst-pop TORORO.list --fst-window-size 2000 --out MAYUGE_TORORO_2000.windowed.weir.txt

# RUN VCFtools (alternative to PIXY for FST calculations) for between treatment comparisons

vcftools --vcf PZQ_POPGEN.vcf --weir-fst-pop PRE_TREATMENT.list --weir-fst-pop GOOD_CLEARERS_MAYUGE.list --out PRE_GOOD_2kb --fst-window-size 2000
vcftools --vcf PZQ_POPGEN.vcf --weir-fst-pop POST_TREATMENT.list --weir-fst-pop GOOD_CLEARERS_MAYUGE.list --out POST_GOOD_2kb --fst-window-size 2000
vcftools --vcf PZQ_POPGEN.vcf --weir-fst-pop PRE_TREATMENT.list --weir-fst-pop POST_TREATMENT.list --out PRE_POST_2kb --fst-window-size 2000

sed 's/$/  POST_GOOD/g' POST_GOOD_2kb.windowed.weir.fst >> fst.windows.2kb.treatment.txt
sed 's/$/  PRE_GOOD/g' PRE_GOOD_2kb.windowed.weir.fst >> fst.windows.2kb.treatment.txt
sed 's/$/  PRE_POST/g' PRE_POST_2kb.windowed.weir.fst >> fst.windows.2kb.treatment.txt
```
___
## 03 - Association <a name="association"></a>
```
# Create phenotype and covariate files
cut -f1,4 /00_METADATA/supplementary_data_10.txt | grep PZQ | sed -e 's/^/0\t/g' | sed 's/Good_clearers/1/g' | sed 's/Post-treatment/2/g' | grep -v 'Pre' > bin_treatment.pheno
cut -f1,7 ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep PZQ | sed -e 's/^/0\t/g' > quant_ERR.pheno
cut -f2,3,4,5,6 ${WORKING_DIR}/06_ANALYSIS/01_STRUCTURE/prunedData.eigenvec | grep -v 'PC' | sed -e 's/^/0\t/g' > pca_covar_4.txt
```
### Binary trait association (Mayuge good vs post-treatment)
```
# Run the association test (use plink1.9)
plink --logistic --bfile ${WORKING_DIR}/06_ANALYSIS/01_STRUCTURE/prunedData --allow-extra-chr --allow-no-sex --out BIN_assoc_covar4_mayuge_maf --adjust --covar pca_covar_4.txt --pheno bin_treatment.pheno --set-missing-var-ids @[# --keep mayuge.list --chr 1,2,3,4,5,6,7 --maf 0.05 --keep treatment.list

# Fix output
cat BIN_assoc_covar4_mayuge_maf.logistic.adjusted | tr -s ' ' | sed 's/^ //g' | sed -e '2,$s/_/ /' | sed 's/SNP/CHR SNP/g' | tr ' ' '\t' > assoc_err_binary.txt
```
### Quantitative trait association (Using host ERR values)
```
# Run the association test (use plink1.9)
plink --linear --bfile ${WORKING_DIR}/06_ANALYSIS/01_STRUCTURE/prunedData --allow-extra-chr --allow-no-sex --out ERR_linear_covar4_mayuge_maf --adjust --covar pca_covar_4.txt --pheno quant_ERR.pheno --all-pheno --set-missing-var-ids @[# --keep treatment.list --maf 0.05 --chr 1,2,3,4,5,6,7

# Fix output
cat ERR_linear_covar4_mayuge_maf.linear.adjusted | tr -s ' ' | sed 's/^ //g' | sed -e '2,$s/_/ /' | sed 's/SNP/CHR SNP/g' | tr ' ' '\t' > assoc_err_linear.txt
```
## 04 - Population bottleneck <a name="bottleneck"></a>
### Site frequency spectra
```
# Get variant sites not in linkage disequilibrium (from plink filtering above)
sed 's/_/ /3' autosomes_unfiltered.prune.in > keep.LD.list

# Randomly subset to 200,000 variants 5x times.
cat keep.LD.list | shuf | head -200000 > keep.LD.A.list
cat keep.LD.list | grep -v -f keep.LD.A.list| shuf | head -200000 > keep.LD.B.list
cat keep.LD.list | grep -v -f keep.LD.A.list | grep -v -f keep.LD.B.list| shuf | head -200000 > keep.LD.C.list
cat keep.LD.list | grep -v -f keep.LD.A.list | grep -v -f keep.LD.B.list| grep -v -f keep.LD.C.list| shuf | head -200000 > keep.LD.D.list
cat keep.LD.list | grep -v -f keep.LD.A.list | grep -v -f keep.LD.B.list| grep -v -f keep.LD.C.list| grep -v -f keep.LD.D.list| shuf | head -200000 > keep.LD.E.list

# Subset VCF to contain only those subset variants, using the vcf that hasn't been filtered by MAF (repeat for subsets B-E), using a list of samples from each school (repeat for other schools). 
vcftools --vcf merged_all_samples.filtered.vcf.FL2.vcf --positions keep.LD.A.list --recode --recode-INFO-all --out LD_pruned.A.allMAF.vcf

# Identify best value for projecting down each population. Using a file of samples names (col 1) and source population (col 2)
easySFS.py -i LD_pruned.A.allMAF.vcf -p pops_file.txt --preview

# Pick best projection value that works for all populations 
easySFS.py -i LD_pruned.A.allMAF.vcf -p pop.list -f --proj 30,30,30,30 -a -o REP_A

# Repeat for each replicate (B-E) and collect SFS values
parallel "cat REP_{}/dadi/Kocoge-30.sfs | head -2 | tail -1 | tr ' ' '\n' | awk '{print \$1,NR,\"Kocoge\"}'" ::: A B C D E >> sfs_30.list
parallel "cat REP_{}/dadi/Bwondha-30.sfs | head -2 | tail -1 | tr ' ' '\n' | awk '{print \$1,NR,\"Bwondha\"}'" ::: A B C D E >> sfs_30.list
parallel "cat REP_{}/dadi/Bugoto-30.sfs | head -2 | tail -1 | tr ' ' '\n' | awk '{print \$1,NR,\"Bugoto\"}'" ::: A B C D E >> sfs_30.list
parallel "cat REP_{}/dadi/Musubi-30.sfs | head -2 | tail -1 | tr ' ' '\n' | awk '{print \$1,NR,\"Musubi\"}'" ::: A B C D E >> sfs_30.list
cat sfs_30.list | tr '\t' ',' > sfs.csv
# And collect residual values
parallel "cat REP_{}/dadi/Kocoge-30.sfs | head -3 | tail -1 | tr ' ' '\n' | awk '{print \$1,NR,\"Kocoge\"}'" ::: A B C D E >> sfs_30_res.list
parallel "cat REP_{}/dadi/Bwondha-30.sfs | head -3 | tail -1 | tr ' ' '\n' | awk '{print \$1,NR,\"Bwondha\"}'" ::: A B C D E >> sfs_30_res.list
parallel "cat REP_{}/dadi/Bugoto-30.sfs | head -3 | tail -1 | tr ' ' '\n' | awk '{print \$1,NR,\"Bugoto\"}'" ::: A B C D E >> sfs_30_res.list
parallel "cat REP_{}/dadi/Musubi-30.sfs | head -3 | tail -1 | tr ' ' '\n' | awk '{print \$1,NR,\"Musubi\"}'" ::: A B C D E >> sfs_30_res.list
cat sfs_30_res.list | tr '\t' ',' > sfs_res.csv
```
### Tajima's D
```
# Calculate Tajima's D for each school subpopulation in 2 kb windows. 
vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FREEZE/PZQ_POPGEN.vcf --keep KOCOGE.list --TajimaD 2000 --out KOCOGE_TAJIMA_D_2000
vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FREEZE/PZQ_POPGEN.vcf --keep MUSUBI.list --TajimaD 2000 --out MUSUBI_TAJIMA_D_2000
vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FREEZE/PZQ_POPGEN.vcf --keep BWONDHA.list --TajimaD 2000 --out BWONDHA_TAJIMA_D_2000
vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FREEZE/PZQ_POPGEN.vcf --keep BUGOTO.list --TajimaD 2000 --out BUGOTO_TAJIMA_D_2000
  
cat KOCOGE_TAJIMA_D_2000.Tajima.D | grep -v Taj | awk '{print $1,$2,$3,$4,"Kocoge"}' | sed 's/ / /g' >> KOCOGE_TAJIMA_D.Tajima.D.2kb.txt
cat BUGOTO_TAJIMA_D_2000.Tajima.D | grep -v Taj | awk '{print $1,$2,$3,$4,"Bugoto"}' | sed 's/ / /g' >> MAYUGE_TAJIMA_D.Tajima.D.2kb.txt
cat BWONDHA_TAJIMA_D_2000.Tajima.D | grep -v Taj | awk '{print $1,$2,$3,$4,"Bwondha"}' | sed 's/ / /g' >> MAYUGE_TAJIMA_D.Tajima.D.2kb.txt
cat MUSUBI_TAJIMA_D_2000.Tajima.D | grep -v Taj | awk '{print $1,$2,$3,$4,"Musubi"}' | sed 's/ / /g' >> MAYUGE_TAJIMA_D.Tajima.D.2kb.txt

  # Calculate Tajima's D for each school subpopulation in 5 kb windows. 
vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FREEZE/PZQ_POPGEN.vcf --keep KOCOGE.list --TajimaD 5000 --out KOCOGE_TAJIMA_D_5000
vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FREEZE/PZQ_POPGEN.vcf --keep MUSUBI.list --TajimaD 5000 --out MUSUBI_TAJIMA_D_5000
vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FREEZE/PZQ_POPGEN.vcf --keep BWONDHA.list --TajimaD 5000 --out BWONDHA_TAJIMA_D_5000
vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FREEZE/PZQ_POPGEN.vcf --keep BUGOTO.list --TajimaD 5000 --out BUGOTO_TAJIMA_D_5000
cat KOCOGE_TAJIMA_D_5000.Tajima.D | grep -v Taj | awk '{print $1,$2,$3,$4,"Kocoge"}' | sed 's/ / /g' >> TD.all.txt
cat BUGOTO_TAJIMA_D_5000.Tajima.D | grep -v Taj | awk '{print $1,$2,$3,$4,"Bugoto"}' | sed 's/ / /g' >> TD.all.txt
cat BWONDHA_TAJIMA_D_5000.Tajima.D | grep -v Taj | awk '{print $1,$2,$3,$4,"Bwondha"}' | sed 's/ / /g' >> TD.all.txt
cat MUSUBI_TAJIMA_D_5000.Tajima.D | grep -v Taj | awk '{print $1,$2,$3,$4,"Musubi"}' | sed 's/ / /g' >> TD.all.txt
```
