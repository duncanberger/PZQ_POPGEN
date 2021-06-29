# Figure 2: *Schistosoma mansoni* population structure

0. [Setup](#setup)
1. [Figure 2A](#figure2a)
2. [Figure 2B](#figure2b)
3. [Figure 2C](#figure2c)
4. [Figure 2D](#figure2d)
5. [Figure 2E](#figure2e)
6. [Merged figure](#figure2f)
## Setup <a name="setup"></a>
```{r}
# Load required packages
library("ggplot2")
library("reshape2")
library("lattice")
library("cowplot")
library("SNPRelate")
library("ape")
library("ggtree")
library("phangorn")

# Load metadata
key <- read.table("supplementary_data_10.txt", header=TRUE, sep="\t", check.names = FALSE, comment.char = "")
```
## Figure 2A: Principal component analysis <a name="figure2a"></a>
```{r}
# Create a theme and palette
PCA_theme <- theme(axis.title=element_text(face="bold",size=9),
                   panel.background=element_blank(),
                   legend.text = element_text(face="bold"),
                   legend.position = "none",
                   axis.text=element_text(face="bold"),
                   panel.border = element_rect(color="#4c4c4c",fill=NA),
                   panel.grid=element_blank(),
                   legend.title=element_blank())
pca_palette <- c("#56B4E9", "#009e73","#CC79A7","#E69f00")

# Load data, produce in STEP X, and merge with metadata
eigenvec <- read.delim("207_PCA.eigenvec", header=TRUE, sep="\t")
eigenval <- read.delim("207_PCA.eigenval", sep="\t", header=FALSE)
eigenvec_merged <- (merge(key, eigenvec, all=TRUE, by.y = "IID", by.x='sample_ID'))

# Calculate contribution of each eigenvalue to total variance
sum_eigenval<-sum(eigenval$V1)
sum_eigenval<-lapply(eigenval$V1,function(x){
  rt<-(x/sum_eigenval)*100
  rt<-round(rt)
  return(rt) 
})

# Plot first two principal components, color points by school
pc1_pc2 <- ggplot(eigenvec_merged, aes((PC1),(PC2))) + 
  geom_point(size=0.75, aes(color=school, fill=school)) +
  xlab(paste0("PC1 (",sum_eigenval[[1]],"%)")) + 
  ylab(paste0("PC2 (",sum_eigenval[[2]],"%)")) + 
  scale_color_manual(values=pca_palette)+
  scale_fill_manual(values=pca_palette)+
  scale_shape_manual(values=c(21,23,22,24)) +
  scale_x_continuous(limits=c(-0.100000001,0.400000001), expand=c(0,0)) + 
  scale_y_continuous(limits=c(-0.200000001,0.200000001), expand=c(0,0)) +
  theme_bw() + PCA_theme 

# Plot first third and fourth principal components, color points by school
pc3_pc4 <- ggplot(eigenvec_merged, aes((PC3),(PC4))) + 
  geom_point(size=0.75, aes(color=school, fill=school)) +
  xlab(paste0("PC3 (",sum_eigenval[[3]],"%)")) + 
  ylab(paste0("PC4 (",sum_eigenval[[4]],"%)")) + 
  scale_color_manual(values=pca_palette)+
  scale_fill_manual(values=pca_palette)+
  scale_shape_manual(values=c(21,23,22,24)) +
  scale_x_continuous(limits=c(-0.400000001,0.600000001), expand=c(0,0)) + 
  scale_y_continuous(limits=c(-0.2000000001,0.800000001), expand=c(0,0)) +
  theme_bw() + PCA_theme
```
## Figure 2C: Mid-point rooted neighbour-joining tree <a name="figure2b"></a>
```{r}
# Load distance matrix
mdist <- as.matrix(read.table("autosomes.mdist", sep="\t", header=TRUE, row.names=1))

# Neighbor-joining tree estimation
nj_tree <- (nj(mdist))

# Plot tree
tree <- ggtree(nj_tree, layout="circular", aes(color=Site)) %<+% key +
  scale_color_manual(values=pca_palette, na.value='grey50')+
  geom_treescale(x=0.02, color='grey50', offset = 0.8,width = 0.025)
```
## Figure 2D: Autosomal nucleotide diversity values <a name="figure2c"></a>
```{r}
# Load nucleotide diversity results
pi_5kb_schools <- read.table("all.pi.pixy.schools.txt", header=TRUE)

# Order by school
pi_5kb_schools$pop = factor(pi_5kb_schools$pop, levels=c('Bugoto','Bwondha','Musubi','Kocoge'))

# Randomly subset dataset (for plotting)
pi_5kb_schools_subset <- sample_n(pi_5kb_schools, 10000)

# Plot figure
pi_all_ps <- ggplot(data=pi_5kb_schools_subset, aes(x=pop, y=log10(avg_pi), fill=pop, color=pop)) + 
  geom_point(position=position_jitterdodge(dodge.width =1,jitter.width = 0.8),alpha=0.3, size=0.0001) +
  geom_boxplot(aes(fill=pop),outlier.alpha = 0.0,notch = TRUE, outlier.colour = "grey35", color="black", alpha=0, width=0.275) +
  theme_bw()+
  scale_fill_manual(values=c("#56B4E9", "#009e73","#E69f00","#CC79A7")) +
  scale_color_manual(values=c("#56B4E9", "#009e73","#E69f00","#CC79A7")) +
  scale_y_continuous(expand=c(0,0), limits=c(-50,50000), breaks=c(-4.0,-3.0,-2.0,-1.0)) +
  coord_cartesian(ylim = c(-4, -1)) +
  xlab("A") +
  labs(y=expression(bold(-log[10]*("\U03C0")))) +
  PCA_theme + theme(legend.position = "none") +
  theme(legend.text = element_text(size=6.5, face = "bold"))
```
## Figure 2E: Pairwise comparisons of absolute (d<sub>XY</sub>) and relative (F<sub>ST</sub>) differentiation <a name="figure2d"></a>
### Calculate F<sub>ST</sub> estimates between each population
```{r}
# Load FST scores for all pairwise comparisons between schools
autosome_5kb_school_fst <- read.table("autosomes.fst.5kb.schools.txt", header=TRUE)

# Remove NA values and convert all negative values to 0
autosome_5kb_school_fst_1 <- subset(autosome_5kb_school_fst, avg_wc_fst!="NaN" & avg_wc_fst!="NA" & no_snps>125)
autosome_5kb_school_fst_1[autosome_5kb_school_fst_1 < 0] <- 0

# Subset for each combination
autosome_5kb_school_fst_1_Bugoto_Kocoge <- subset(autosome_5kb_school_fst_1, pop1=="Bugoto" & pop2=="Kocoge")
autosome_5kb_school_fst_1_Bwondha_Kocoge <- subset(autosome_5kb_school_fst_1, pop2=="Bwondha" & pop1=="Kocoge")
autosome_5kb_school_fst_1_Musubi_Kocoge <- subset(autosome_5kb_school_fst_1, pop1=="Musubi" & pop2=="Kocoge")
autosome_5kb_school_fst_1_Musubi_Bugoto <- subset(autosome_5kb_school_fst_1, pop1=="Musubi" & pop2=="Bugoto")
autosome_5kb_school_fst_1_Musubi_Bwondha <- subset(autosome_5kb_school_fst_1, pop1=="Musubi" & pop2=="Bwondha")
autosome_5kb_school_fst_1_Bugoto_Bwondha <- subset(autosome_5kb_school_fst_1, pop1=="Bugoto" & pop2=="Bwondha")

# Calculate mean and median (example for one comparison only)
median(autosome_5kb_school_fst_1_Bugoto_Bwondha$avg_wc_fst)
mean(autosome_5kb_school_fst_1_Bugoto_Bwondha$avg_wc_fst)
bstrap_means <- c()
bstrap_medians <- c()

# Calculate bootstrap medians, mean and quantiles (example for one comparison only)
for (i in 1:100) { 
  bstrap_medians <- c(bstrap_medians,median(sample(autosome_5kb_school_fst_1_Bugoto_Bwondha$avg_wc_fst,size=length(autosome_5kb_school_fst_1_Bugoto_Bwondha$avg_wc_fst),replace=TRUE)))
  bstrap_means <- c(bstrap_means,mean(sample(autosome_5kb_school_fst_1_Bugoto_Bwondha$avg_wc_fst,size=length(autosome_5kb_school_fst_1_Bugoto_Bwondha$avg_wc_fst),replace=TRUE)))
}
quantile(bstrap_medians,c(0.025,0.975))
```
### Calculate d<sub>XY</sub> estimates between each population
```{r}
# Load dXY scores for all pairwise comparisons between schools
autosome_5kb_schools_dxy <- read.table("autosomes.dxy.5kb.schools.txt", header=TRUE)

# Remove NA values
autosome_5kb_schools_dxy_1 <- subset(autosome_5kb_schools, avg_dxy!="NaN" & no_snps>125)

# Subset for each combination
autosome_5kb_schools_dxy_1_Bugoto_Kocoge <- subset(autosome_5kb_schools_dxy_1, pop1=="Bugoto" & pop2=="Kocoge")
autosome_5kb_schools_dxy_1_Bwondha_Kocoge <- subset(autosome_5kb_schools_dxy_1, pop2=="Bwondha" & pop1=="Kocoge")
autosome_5kb_schools_dxy_1_Musubi_Kocoge <- subset(autosome_5kb_schools_dxy_1, pop1=="Musubi" & pop2=="Kocoge")
autosome_5kb_schools_dxy_1_Musubi_Bugoto <- subset(autosome_5kb_schools_dxy_1, pop1=="Musubi" & pop2=="Bugoto")
autosome_5kb_schools_dxy_1_Musubi_Bwondha <- subset(autosome_5kb_schools_dxy_1, pop1=="Musubi" & pop2=="Bwondha")
autosome_5kb_schools_dxy_1_Bugoto_Bwondha <- subset(autosome_5kb_schools_dxy_1, pop1=="Bugoto" & pop2=="Bwondha")

# Calculate mean and median (example for one comparison only)
median(autosome_5kb_schools_dxy_1_Bugoto_Bwondha$avg_dxy)
mean(autosome_5kb_schools_dxy_1_Bugoto_Bwondha$avg_dxy)
bstrap_means <- c()
bstrap_medians <- c()

# Calculate bootstrap medians, mean and quantiles (example for one comparison only)
for (i in 1:100) { 
  bstrap_medians <- c(bstrap_medians,median(sample(autosome_5kb_schools_dxy_1_Bugoto_Bwondha$avg_dxy,size=length(autosome_5kb_schools_dxy_1_Bugoto_Bwondha$avg_dxy),replace=TRUE)))
  bstrap_means <- c(bstrap_means,mean(sample(autosome_5kb_schools_dxy_1_Bugoto_Bwondha$avg_dxy,size=length(aautosome_5kb_schools_dxy_1_Bugoto_Bwondha$avg_dxy),replace=TRUE)))
}
quantile(bstrap_medians,c(0.05,0.95))

# At present these values need to be manually added to figure, as ggplot/cowplot don't combine well with tables. 
```
## Figure 2F: Admixture <a name="figure2e"></a>
```{r}
# Set theme 
Admixture_theme <- theme(panel.grid=element_blank(), 
                         axis.ticks.x = element_blank(),
                         axis.text.y = element_text(face="bold", color="black"),
                         axis.text.x = element_blank(),
                         strip.text=element_text(face="bold"),
                         axis.title.y=element_text(face="bold",size=9),
                         panel.background = element_blank(),
                         strip.background = element_blank(),
                         panel.border = element_rect(color="black",fill=NA))

# Load data
admix <- read.table("admixture_all.txt", sep="\t", header=FALSE)

# Fix columns, merge with metadata and order by school
admix_summary <- (melt(admix,id.vars = c("V1","V2")))
admix_summary_merged <- (merge(key, admix_summary, all=TRUE, by.y = "V2", by.x='sample_ID'))
admix_summary_merged$g_order = factor(admix_summary_merged$school, levels=c('Bugoto','Bwondha','Musubi','Kocoge'))

#Plot for values of K (2-4)
admixture2 <- ggplot(data=subset(admix_summary_merged, V1==2)) +
  geom_bar(aes(x=sample_ID, y=as.numeric(value), fill=variable, color=variable), 
           width=1, show.legend=F,stat="identity")  + 
  facet_grid(.~g_order,scales="free_x", space = "free_x") + 
  xlab("") + ylab("") +
  scale_fill_manual(values=c("#dadaeb","#4a1486","#807dba","#4a1486","#D55E00","#f4a582","#F0E442","#016450")) +
  scale_color_manual(values=c("#dadaeb","#4a1486","#807dba","#4a1486","#D55E00","#f4a582","#F0E442","#016450")) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() + Admixture_theme
  
admixture3 <- ggplot(data=subset(admix_summary_merged, V1==3)) +
  geom_bar(aes(x=sample_ID, y=as.numeric(value), fill=variable, color=variable), 
           width=1, show.legend=F,stat="identity")  + 
  facet_grid(.~g_order,scales="free_x", space = "free_x") + 
  xlab("") + 
  ylab("") +
  scale_fill_manual(values=c("#807dba","#4a1486","#dadaeb")) +
  scale_color_manual(values=c("#807dba","#4a1486","#dadaeb")) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() + Admixture_theme

admixture4 <- ggplot(data=subset(admix_summary_merged, V1==4)) +
  geom_bar(aes(x=sample_ID, y=as.numeric(value), fill=variable, color=variable), 
           width=1, show.legend=F,stat="identity")  + 
  facet_grid(.~g_order,scales="free_x", space = "free_x") + 
  xlab("") + 
  ylab("") +
  scale_fill_manual(values=c("#9e9ac8","#807dba","#dadaeb","#4a1486")) +
  scale_color_manual(values=c("#9e9ac8","#807dba","#dadaeb","#4a1486")) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() + Admixture_theme
```
## Merged figure <a name="figure2f"></a>
```{r}
# Merge parts of figure together
top <-plot_grid(pc1_pc2,pc3_pc4, nrow=2, align="h", labels=c('A',''))
top_2 <- plot_grid(top, tree, nrow=1, rel_widths =c(0.75,1),rel_heights = c(1,2), labels=c('','B'))
middle <- plot_grid(pi_all_ps,"",nrow=1, labels=c('C','D'), rel_widths =c(0.75,1))
middle_2 <- plot_grid(top_2,middle,nrow=2,rel_heights = c(1,0.5), labels=c('F',''))
bottom<- plot_grid(admixture2,admixture3, admixture4, nrow=3, align="v",labels=c('E','',''))
plot_grid(middle_2,bottom, nrow=2, align="v", rel_heights =c(0.6,0.375),  rel_widths =c(1,0.6))
