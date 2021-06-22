# Figure 5: The impact of a single round of praziquantel treatment on Mayuge district *Schistosoma mansoni* populations. 

0. [Setup](#setup)
1. [Figure 5A](#figure5a)
2. [Figure 5B](#figure5b)
3. [Figure 5C](#figure5c)
4. [Figure 5D](#figure5d)
5. [Figure 5E](#figure5e)
6. [Merged figure](#merge)
## Setup <a name="setup"></a>
```{r}
# Load required packages
library("ggplot2")
library("reshape2")
library("lattice")
library("cowplot")
library("scales")
library("dplyr")

# Load metadata
key <- read.table("supplementary_data_10.txt", header=TRUE, sep="\t", check.names = FALSE, comment.char = "")
scaleFUN <- function(x) sprintf("%.1f", x)

```
## Figure 5A: Nucleotide diversity estimates <a name="figure5a"></a>
```{r}
#Create theme
pi_theme <- theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  legend.title = element_blank(),
                  panel.grid = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.y = element_text(face="bold", color="black"),
                  axis.text.x = element_text(face="bold", color="black"),
                  axis.title.y = element_text(face="bold", color="black"),
                  axis.ticks.x=element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(color="black"),
                  legend.text = element_text(face="bold", color="black"))

# Load data
treatment_5kb_pi <- read.table("all.pi.treat.fix.txt", header=TRUE, sep='\t')

# Order by treatment stage
treatment_5kb_pi$pop = factor(treatment_5kb_pi$pop, levels=c('Good Clearers','Poor clearers: Pre-treatment','Poor clearers: Post-treatment'))

# Subset for plotting
treatment_5kb_pi_2 <- sample_n(treatment_5kb_pi, 15000)

# Plot boxplot
pi_all_a <- ggplot(data=treatment_5kb_pi_2, aes(x=pop, y=log10(avg_pi), fill=pop, color=pop)) + 
  geom_point(position=position_jitterdodge(dodge.width =1,jitter.width = 0.8),alpha=0.3, size=0.0001) +
  geom_boxplot(aes(fill=avg_pi),outlier.alpha = 0.0,notch = TRUE, outlier.colour = "grey35", color="black", alpha=0, width=0.325) +
  theme_bw() +
  scale_fill_manual(values=c("#785EF0","#DC267F","#FE6100","#02818a")) +
  scale_color_manual(values=c("#785EF0","#DC267F","#FE6100","#02818a")) +
  scale_y_continuous(expand=c(0,0), breaks=c(-4.00,-3.00,-2.00,-1.00),labels=scaleFUN) +
  coord_cartesian(ylim = c(-4, -1)) +
  labs(y=expression(bold(-log[10]*("\U03C0")))) +
  PCA_theme + theme(legend.position = "none") +
  xlab("") + 
  theme(legend.text = element_text(size=6.5, face = "bold"),
        axis.text = element_text(color = "black"))

# Plot density plot
pi_all_b <- ggplot(data=treatment_5kb_pi_2, aes(x=log10(avg_pi), fill=pop, color=pop)) +
  geom_density(alpha=0.3) + 
  theme_nothing()+
  xlab("") + ylab("") +
  scale_fill_manual(values=c("#785EF0","#DC267F","#FE6100","#02818a")) +
  scale_color_manual(values=c("#785EF0","#DC267F","#FE6100","#02818a")) +
  scale_y_continuous(expand=c(0,0), breaks=c(-4.00,-3.00,-2.00,-1.00),labels=scaleFUN) +
  PCA_theme + theme(legend.position = "none") +
  coord_cartesian(xlim = c(-4, -1)) +
  theme(legend.text = element_text(face = "bold"),
        axis.text =element_blank(), panel.border = element_blank(), axis.ticks = element_blank()) + 
  coord_flip()

# Combine plots
pi_all <- plot_grid(pi_all_a,pi_all_b, nrow=1, align="h", rel_widths = c(1,0.2))

```
## Figure 5B:  <a name="figure5b"></a>
### Calculate median F<sub>ST</sub> 
```{r}
# Read in data
autosome_5kb_treatment_fst <- read.table("autosomes.fst.5kb.treatment.txt", header=TRUE)

# Remove NA values and set all negative values to zero
autosome_5kb_treatment_fst_1 <- subset(autosome_5kb_treatment_fst, avg_wc_fst!="NaN" & avg_wc_fst!="NA")
autosome_5kb_treatment_fst_1[autosome_5kb_fst_1 < 0] <- 0

# Subset each target population
autosome_5kb_treatment_fst_1_POST_PRE <- subset(autosome_5kb_fst_1, pop1=="Post-treatment" & pop2=="Pre-treatment")
autosome_5kb_treatment_fst_1_Good_PRE <- subset(autosome_5kb_fst_1, pop1=="Good_clearers" & pop2=="Pre-treatment")
autosome_5kb_treatment_fst_1_Good_POST <- subset(autosome_5kb_fst_1, pop1=="Good_clearers" & pop2=="Post-treatment")

# Calculate mean and median FST values for each treatment group (example for one comparison only)
median(autosome_5kb_treatment_fst_1_Good_POST$avg_wc_fst)
mean(autosome_5kb_treatment_fst_1_Good_POST$avg_wc_fst)
bstrap_means <- c()
bstrap_medians <- c()

# Calculate bootstrap medians, mean and quantiles (example for one comparison only)
for (i in 1:100) { 
  bstrap_medians <- c(bstrap_medians,median(sample(autosome_5kb_treatment_fst_1_Good_POST$avg_wc_fst,size=length(autosome_5kb_treatment_fst_1_Good_POST$avg_wc_fst),replace=TRUE)))
  bstrap_means <- c(bstrap_means,mean(sample(autosome_5kb_treatment_fst_1_Good_POST$avg_wc_fst,size=length(autosome_5kb_treatment_fst_1_Good_POST$avg_wc_fst),replace=TRUE)))
}
quantile(bstrap_medians,c(0.025,0.975))
```
### Calculate median d<sub>XY</sub>
```{r}
# Read in data
autosome_5kb_treatment_dxy <- read.table("autosomes.dxy.5kb.treatment.txt", header=TRUE)
# Remove NA values and set all negative values to zero
autosome_5kb_treatment_dxy_1 <- subset(autosome_5kb_dxy, avg_dxy!="NaN")

# Subset each target population
autosome_5kb_treatment_dxy_1_POST_PRE <- subset(autosome_5kb_dxy_1, pop1=="Post-treatment" & pop2=="Pre-treatment")
autosome_5kb_treatment_dxy_1_Good_PRE <- subset(autosome_5kb_dxy_1, pop1=="Good_clearers" & pop2=="Pre-treatment")
autosome_5kb_treatment_dxy_1_Good_POST <- subset(autosome_5kb_dxy_1, pop1=="Good_clearers" & pop2=="Post-treatment")

# Calculate mean and median dXY values for each treatment group (example for one comparison only)
median(autosome_5kb_treatment_dxy_1_Good_POST$avg_dxy)
mean(autosome_5kb_treatment_dxy_1_Good_POST$avg_dxy)
bstrap_means <- c()
bstrap_medians <- c()

# Calculate bootstrap medians, mean and quantiles (example for one comparison only)
for (i in 1:100) { 
  bstrap_medians <- c(bstrap_medians,median(sample(autosome_5kb_treatment_dxy_1_Good_POST$avg_dxy,size=length(autosome_5kb_treatment_dxy_1_Good_POST$avg_dxy),replace=TRUE)))
  bstrap_means <- c(bstrap_means,mean(sample(autosome_5kb_treatment_dxy_1_Good_POST$avg_dxy,size=length(autosome_5kb_treatment_dxy_1_Good_POST$avg_dxy),replace=TRUE)))
}
quantile(bstrap_medians,c(0.05,0.95))
```
## Figure 5C:  <a name="figure5c"></a>
```{r}
# Create plot theme
fst_theme <- theme(legend.position="none",panel.grid = element_blank(),
                     axis.text.y=element_text(face="bold", color="black"),
                     axis.text.x=element_text(face="bold", color="black"),
                     axis.title.y = element_text(face="bold", color="black"),
                     axis.ticks.x=element_blank(), 
                   strip.background = element_blank(),
                   panel.border = element_rect(color="black"),
                     axis.title.x=element_blank())
# Read in data
fst_treatment_2k <- read.table('fst.windows.2kb.treatment.txt', header=TRUE)
fst_treatment_2kb <- merge(fst_treatment_2k, windows, by.x=c("CHROM","BIN_START"), by.y=c("CHROM","BIN_START"))

# Order by chromosome and position
fst_treatment_2kb <- fst_treatment_2kb[order(fst_treatment_2kb$CHROM, fst_treatment_2kb$BIN_START),]
fst_treatment_2kb$ID <- seq.int(nrow(fst_treatment_2kb))
names(fst_treatment_2kb) <- c("CHROM", "BIN_START","BIN_END", "N_VARIANTS","WEIGHTED_FST","MEAN_FST","PRE_POST","BIN_STOP","SNP_COUNT","ID")

# Order placement of axis labels
axisdf = fst_treatment_2kb %>% group_by(CHROM) %>% dplyr::summarize(center=( max(ID) + min(ID) ) / 2 )

# Plot good clearers vs post-treatment 
GvsPO <- ggplot(subset(fst_treatment_2kb,SNP_COUNT>=20 & PRE_POST=="POST_GOOD"), aes(x=ID, y=as.numeric(WEIGHTED_FST))) +
  geom_point( aes(color=as.factor(CHROM)),size=0.0001) +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center, expand=c(0,0) ) +
  scale_y_continuous(expand=c(0,0), limits=c(0.0,0.2), breaks=c(0,0.1,0.2)) +
  theme_bw() + 
  labs(y=expression(bold("Median "*F[ST])))+
  fst_theme +
  theme(axis.text.x=element_blank())

# Plot good clearers vs pre-treatment 
GvsPR <- ggplot(subset(fst_treatment_2kb,SNP_COUNT>=20 & PRE_POST=="PRE_GOOD"), aes(x=ID, y=as.numeric(WEIGHTED_FST))) +
  geom_point( aes(color=as.factor(CHROM)),size=0.0001) +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center, expand=c(0,0) ) +
  scale_y_continuous(expand=c(0,0), limits=c(0.0,0.2), breaks=c(0,0.1,0.2)) +
  theme_bw() + 
  labs(y=expression(bold(F[ST])))+
  fst_theme+
  theme(axis.text.x=element_blank())

# Plot pre-treatment vs post-treatment 
PRvsPO <- ggplot(subset(fst_treatment_2kb,SNP_COUNT>=20 & PRE_POST=="PRE_POST"), aes(x=ID, y=as.numeric(WEIGHTED_FST))) +
  geom_point( aes(color=as.factor(CHROM)),size=0.0001)+
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center, expand=c(0,0) ) +
  scale_y_continuous(expand=c(0,0), limits=c(0.0,0.2), breaks=c(0,0.1,0.2)) +
  theme_bw() +
  labs(y=expression(bold("Weighted "*F[ST])))+
  fst_theme +
  theme(axis.text.x=element_blank()) 
```
## Figure 5D: Logistic regression genome-wide association test <a name="figure5d"></a>
```{r}
# Create theme
err_theme <- theme(legend.title = element_blank(),
                       legend.position = "none",
                       panel.grid = element_blank(),
                       strip.text.x = element_text(size=10, face="bold"),
                       axis.text.y = element_text(face="bold", color="black", size=10),
                       axis.text.x = element_text(face="bold", color="black", size=10),
                       axis.title.y = element_text(face="bold", color="black", size=10),
                       axis.title.x = element_text(face="bold", color="black", size=10),
                       strip.background = element_blank())
                       
# Read in data
assoc_BIN_ALL <- read.table("assoc_err_binary.txt", header=TRUE, sep='\t')

# Add a color variable for each chromosome
assoc_BIN_ALL$COLOR<- assoc_BIN_ALL$CHR
c$COLOR <- ifelse(assoc_BIN_ALL$CHR == 1,"A",
                              ifelse(assoc_BIN_ALL$CHR == 3,"A",
                                     ifelse(assoc_BIN_ALL$CHR == 5,"A",
                                            ifelse(assoc_BIN_ALL$CHR == 7,"A", no="B"))))

# Calculate a multiple correction cutoff value (Bonferroni correction)
bin_cutoff <- 0.05/(as.numeric(nrow(assoc_BIN_ALL)))

# Color variants above this threshold
assoc_BIN_ALL$COLOR[abs(assoc_BIN_ALL$FDR_BY)<=(bin_cutoff)] <- "RED"
# Remove a the least significant variants
assoc_BIN_ALL2 <- subset(assoc_BIN_ALL, UNADJ<0.9999)

# Order by chromosome and window, then assign a value to each window (for plotting x-axis)
assoc_BIN_ALL2 <- assoc_ERR_ALL[order(assoc_BIN_ALL2$CHR, assoc_BIN_ALL2$SNP),]
assoc_BIN_ALL2$ID <- seq.int(nrow(assoc_BIN_ALL2))

# Order chromosome labels
axisdf = assoc_BIN_ALL2 %>% group_by(CHR) %>% summarize(center=(max(as.numeric(ID)) + min(as.numeric(ID)) ) / 2 )

# Plot
ERR_BIN_UNADJ <- ggplot(assoc_BIN_ALL2, aes(x=as.numeric(ID), y=-log10(as.numeric(UNADJ)))) +
  geom_point(aes(color=as.factor(COLOR)),size=0.5) +
  scale_color_manual(values = rep(c("grey75", "grey40","red"))) +
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,10), breaks=c(0,2.5,5,7.5,10))+
  theme_bw() +
  labs(y=expression(bold(-log[10]*("p")))) +
  geom_hline(yintercept=-log10(bin_cutoff), linetype="dashed") +
  annotate("text",x=150000, y=8,parse=TRUE,size=3, label='bold("Good clearers vs Post-treatment")') +
  fst_theme +
  theme(axis.text.x = element_blank())
```
## Figure 5E: Linear regression genome-wide association test  <a name="figure5e"></a>
```{r}
# Read in data
assoc_ERR_ALL <- read.table("assoc_err_linear.txt", header=TRUE, sep='\t')

# Add a color variable for each chromosome
assoc_ERR_ALL$COLOR<- assoc_ERR_ALL$CHR
assoc_ERR_ALL$COLOR <- ifelse(assoc_ERR_ALL$CHR == 1,"A",
                              ifelse(assoc_ERR_ALL$CHR == 3,"A",
                                     ifelse(assoc_ERR_ALL$CHR == 5,"A",
                                            ifelse(assoc_ERR_ALL$CHR == 7,"A", no="B"))))

# Calculate a multiple correction cutoff value (Bonferroni correction)
err_cutoff <- 0.05/(as.numeric(nrow(assoc_ERR_ALL)))

# Color variants above this threshold
assoc_ERR_ALL$COLOR[abs(assoc_ERR_ALL$UNADJ)<=((err_cutoff))] <- "RED"

# Remove a the least significant variants
assoc_ERR_ALL2 <- subset(assoc_ERR_ALL, UNADJ<0.9999)

# Order by chromosome and window, then assign a value to each window (for plotting x-axis)
assoc_ERR_ALL2 <- assoc_ERR_ALL[order(assoc_ERR_ALL2$CHR, assoc_ERR_ALL2$SNP),]
assoc_ERR_ALL2$ID <- seq.int(nrow(assoc_ERR_ALL2))

# Order chromosome labels
axisdf = assoc_ERR_ALL2 %>% group_by(CHR) %>% summarize(center=(max(as.numeric(ID)) + min(as.numeric(ID)) ) / 2 )

# Plot
ERR_LIN_UNADJb <- ggplot(assoc_ERR_ALL2, aes(x=as.numeric(ID), y=-log10(as.numeric(UNADJ)))) +
  geom_point(aes(color=as.factor(COLOR)),size=0.5) +
  scale_color_manual(values = rep(c("grey75", "grey40","red"))) +
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,10), breaks=c(0,2.5,5,7.5,10))+
  theme_bw() +
  labs(y=expression(bold(-log[10]*("p")))) +
  annotate("text",x=152181224, y=8,parse=TRUE,size=3, label='bold("Host egg reduction rate")') +
  geom_hline(yintercept=-log10(err_cutoff), linetype="dashed") +
  fst_theme +
  theme(axis.text.x = element_blank())
```
## Merge figures  <a name="merge"></a>
```{r}
top_row <- plot_grid(pi_all,"",rel_widths = c(2,1), nrow=1, labels = c('A', 'B'))
bottom_row <- plot_grid(GvsPR,GvsPO,PRvsPO,label,ERR_BIN_UNADJ,ERR_LIN_UNADJb,label, 
          labels = c('C', '','','','D','E',''), align = 'v',
          rel_widths = c(1,1,1,1,1,1,1,1),
          rel_heights = c(1,1,1,0.2,1,1,0.2), nrow=7, ncol=1)

plot_grid(top_row,bottom_row, axis="r", nrow = 2, rel_heights = c(0.185,0.54), align = 'v')
```
