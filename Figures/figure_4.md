# Figure 4: Signatures of recent selection

0. [Setup](#setup)
1. [Figure 4A](#figure4a)
2. [Figure 4B](#figure4b)
3. [Figure 4D](#figure4d)
4. [Figure 4E](#figure4e)
6. [Plot individual figures](#pindv)
7. [Merge figures](#merge)
## Project setup <a name="setup"></a>
```{r}
# Load required packages
library("DescTools")
library("dplyr")
library("cowplot")
library("ggplot2")
library("reshape2")
library("lattice")

# Load metadata
key <- read.table("supplementary_table_2.txt", header=TRUE, sep="\t", check.names = FALSE, comment.char = "")
```
## Setup <a name="figure2a"></a>
```{r}
# Load data
windows <- read.table("Sm_v7_nohap.fa.bed", header=FALSE)
names(windows) <- c("CHROM","BIN_START", "BIN_STOP")

# Set ggplot themes and palettes
selection_theme <- theme(
  legend.position="none",
  panel.grid = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y=element_text(face="bold", color="black"),
  axis.title.y = element_text(face="bold", color="black"),
  axis.title.x = element_blank(),
  axis.ticks.x=element_blank())

scaleFUN <- function(x) sprintf("%.1f", x)
selection_colors <- rep(c("grey75", "grey40","#ef3b2c"))
selection_colors2 <- rep(c("grey75", "grey40","#19abff"))
```
## Figure 4A: iHS <a name="figure3a"></a>
```{r}
# Read in data
ihs_noko <- read.table("mayuge.ihs.out.100bins.norm.all")

#Create windows based on site locations
ihs_noko$START <- (RoundTo(ihs_noko$V3, multiple = 25000, FUN = floor)+1)
ihs_noko$STOP <- (RoundTo(ihs_noko$V3, multiple = 25000, FUN = ceiling))

# Calculate median iHS scores for each window and calculate number of variants per window
ihs_noko_median <- aggregate((abs(V8))~V1+START, ihs_noko, median)
ihs_noko_length <- aggregate((abs(V8))~V1+START, ihs_noko, FUN=length)

# Merge median iHS and variant counts per window 
ihs_noko_summary <- merge(ihs_noko_median, ihs_noko_length, by.x = c("V1","START"), by.y = c("V1","START"))

# Rename table headers
names(ihs_noko_summary) <- c("CHROM", "BIN_START","MEDIAN_ihs_noko", "N_VARIANTS_ihs_noko")

# Identify highest 0.5% of median iHS scores
tp_subset <- ihs_noko_summary[ihs_noko_summary$MEDIAN_ihs_noko >= quantile(ihs_noko_summary$MEDIAN_ihs_noko,.9950),]
ihsnoko_cutoff <- as.numeric(head(tp_subset[order(tp_subset$MEDIAN_ihs_noko),],1)[3])

# Group windows into 3 categories based on chromosome number and whether windows fall into highest 0.5% of values
ihs_noko_summary$COLOR_ihs_noko <- ihs_noko_summary$CHROM
ihs_noko_summary$COLOR_ihs_noko <- ifelse(ihs_noko_summary$CHROM == 1,"A",
                                          ifelse(ihs_noko_summary$CHROM == 3,"A",
                                                 ifelse(ihs_noko_summary$CHROM == 5,"A",
                                                        ifelse(ihs_noko_summary$CHROM == 7,"A", no="B"))))
ihs_noko_summary$COLOR_ihs_noko[ihs_noko_summary$MEDIAN_ihs_noko>(ihsnoko_cutoff)] <- "TP"
```
## Figure 4B: CLR <a name="figure3b"></a>
```{r}
# Read in data
clr_MAYUGE <- read.table("all.mayuge.swd", header=FALSE)

#Create windows based on site locations
clr_MAYUGE$START <- (RoundTo(clr_MAYUGE$V1, multiple = 25000, FUN = floor)+1)
clr_MAYUGE$STOP <- (RoundTo(clr_MAYUGE$V1, multiple = 25000, FUN = ceiling))

# Merge median iHS and variant counts per window 
clr_MAYUGE_median <- aggregate((V2)~V4+START, clr_MAYUGE, median)
clr_MAYUGE_length <- aggregate((V2)~V4+START, clr_MAYUGE, FUN=length)

# Merge median CLR and variant counts per window 
clr_MAYUGE_summary <- merge(clr_MAYUGE_median, clr_MAYUGE_length, by.x = c("V4","START"), by.y = c("V4","START"))

# Rename table headers
names(clr_MAYUGE_summary) <- c("CHROM", "BIN_START","MEDIAN_clr_MAYUGE", "N_WIN_clr_MAYUGE")

# Identify highest 0.5% of median CLR scores
tp_subset <- clr_MAYUGE_summary[clr_MAYUGE_summary$MEDIAN_clr_MAYUGE >= quantile(clr_MAYUGE_summary$MEDIAN_clr_MAYUGE,.9950),]
clr_MAYUGE_cutoff <- as.numeric(head(tp_subset[order(tp_subset$MEDIAN_clr_MAYUGE),],1)[3])

# Group windows into 3 categories based on chromosome number and whether windows fall into highest 0.5% of values
clr_MAYUGE_summary$COLOR_clr_MAYUGE <- clr_MAYUGE_summary$CHROM
clr_MAYUGE_summary$COLOR_clr_MAYUGE<- ifelse(clr_MAYUGE_summary$CHROM == 1,"A",
                                          ifelse(clr_MAYUGE_summary$CHROM == 3,"A",
                                                 ifelse(clr_MAYUGE_summary$CHROM == 5,"A",
                                                        ifelse(clr_MAYUGE_summary$CHROM == 7,"A", no="B"))))
clr_MAYUGE_summary$COLOR_clr_MAYUGE[clr_MAYUGE_summary$MEDIAN_clr_MAYUGE>(clr_MAYUGE_cutoff)] <- "TP"
```
## Figure 4D: XP-EHH <a name="figure3d"></a>
```{r}
# Read in data
xpehh_allvsko <- read.table("xpehh.out.norm.all", header=TRUE)

# Create windows based on site locations
xpehh_allvsko_median <- aggregate(as.numeric(normxpehh)~CHROM+START, xpehh_allvsko, median)
xpehh_allvsko_length <- aggregate(as.numeric(normxpehh)~CHROM+START, xpehh_allvsko, FUN=length)

# Merge median iHS and variant counts per window 
xpehh_allvsko$START <- (RoundTo(as.numeric(xpehh_allvsko$pos), multiple = 25000, FUN = floor)+1)
xpehh_allvsko$STOP <- (RoundTo(as.numeric(xpehh_allvsko$pos), multiple = 25000, FUN = ceiling))

# Merge median XP-EHH and variant counts per window 
xpehh_summary <- merge(xpehh_allvsko_median, xpehh_allvsko_length, by.x = c("CHROM","START"), by.y = c("CHROM","START"))

# Rename table headers
names(xpehh_summary) <- c("CHROM", "BIN_START","MEDIAN_XPEHH", "N_VARIANTS_XPEHH")

# Identify highest 0.5% of median XP-EHH scores
tp_subset <- xpehh_summary[xpehh_summary$MEDIAN_XPEHH >= quantile(xpehh_summary$MEDIAN_XPEHH,.9950),]
xpehh_cutoff <- as.numeric(head(tp_subset[order(tp_subset$MEDIAN_XPEHH),],1)[3])

# Group windows into 3 categories based on chromosome number and whether windows fall into highest 0.5% of values
xpehh_summary$COLOR_XPEHH <- xpehh_summary$CHROM
xpehh_summary$COLOR_XPEHH <- ifelse(xpehh_summary$CHROM == 1,"A",
                                    ifelse(xpehh_summary$CHROM == 3,"A",
                                           ifelse(xpehh_summary$CHROM == 5,"A",
                                                  ifelse(xpehh_summary$CHROM == 7,"A", no="B"))))
xpehh_summary$COLOR_XPEHH[xpehh_summary$MEDIAN_XPEHH>(xpehh_cutoff)] <- "TP"
```
## Figure 4E: F<sub>ST</sub> <a name="figure3e"></a>
```{r}
# Read in data
fst_district <- read.table("district.25kb.fst.txt", header=TRUE)

# Identify highest 0.5% of median CLR scores
fst_district_median <- fst_district[fst_district$avg_wc_fst  >= quantile(fst_district$avg_wc_fst,.9950),]
FST_cutoff <- as.numeric(head(fst_district_median[order(fst_district_median$avg_wc_fst ),],1)[6])

# Group windows into 3 categories based on chromosome number and whether windows fall into highest 0.5% of values
fst_district$COLOR_FST <- fst_district$chromosome
fst_district$COLOR_FST <- ifelse(fst_district$chromosome == 1,"A",
                                 ifelse(fst_district$chromosome == 3,"A",
                                        ifelse(fst_district$chromosome == 5,"A",
                                               ifelse(fst_district$chromosome == 7,"A", no="B"))))
fst_district$COLOR_FST[fst_district$avg_wc_fst>FST_cutoff] <- "TP"
```
## Merge and combine datasets <a name="figure3e"></a>
```{r}
# Merge datasets together, one by one
summary_A <- merge(windows,ihs_noko_summary, by.x = c("CHROM","BIN_START"),by.y = c("CHROM","BIN_START"), all.x = TRUE, all.y=TRUE)
summary_B <- merge(summary_A, xpehh_summary,by.x = c("CHROM","BIN_START"),by.y = c("CHROM","BIN_START"), all.x = TRUE, all.y=TRUE)
summary_C <- merge(summary_B, fst_district ,by.x = c("CHROM","BIN_START"),by.y = c("chromosome","window_pos_1"), all.x = TRUE, all.y=TRUE)
summary_D <- merge(summary_C, clr_MAYUGE_summary ,by.x = c("CHROM","BIN_START"),by.y = c("CHROM","BIN_START"), all.x = TRUE, all.y=TRUE)
summary_E <- merge(summary_D, snpden_length,by.x = c("CHROM","BIN_START","BIN_STOP"),by.y = c("CHROM","BIN_START","BIN_STOP"), all.x = TRUE, all.y=TRUE)

# Order chromosomes by chromosomes and position
summary_E<- summary_f1[order(summary_E$CHROM, summary_E$BIN_START),]

#Create a consistent position value for each window
summary_E$ID <- seq.int(nrow(summary_E))

# Identify continuous regions of elevated CLR and iHS scores and summarise these in a dataframe
summary_G <- subset(summary_E, COLOR_clr_MAYUGE=="TP" | COLOR_ihs_noko=="TP") %>% 
  filter(DENSITY>100) %>%
  group_by(CHROM, grp = cumsum(c(1, diff(ID) > 12))) %>% 
  filter(n() > 2)
summary_H <- summary_G %>% group_by(grp, CHROM) %>% summarise(ID_START = min(ID), ID_STOP=max(ID), BIN_START=min(BIN_START), BIN_STOP=max(BIN_STOP))

# Identify continuous regions of elevated XP-EHH and FST scores and summarise these in a dataframe
summary_G2 <- subset(summary_f1, COLOR_FST=="TP" | COLOR_XPEHH=="TP") %>% 
  filter(DENSITY>100) %>%
  group_by(CHROM, grp = cumsum(c(1, diff(ID) > 12))) %>% 
  filter(n() > 2)
summary_H2 <- summary_G2 %>% group_by(grp, CHROM) %>% summarise(ID_START = min(ID), ID_STOP=max(ID), BIN_START=min(BIN_START), BIN_STOP=max(BIN_STOP))
```
## Plot merged datasets in individual figures <a name="pindv"></a>
```{r}
# Order axis labels on chromosomes for each plot
axisdf = summary_f1 %>% group_by(CHROM) %>% summarize(center=( max(ID) + min(ID) ) / 2 )

# Plot iHS - figure 4A
IHS_NOKO <- ggplot(subset(summary_f1,N_VARIANTS_ihs_noko >25 & DENSITY>100), aes(x=ID, y=abs(MEDIAN_ihs_noko))) +
  geom_point( aes(color=as.factor(COLOR_ihs_noko)),size=0.3) +
  annotate("text",x=13500, y=4.6, label="Mayuge",fontface = "bold") +
  ylab("Median | iHS |") +
  scale_color_manual(values = selection_colors) +
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,5), breaks=c(0.0,2.5,5.0)) +
  theme_bw() +
  selection_theme

# Plot CLR - figure 4B
CLR_MAYUGE <- ggplot(subset(summary_f1, DENSITY>100), aes(x=ID, y=log10(MEDIAN_clr_MAYUGE))) +
  geom_point( aes(color=as.factor(COLOR_clr_MAYUGE)),size=0.3) +
  annotate("text",x=13500, y=4, label="Mayuge",fontface = "bold") +
  labs(y=expression(bold(-log[10]*("CLR")))) +
  scale_color_manual(values = selection_colors) +
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(-2,6), breaks=c(-2,0,2,4,6)) +
  theme_bw() +
  selection_theme

# Plot selected regions - figure 4B
HIGHLIGHT <- ggplot(summary_H) + 
  geom_rect(aes(xmin=(ID_START-0.5), ymin=0, xmax=(as.numeric(ID_STOP+0.5)), ymax=1), fill="#ef3b2c")+
  geom_rect(aes(xmin=3602, ymin=0, xmax=3602, ymax=1), fill=NA, color="black", size=0.2) +
  geom_rect(aes(xmin=5552, ymin=0, xmax=5552, ymax=1), fill=NA, color="black", size=0.2) +
  geom_rect(aes(xmin=7595, ymin=0, xmax=7595, ymax=1), fill=NA, color="black", size=0.2) +
  geom_rect(aes(xmin=9512, ymin=0, xmax=9512, ymax=1), fill=NA, color="black", size=0.2) +
  geom_rect(aes(xmin=10536, ymin=0, xmax=10536, ymax=1), fill=NA, color="black", size=0.2) +
  geom_rect(aes(xmin=11545, ymin=0, xmax=11545, ymax=1), fill=NA, color="black", size=0.2) +
  geom_rect(aes(xmin=12331, ymin=0, xmax=12331, ymax=1), fill=NA, color="black", size=0.2) +
  geom_rect(aes(xmin=0, ymin=0, xmax=15911, ymax=1), fill=NA, color="black", size=0.8) +
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,1)) + ylab(expression(bold(paste("Selected \n regions")))) +
  theme(axis.title.x = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.y = element_text(face="bold", size = 10,angle = 0, hjust = 1, vjust=1),
        axis.text.x = element_text(face="bold", size = 10),
        axis.ticks =element_blank(),
        axis.line = element_blank())

# Plot XP-EHH - figure 4D
XPEHH <- ggplot(subset(summary_f1,N_VARIANTS_XPEHH >25 & DENSITY>100), aes(x=ID, y=(MEDIAN_XPEHH))) +
  geom_point( aes(color=as.factor(COLOR_XPEHH)),size=0.3) +
  annotate("text",x=13500, y=8.0, label=paste0("Mayuge vs Tororo"), fontface = "bold") +
  ylab("XP-EHH") +
  scale_color_manual(values = selection_colors2) +
  scale_x_continuous(label = axisdf$CHROM, breaks= axisdf$center, expand=c(0,0) ) +
  scale_y_continuous(expand=c(0,0), limits=c(-10.0,10.0), breaks=c(-10.0,0.0,10.0),labels=scaleFUN)+
  theme_bw() + 
  selection_theme 

# Plot FST - figure 4E
FST <- ggplot(subset(summary_f1,no_snps>300 & DENSITY>100), aes(x=ID, y=(avg_wc_fst))) +
  geom_point( aes(color=as.factor(COLOR_FST)),size=0.3) +
  #  annotate("text",x=13500, y=0.65, label=paste0("Mayuge vs Tororo"), fontface = "bold") +
  annotate("text",x=13500, y=0.7, label=paste0("Mayuge vs Tororo"), fontface = "bold") +
  labs(y=expression(bold(F[ST])))+
  scale_color_manual(values = selection_colors2) +
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center, expand=c(0,0) ) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.8), breaks=c(0.0,0.2,0.4,0.6,0.8))+
  theme_bw() +
  selection_theme

# Plot selected regions - figure 4F
HIGHLIGHT2 <- ggplot(summary_H2) + 
  geom_rect(aes(xmin=(ID_START-0.5), ymin=0, xmax=(as.numeric(ID_STOP+0.5)), ymax=1), fill="#19abff")+
  geom_rect(aes(xmin=3602, ymin=0, xmax=3602, ymax=1), fill=NA, color="black", size=0.2) +
  geom_rect(aes(xmin=5552, ymin=0, xmax=5552, ymax=1), fill=NA, color="black", size=0.2) +
  geom_rect(aes(xmin=7595, ymin=0, xmax=7595, ymax=1), fill=NA, color="black", size=0.2) +
  geom_rect(aes(xmin=9512, ymin=0, xmax=9512, ymax=1), fill=NA, color="black", size=0.2) +
  geom_rect(aes(xmin=10536, ymin=0, xmax=10536, ymax=1), fill=NA, color="black", size=0.2) +
  geom_rect(aes(xmin=11545, ymin=0, xmax=11545, ymax=1), fill=NA, color="black", size=0.2) +
  geom_rect(aes(xmin=12331, ymin=0, xmax=12331, ymax=1), fill=NA, color="black", size=0.2) +
  geom_rect(aes(xmin=0, ymin=0, xmax=15911, ymax=1), fill=NA, color="black", size=0.8) +
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,1)) + ylab(expression(bold(paste("Selected \n regions")))) +
  theme(axis.title.x = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.y = element_text(face="bold", size = 10,angle = 0, hjust = 1, vjust=1),
        axis.text.x = element_text(face="bold", size = 10),
        axis.ticks =element_blank(),
        axis.line = element_blank())
```
## Merge figures <a name="merge"></a>
```{r}
plot_grid(IHS_NOKO,CLR_MAYUGE,HIGHLIGHT,XPEHH,FST,HIGHLIGHT2, align = "v",nrow = 6,rel_heights = c(1,1,0.325,1,1,0.325),labels=c('A','B','C','D','E','F'))
```
