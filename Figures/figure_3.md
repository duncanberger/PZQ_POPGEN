# Figure 3: Signatures of recent selection

0. [Setup](#setup)
1. [Figure 3A](#figure3a)
2. [Figure 3B](#figure3b)
3. [Figure 3C](#figure3c)
4. [Figure 3D](#figure3d)
5. [Figure 3E](#figure3e)
6. [Merged figure](#figure3f)
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
?selection_colors <- rep(c("grey75", "grey40","#ef3b2c"))
?selection_colors2 <- rep(c("grey75", "grey40","#19abff"))
```
## Figure 3A: iHS <a name="figure3a"></a>
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
ihs_noko_summary$COLOR_ihs_noko[ihs_noko_summary$N_VARIANTS_XPEHH<25] <- "NA"
```
## Figure 3B: CLR <a name="figure3a"></a>
```{r}
# Read in data
clr_MAYUGE <- read.table("all.mayuge.swd", header=FALSE)




xpehh_allvsko <- read.table("ALL.MAYUGEvsTORORO.xpehh.xpehh.out.norm.fix", header=TRUE)
fst_district <- read.table("district.25kb.fst.txt", header=TRUE)

