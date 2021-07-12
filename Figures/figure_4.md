# Figure 4: Signatures of recent selection. 
0. [Setup](#setup)
1. [Figure 4A](#figure4a)
2. [Figure 4B](#figure4b)
3. [Figure 4C](#figure4c)
4. [Figure 4D](#figure4d)
5. [Figure 4E](#figure4e)
6. [Merge figures](#merge)
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
key <- read.table("supplementary_data_10.txt", header=TRUE, sep="\t", check.names = FALSE, comment.char = "")
```
## Setup <a name="figure2a"></a>
```{r}
# Load data
windows <- read.table("sm_2000.bed", header=FALSE)
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
## Figure 4A: Mayuge iHS <a name="figure4a"></a>
```{r}
# Read in data
ihs_noko <- read.table("ALL.MAYUGE.IHS.ihs.out.100bins.norm.txt")

#Create windows based on site locations
ihs_noko$START <- (RoundTo(ihs_noko$V3, multiple = 2000, FUN = floor)+1)
ihs_noko$STOP <- (RoundTo(ihs_noko$V3, multiple = 2000, FUN = ceiling))

# Calculate median iHS scores for each window and calculate number of variants per window
ihs_noko_median <- aggregate((abs(V8))~V1+START, ihs_noko, median)
ihs_noko_length <- aggregate((abs(V8))~V1+START, ihs_noko, FUN=length)

# Merge median iHS and variant counts per window 
ihs_noko_summary0 <- merge(ihs_noko_median, ihs_noko_length, by.x = c("V1","START"), by.y = c("V1","START"))
names(ihs_noko_summary0) <- c("CHROM", "BIN_START","MEDIAN_ihs_noko", "N_VARIANTS_ihs_noko")

# Merge window bed file (in case of missing windows)
ihs_noko_summary1 <- merge(ihs_noko_summary0, windows, by.x = c("CHROM","BIN_START"), by.y = c("CHROM","BIN_START"), all.x = TRUE, all.y = TRUE)

# Order windows and create x-axis coordinates
ihs_noko_summary1<- ihs_noko_summary1[order(ihs_noko_summary1$CHROM, ihs_noko_summary1$BIN_START),]
ihs_noko_summary1$ID <- seq.int(nrow(ihs_noko_summary1))

# Remove windows with less than 20 variants per window
ihs_noko_summary <- subset(ihs_noko_summary1, SNP_COUNT>=20)

# Identify highest 0.25% of median iHS scores and set a cutoff
tp_subset <- ihs_noko_summary[ihs_noko_summary$MEDIAN_ihs_noko >= quantile(ihs_noko_summary$MEDIAN_ihs_noko,.9975, na.rm = TRUE),]
ihsnoko_cutoff <- as.numeric(head(tp_subset[order(tp_subset$MEDIAN_ihs_noko),],1)[3])

# Group windows into 3 categories based on chromosome number and whether windows fall into highest 0.25% of values
ihs_noko_summary$COLOR_ihs_noko <- ihs_noko_summary$CHROM
ihs_noko_summary$COLOR_ihs_noko <- ifelse(ihs_noko_summary$CHROM == 1,"A",
                                          ifelse(ihs_noko_summary$CHROM == 3,"A",
                                                 ifelse(ihs_noko_summary$CHROM == 5,"A",
                                                        ifelse(ihs_noko_summary$CHROM == 7,"A", no="B"))))

# Subset and group regions of selection and combine into a final table
test <- subset(ihs_noko_summary, MEDIAN_ihs_noko>ihsnoko_cutoff) %>% 
  group_by(CHROM, grp = cumsum(c(1, diff(BIN_START) > 300000))) %>% 
  filter(n() > 4) %>% select(CHROM,BIN_START)
aggregate(test$BIN_START, by = list(test$grp, test$CHROM), max)
test$COLOR_ihs_noko <- "SELEC"
ihs_noko_summary_2 <- merge(ihs_noko_summary, test, by.x = c("CHROM","BIN_START"), by.y = c("CHROM","BIN_START"), all.x = TRUE, all.y = TRUE)
ihs_noko_summary_2$mycol <- coalesce(ihs_noko_summary_2$COLOR_ihs_noko.y, ihs_noko_summary_2$COLOR_ihs_noko.x)

# Set axis labels for windows
axisdf = ihs_noko_summary_2 %>% group_by(CHROM) %>% summarize(center=( max(ID) + min(ID) ) / 2 )

MAYUGE_IHS <- ggplot(subset(ihs_noko_summary_2,SNP_COUNT>=20), aes(x=ID, y=abs(MEDIAN_ihs_noko))) +
  geom_point( aes(color=as.factor(mycol)),size=0.0001) +
  ylab("Median | iHS |") +
  scale_color_manual(values = selection_colors) +
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,6), breaks=c(0.0,3,6.0)) +
  theme_bw() +
  selection_theme
```
## Figure 4B: Tororo iHS <a name="figure4b"></a>
```{r}
# Import data
ihs_ko <- read.table("ALL.TORORO.IHS.ihs.out.100bins.norm.txt")

#Create windows based on site locations
ihs_ko$START <- (RoundTo(ihs_ko$V3, multiple = 2000, FUN = floor)+1)
ihs_ko$STOP <- (RoundTo(ihs_ko$V3, multiple = 2000, FUN = ceiling))

# Calculate median iHS scores for each window and calculate number of variants per window
ihs_ko_median <- aggregate((abs(V8))~V1+START, ihs_ko, median)
ihs_ko_length <- aggregate((abs(V8))~V1+START, ihs_ko, FUN=length)

# Merge median iHS and variant counts per window 
ihs_ko_summary0 <- merge(ihs_ko_median, ihs_ko_length, by.x = c("V1","START"), by.y = c("V1","START"))
names(ihs_ko_summary0) <- c("CHROM", "BIN_START","MEDIAN_ihs_ko", "N_VARIANTS_ihs_ko")

# Merge window bed file (in case of missing windows)
ihs_ko_summary7 <- merge(ihs_ko_summary0, windows, by.x = c("CHROM","BIN_START"), by.y = c("CHROM","BIN_START"), all.x = TRUE, all.y = TRUE)

# Order windows and create x-axis coordinates
ihs_ko_summary7<- ihs_ko_summary7[order(ihs_ko_summary7$CHROM, ihs_ko_summary7$BIN_START),]
ihs_ko_summary7$ID <- seq.int(nrow(ihs_ko_summary7))

# Remove windows with less than 20 variants per window
ihs_ko_summary <- subset(ihs_ko_summary7, SNP_COUNT>=20)

# Identify highest 0.25% of median iHS scores
tp_subset <- ihs_ko_summary[ihs_ko_summary$MEDIAN_ihs_ko >= quantile(ihs_ko_summary$MEDIAN_ihs_ko,.9975, na.rm=TRUE),]
ihsko_cutoff <- as.numeric(head(tp_subset[order(tp_subset$MEDIAN_ihs_ko),],1)[3])

# Group windows into 3 categories based on chromosome number and whether windows fall into highest 0.25% of values
ihs_ko_summary$COLOR_ihs_ko <- ihs_ko_summary$CHROM
ihs_ko_summary$COLOR_ihs_ko <- ifelse(ihs_ko_summary$CHROM == 1,"A",
                                          ifelse(ihs_ko_summary$CHROM == 3,"A",
                                                 ifelse(ihs_ko_summary$CHROM == 5,"A",
                                                        ifelse(ihs_ko_summary$CHROM == 7,"A", no="B"))))

TORORO_IHS<- ggplot(subset(ihs_ko_summary,SNP_COUNT>=20), aes(x=ID, y=abs(MEDIAN_ihs_ko))) +
  geom_point( aes(color=as.factor(CHROM)),size=0.0001) +
  ylab("Median | iHS |") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,6), breaks=c(0.0,3.0,6.0)) +
  theme_bw() +
  selection_theme
```
## Figure 4C: XP-EHH <a name="figure3c"></a>
```{r}
# Read in data
xpehh_allvsko <- read.table("ALL.MAYUGEvsTORORO.xpehh.xpehh.out.norm.txt", header=TRUE)

# Create windows based on site locations
xpehh_allvsko$START <- (RoundTo(as.numeric(xpehh_allvsko$pos), multiple = 2000, FUN = floor)+1)
xpehh_allvsko$STOP <- (RoundTo(as.numeric(xpehh_allvsko$pos), multiple = 2000, FUN = ceiling))
xpehh_allvsko_median <- aggregate(as.numeric(normxpehh)~CHROM+START, xpehh_allvsko, median)
xpehh_allvsko_length <- aggregate(as.numeric(normxpehh)~CHROM+START, xpehh_allvsko, FUN=length)

# Merge median XP-EHH and variant counts per window 
xpehh_summary0 <- merge(xpehh_allvsko_median, xpehh_allvsko_length, by.x = c("CHROM","START"), by.y = c("CHROM","START"))
names(xpehh_summary0) <- c("CHROM", "BIN_START","MEDIAN_XPEHH", "N_VARIANTS_XPEHH")
xpehh_summary1 <- merge(xpehh_summary0, windows, by.x = c("CHROM","BIN_START"), by.y = c("CHROM","BIN_START"), all.x = TRUE, all.y = TRUE)

# Order windows and create x-axis coordinates
xpehh_summary1<- xpehh_summary1[order(xpehh_summary1$CHROM, xpehh_summary1$BIN_START),]
xpehh_summary1$ID <- seq.int(nrow(xpehh_summary1))

# Remove windows with less than 20 variants per window
xpehh_summary <- subset(xpehh_summary1, SNP_COUNT>=20)

# Identify highest 0.25% of median XP-EHH scores
tp_subset <- xpehh_summary[xpehh_summary$MEDIAN_XPEHH >= quantile(xpehh_summary$MEDIAN_XPEHH,.9975, na.rm=TRUE),]
xpehh_cutoff <- as.numeric(head(tp_subset[order(tp_subset$MEDIAN_XPEHH),],1)[3])

# Group windows into 3 categories based on chromosome number and whether windows fall into highest 0.25% of values
xpehh_summary$COLOR_XPEHH <- xpehh_summary$CHROM
xpehh_summary$COLOR_XPEHH <- ifelse(xpehh_summary$CHROM == 1,"A",
                                            ifelse(xpehh_summary$CHROM == 3,"A",
                                                   ifelse(xpehh_summary$CHROM == 5,"A",
                                                          ifelse(xpehh_summary$CHROM == 7,"A", no="B"))))

# Subset and group regions of selection and combine into a final table
test_XPEHH <- subset(xpehh_summary, xpehh_summary$MEDIAN_XPEHH>xpehh_cutoff) %>% 
  group_by(CHROM, grp = cumsum(c(1, diff(BIN_START) > 150000))) %>% 
  filter(n() > 4) %>% select(CHROM,BIN_START)
aggregate(test_XPEHH$BIN_START, by = list(test_XPEHH$grp, test_XPEHH$CHROM), min)
test_XPEHH$COLOR_XPEHH <- "SELEC"
xpehh_summary_2 <- merge(xpehh_summary, test_XPEHH, by.x = c("CHROM","BIN_START"), by.y = c("CHROM","BIN_START"), all.x = TRUE, all.y = TRUE)
xpehh_summary_2$mycol <- coalesce(xpehh_summary_2$COLOR_XPEHH.y, xpehh_summary_2$COLOR_XPEHH.x)

# Set axis labels for windows
axisdf = xpehh_summary_2 %>% group_by(CHROM) %>% summarize(center=( max(ID) + min(ID) ) / 2 )

# Plot figure
XPEHH <- ggplot(subset(xpehh_summary_2,SNP_COUNT>20), aes(x=ID, y=(MEDIAN_XPEHH))) +
  geom_point( aes(color=as.factor(mycol)),size=0.0001) +
  ylab("Median XP-EHH") +
  scale_color_manual(values = selection_colors) +
  scale_x_continuous(label = axisdf$CHROM, breaks= axisdf$center, expand=c(0,0) ) +
  scale_y_continuous(expand=c(0,0), limits=c(-10.0,10.0), breaks=c(-10.0,0.0,10.0),labels=scaleFUN)+
  theme_bw() + 
  selection_theme 
```
## Figure 4D: F<sub>ST</sub> <a name="figure4d"></a>
```{r}
# Read in data
fst_district <- read.table("MAYUGE_TORORO_2000.windowed.weir.txt", header=TRUE)

# Merge data into windows 
fst_district_summary <- merge(fst_district, windows, by.x = c("CHROM","BIN_START"), by.y = c("CHROM","BIN_START"), all.x = TRUE, all.y = TRUE)
fst_district_summary<- fst_district_summary[order(fst_district_summary$CHROM, fst_district_summary$BIN_START),]

# Order windows by x-axis coordinate
fst_district_summary$ID <- seq.int(nrow(fst_district_summary))

# Remove windows with less than 20 variants per window
fst_district_2 <- subset(fst_district_summary, SNP_COUNT>=20)

# Identify highest 0.25% of median FST scores
fst_district_median <- fst_district_2[fst_district_2$WEIGHTED_FST  >= quantile(fst_district_2$WEIGHTED_FST,.9975),]
FST_cutoff <- as.numeric(head(fst_district_median[order(fst_district_median$WEIGHTED_FST ),],1)[5])

# Group windows into 3 categories based on chromosome number and whether windows fall into highest 0.25% of values
fst_district_2$COLOR_FST <- fst_district_2$CHROM
fst_district_2$COLOR_FST <- ifelse(fst_district_2$CHROM == 1,"A",
                                 ifelse(fst_district_2$CHROM == 3,"A",
                                        ifelse(fst_district_2$CHROM == 5,"A",
                                               ifelse(fst_district_2$CHROM == 7,"A", no="B"))))

test_FST <- subset(fst_district_2, fst_district_2$WEIGHTED_FST>FST_cutoff) %>% 
  group_by(CHROM, grp = cumsum(c(1, diff(BIN_START) > 150000))) %>% 
  filter(n() > 4) %>% select(CHROM,BIN_START)
aggregate(test_FST$BIN_START, by = list(test_FST$grp, test_FST$CHROM), min)
test_FST$COLOR_FST <- "SELEC"
fst_district_3 <- merge(fst_district_2, test_FST, by.x = c("CHROM","BIN_START"), by.y = c("CHROM","BIN_START"), all.x = TRUE, all.y = TRUE)
fst_district_3$mycol <- coalesce(fst_district_3$COLOR_FST.y, fst_district_3$COLOR_FST.x)

# Set x-axis labels
axisdf = fst_district_3 %>% group_by(CHROM) %>% summarize(center=( max(ID) + min(ID) ) / 2 )

FST <- ggplot(subset(fst_district_3,SNP_COUNT>=20), aes(x=ID, y=WEIGHTED_FST)) +
  geom_point( aes(color=as.factor(mycol)),size=0.0003) +
  labs(y=expression(bold(Median*" "*F[ST])))+
  scale_color_manual(values = selection_colors) +
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center, expand=c(0,0) ) +
  scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks=c(0.0,0.2,0.4,0.6,0.8,1))+
  theme_bw() +
  selection_theme
```
## Figure 4E: Ratio of nucleotide diversity  <a name="figure4e"></a>
```{r}
# Import data
PI_TORORO <- read.table("TORORO_PI.windowed.pi", header=TRUE, na.strings = "nan")
PI_MAYUGE <- read.table("MAYUGE_PI.windowed.pi", header=TRUE, na.strings = "nan")

# Merge datasets for eacch district then merge with set 2kb windows (to check for missingness)
PI_merged0 <- merge(PI_MAYUGE, PI_TORORO, by.x = c("CHROM","BIN_START"), by.y = c("CHROM","BIN_START"), all.x = TRUE, all.y = TRUE)
PI_merged <- merge(PI_merged0, windows,by.x = c("CHROM","BIN_START"), by.y = c("CHROM","BIN_START"), all.x = TRUE, all.y = TRUE)

# Order by chromosome and coordinate
PI_merged <- PI_merged[order(PI_merged$CHROM, PI_merged$BIN_START),]
PI_merged$ID <- seq.int(nrow(PI_merged))

# Color chromosomes
PI_merged$COLOR_RPI <- PI_merged$CHROM
PI_merged$COLOR_RPI <- ifelse(PI_merged$CHROM == 1,"A",
                                      ifelse(PI_merged$CHROM == 3,"A",
                                             ifelse(PI_merged$CHROM == 5,"A",
                                                    ifelse(PI_merged$CHROM == 7,"A", no="B"))))

# Set x-axis order
axisdf = PI_merged %>% dplyr::group_by(CHROM) %>% dplyr::summarize(center=( max(ID) + min(ID) ) / 2 )

# Plot
PI_TM <- ggplot(subset(PI_merged,N_VARIANTS.x>=20), aes(x=ID, y=log10(PI.x/PI.y))) +
  geom_point( aes(color=as.factor(CHROM)),size=0.0001) +
  labs(y=expression(bold(log[10]("\U03C0"[Mayuge]*" "*"/"*""*" \U03C0"[Tororo])))) +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(-3,3), breaks=c(-3,-2,-1,0,1,2,3)) +
  theme_bw() +
  selection_theme
```
## Merge figures <a name="merge"></a>
```{r}
plot_grid(MAYUGE_IHS,TORORO_IHS,XPEHH,FST,PI_TM,nrow=5, align = "v", rel_heights = c(1,1,1,1,1),labels=c("A","B","C","D","E"),label_y = 1.0)
```
