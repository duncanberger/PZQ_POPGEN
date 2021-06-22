# Supplementary Figure 2: Variant quality control
```{r}
# Import data
summary_vQC <- read.table("cohort.genotyped.txt", header=TRUE)
site_QC <- read.table("hard_filtered_filtindv.lmiss.txt", header=TRUE)
sample_QC <- read.table("hard_filtered.imiss.txt", header=TRUE)

# Set ggplot theme
qc_theme <-   theme( legend.position="none",panel.grid = element_blank(), 
                     axis.title.y = element_text(face="bold", color="black", size=10),
                     axis.text.x=element_text(face="bold", color="black", size=8),
                     axis.text.y=element_text(face="bold", color="black", size=8),
                     axis.title.x = element_text(face="bold", color="black", size=10),  
                     panel.border = element_rect(color="black",fill = NA),
                     panel.background = element_blank()) 
```
## Variant QC plots
```{r}
# Designate target chromosomes (ones that are actually going to be analysed)
target_chr <- c("SM_V7_1","SM_V7_2","SM_V7_3","SM_V7_4","SM_V7_5","SM_V7_6","SM_V7_7","SM_V7_ZW")
summary_vQC_subset <- summary_vQC[summary_vQC$CHROM %in% target_chr,]

QD <- ggplot(data=summary_vQC) + 
  geom_histogram(aes(x=QD, fill=TYPE), bins=250) +
  geom_vline(xintercept = 2, linetype="dashed", color="grey50") +
  scale_y_continuous(expand=c(0,0),limits=c(0,2.6e+05), labels = function(x) format(x, scientific = TRUE)) +
  scale_x_continuous(expand=c(0,0), limits=c(0,40)) +
  scale_fill_manual(values=c("#ca0020","#ca0020","lightblue")) +
  xlab("Quality by Depth (QD)") + 
  ylab("") +
  qc_theme

FS <- ggplot(data=summary_vQC) + 
  geom_histogram(aes(x=FS, fill=TYPE), bins=100) +
  geom_vline(xintercept = 60, linetype="dashed", color="grey50") +
  scale_y_continuous(expand=c(0,0),limits=c(0,4e+06), labels = function(x) format(x, scientific = TRUE)) +
  scale_x_continuous(expand=c(0,0), limits=c(0,1000)) +
  xlab("FisherStrand (FS)") + 
  scale_fill_manual(values=c("#ca0020","#ca0020","lightblue")) +
  ylab("") +
  qc_theme

MQ <- ggplot(data=summary_vQC) + 
  geom_histogram(aes(x=MQ, fill=TYPE), bins=250) +
  geom_vline(xintercept = 40, linetype="dashed", color="grey50") +
  scale_y_continuous(expand=c(0,0),limits=c(0,2e+07), labels = function(x) format(x, scientific = TRUE)) +
  scale_x_continuous(expand=c(0,0), limits=c(0,75)) +
  scale_fill_manual(values=c("#ca0020","#ca0020","lightblue")) +
  xlab("RMSMappingQuality (MQ)") + 
  ylab("") +
  qc_theme

MQRankSum <- ggplot(data=summary_vQC) + 
  geom_histogram(aes(x=MQRankSum, fill=TYPE), bins=100) +
  geom_vline(xintercept = -12.5, linetype="dashed", color="grey50") +
  scale_y_continuous(expand=c(0,0),limits=c(0,2.0e+07), labels = function(x) format(x, scientific = TRUE)) +
  scale_x_continuous(expand=c(0,0), limits=c(-30,30)) +
  xlab("MappingQualityRankSumTest (MQRankSum)") + 
  scale_fill_manual(values=c("#ca0020","#ca0020","lightblue")) +
  ylab("") + 
  qc_theme

ReadPosRankSum <- ggplot(data=summary_vQC) +
  geom_histogram(aes(x=ReadPosRankSum, fill=TYPE), bins=100) +
  geom_vline(xintercept = -8, linetype="dashed", color="grey50") +
  scale_y_continuous(expand=c(0,0), limits=c(0,2e+07),labels = function(x) format(x, scientific = TRUE)) +
  scale_x_continuous(expand=c(0,0), limits=c(-50,50)) +
  xlab("ReadPosRankSumTest (ReadPosRankSum)") + 
  scale_fill_manual(values=c("#ca0020","#ca0020","lightblue")) +
  ylab("") +
  qc_theme

SOR <- ggplot(data=summary_vQC) +
  geom_histogram(aes(x=SOR, fill=TYPE), bins=100) +
  geom_vline(xintercept = 3, linetype="dashed", color="grey50") +
  scale_y_continuous(expand=c(0,0),limits=c(0,6e+06), labels = function(x) format(x, scientific = TRUE)) +
  scale_x_continuous(expand=c(0,0), limits=c(0,20)) +
  scale_fill_manual(values=c("#ca0020","#ca0020","lightblue")) +
  xlab("StrandOddsRatio (SOR)") + 
  ylab("") +  
  qc_theme

sample <- ggplot(data=sample_QC) +
  geom_histogram(aes(x=F_MISS), bins=100, fill="#af8dc3") +
  geom_vline(xintercept = 0.55, linetype="dashed", color="grey50") +
  scale_y_continuous(expand=c(0,0), limits=c(0,50),labels = function(x) format(x, scientific = TRUE)) +
  scale_x_continuous(expand=c(0,0), limits=c(0,1)) +
  xlab("Per-sample missingness") + 
  ylab("") + 
  qc_theme

site <- ggplot(data=site_QC) +
  geom_histogram(aes(x=F_MISS), bins=100, fill="#af8dc3") +
  geom_vline(xintercept = 0.1, linetype="dashed", color="grey50") +
  scale_y_continuous(expand=c(0,0), limits=c(0,3e+06),labels = function(x) format(x, scientific = TRUE)) +
  scale_x_continuous(expand=c(0,0), limits=c(0,1)) +
  xlab("Per-site missingness") + 
  ylab("") + 
  qc_theme

# Merge individual figures
plot_grid(QD,FS,MQ,MQRankSum,ReadPosRankSum,SOR,sample, site,align = "v", nrow = 8, labels = c('A', 'B','C','D','E','F','G','H'))
```
