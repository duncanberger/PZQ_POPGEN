# Supplementary figure 14: Selection statistics over candidate regions of selection
## Load themes
```{r}
selection_theme <- theme(
  legend.position="none",
  panel.grid = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y=element_text(face="bold", color="black"),
  axis.title.y = element_text(face="bold", color="black", size=9),
  axis.title.x = element_blank(),
  axis.ticks.x=element_blank())
```
## Load data
```{r}
# Load FST values 
fst_sites <- read.table("MAYUGE_TORORO_sites.weir.fst", header=TRUE)
# Load bed file of gene coordinates 
gff <- read.table("genes.short.txt", header=FALSE)
# Load per-site Pi values
pi_MAY <- read.table("MAYUGE_SITE_PI.sites.pi", header=TRUE)
pi_TOR <- read.table("KOCOGE_SITE_PI.sites.pi", header=TRUE)
# Load coverage in 2kb windows
site_cov_ALL <- read.table("all.name.cov", header=FALSE)
# Load median per sample coverage
sample_cov <- read.table("med.cov.txt", header=TRUE)
# Load metadata
scmeta <- read.table("schools.list", header=FALSE)
```
## Merge datasets
```{r}
merged_cov <- merge(site_cov_ALL,sample_cov, by.x=c("V7"), by.y=c("ID"))
merged_cov2 <- merge(merged_cov,scmeta, by.x=c("V7"), by.y=c("V2"))
pi_site_both <- merge(pi_MAY,pi_TOR, by.x=c("CHROM", "POS"), by.y=c("CHROM", 'POS'), all.x=FALSE, all.y=FALSE)
```
## Group data
```{r}
df.cov.summary <- merged_cov2 %>%
  group_by(V1.x,V2,V3,V1.y) %>%
  dplyr::summarise(
    sd = sd((V4/MC), na.rm = TRUE),
    len = median((V4/MC))
  )
df.cov.summary
```
## Plot FST
```{r}
fst_sites$COLOR <- ifelse(fst_sites$CHROM == "SM_V7_1","A",
                          ifelse(fst_district_2$CHROM == "SM_V7_3","A",
                                 ifelse(fst_district_2$CHROM == "SM_V7_5","A",
                                        ifelse(fst_district_2$CHROM == "SM_V7_7","A", no="B"))))

F_SITES_2 <- ggplot() + 
  geom_point(data=subset(fst_sites,CHROM=="SM_V7_4" & POS>2000000 & POS<3000000), 
             aes(x=POS/1000000, y=WEIR_AND_COCKERHAM_FST, color=COLOR), size=0.1,se=FALSE) +
  scale_x_continuous(expand=c(0,0), limits=c(2, 3)) +
  scale_y_continuous(expand=c(0,0), limits=c(0, 1)) + theme_bw() + PCA_theme +
  labs(y=expression(bold(F[ST]))) +
  scale_color_manual(values = c("grey40")) +
  theme(
    panel.grid = element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=7, color="black", face="bold"),
    axis.title.x=element_blank())
```
## Plot XP-EHH
```{r}
xpehh_allvsko$COLOR <- ifelse(xpehh_allvsko$CHROM == "1","A",
                              ifelse(xpehh_allvsko$CHROM == "3","A",
                                     ifelse(xpehh_allvsko$CHROM == "5","A",
                                            ifelse(xpehh_allvsko$CHROM == "7","A", no="B"))))

XPEHH_SITES_2 <- ggplot() + 
  geom_point(data=subset(xpehh_allvsko,CHROM=="4" & pos>2000000 & pos<3000000), 
             aes(x=pos/1000000, y=normxpehh, color=COLOR), size=0.1,se=FALSE) +
  scale_x_continuous(expand=c(0,0), limits=c(2,3)) +
  scale_y_continuous(expand=c(0,0), limits=c(-15, 15)) + theme_bw() + PCA_theme +
  labs(y=expression(bold("XP-EHH"))) +
  scale_color_manual(values = c("grey40")) +
  theme(
    panel.grid = element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=7, color="black", face="bold"),
    axis.title.x=element_blank())
```
## Plot iHS (Mayuge district)
```{r}
ihs_noko$COLOR <- ifelse(ihs_noko$V1 == "1","A",
                         ifelse(ihs_noko$V1 == "3","A",
                                ifelse(ihs_noko$V1 == "5","A",
                                       ifelse(ihs_noko$V1 == "7","A", no="B"))))
ihs_noko$COLOR <- ifelse(ihs_noko$V1=="4" & ihs_noko$V3>2168001 & ihs_noko$V3<2732001, "C","D")

	
IHS_SITES_2 <- ggplot() + 
  geom_point(data=subset(ihs_noko,V1=="4" & V3>2000000 & V3<3000000), 
             aes(x=V3/1000000, y=abs(V8), color=COLOR), size=0.1,se=FALSE) +
  scale_x_continuous(expand=c(0,0), limits=c(2,3)) +
  scale_y_continuous(expand=c(0,0), limits=c(0, 10)) + theme_bw() + PCA_theme +
  labs(y=expression(bold("|iHS| (Mayuge)"))) +
  scale_color_manual(values = c("#ef3b2c","grey40","grey40")) +
  theme(
    panel.grid = element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=7, color="black", face="bold"),
    axis.title.x=element_blank())
```
## Plot iHS (Tororo district)
```{r}
ihs_ko$COLOR <- ifelse(ihs_ko$V1 == "1","A",
                       ifelse(ihs_ko$V1 == "3","A",
                              ifelse(ihs_ko$V1 == "5","A",
                                     ifelse(ihs_ko$V1 == "7","A", no="B"))))

IHS_SITES_2B <- ggplot() + 
  geom_point(data=subset(ihs_ko,V1=="4" & V3>2000000 & V3<3000000), 
             aes(x=V3/1000000, y=abs(V8), color=COLOR), size=0.1,se=FALSE) +
  scale_x_continuous(expand=c(0,0), limits=c(2, 3)) +
  scale_y_continuous(expand=c(0,0), limits=c(0, 10)) + theme_bw() + PCA_theme +
  labs(y=expression(bold("|iHS| (Tororo)"))) +
  scale_color_manual(values = c("grey40","grey40")) +
  theme(
    panel.grid = element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=7, color="black", face="bold"),
    axis.title.x=element_blank())
```
## Plot Pi Mayuge / Pi Tororo
```{r}
PI_BOTH_2 <- ggplot() + 
  geom_point(data=subset(pi_site_both,CHROM=="SM_V7_4" & POS>2000000 & POS<3000000), 
             aes(x=POS/1000000, y=log10(PI.x/PI.y)), color="grey40", size=0.1) +
  scale_x_continuous(expand=c(0,0), limits=c(2,3)) +
  scale_y_continuous(expand=c(0,0), limits=c(-3, 3)) + 
  theme_bw() + PCA_theme +
  labs(y=expression(bold(log[10]("\U03C0"[Mayuge]*""*"/"*""*"\U03C0"[Tororo])))) +
  scale_color_manual(values = c("grey40")) +
  theme(
    panel.grid = element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=7, color="black", face="bold"),
    axis.title.x=element_blank())
```
## Plot coverage across regions (in 2 kb windows)
```{r}
cov_site_2 <- ggplot() + 
  geom_pointrange(data=subset(df.cov.summary,V1.x=="SM_V7_4" & V2>2000000 & V2<3000000), 
                  aes(x=V2/1000000, ymax=(len+sd), ymin=(len-sd), y=len), size=0.1, shape=21, fatten = 0.1, color="grey40") +
  scale_x_continuous(expand=c(0,0), limits=c(2,3)) +
  scale_y_continuous(expand=c(0,0), limits=c(-30, 30)) + 
  theme_bw() + PCA_theme +
  labs(y=expression(bold("Relative coverage"))) +
  scale_color_manual(values = c("#ef3b2c","grey40")) +
  theme(
    panel.grid = element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=7, color="black", face="bold"),
    axis.title.x=element_blank())
```
## Plot genes
```{r}
GENES_SITES_2 <- ggplot(data=subset(gff,V1=="SM_V7_4" & V2>2000000 & V2<3000000), 
                        aes(xmin = V2/1000000, xmax = V3/1000000, y=1, forward=V4)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm"), fill="grey75", color="black") +
  scale_x_continuous(expand=c(0,0), limits=c(2, 3)) +
  scale_y_continuous(expand=c(0,0), limits=c(0.9899, 1.0100001)) +
  theme_bw() +
  ylab("") +
  xlab("Position on Chromosome 4") +
  theme(
    legend.position="none",
    panel.grid = element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank(),
    axis.text.x=element_text(face="bold", color="black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face="bold", color="black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    panel.background = element_blank())
```
## Merge plots
```{r}
plot_grid(F_SITES_2,XPEHH_SITES_2,IHS_SITES_2,IHS_SITES_2B,PI_BOTH_2,cov_site_2,GENES_SITES_2, ncol=1, align="v")
```




