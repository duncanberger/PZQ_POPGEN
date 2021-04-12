# Supplementary figure 12: Nucleotide diversity for each district population
## Set ggplot themes and palettes
```{r}
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
## Import data (Pi values and a file of 2 kb windows across each chromosome (to account for any potential missing windows in the Pi files))
```{r}
PI_TORORO <- read.table("TORORO_PI.windowed.pi", header=TRUE, na.strings = "nan")
PI_MAYUGE <- read.table("MAYUGE_PI.windowed.pi", header=TRUE, na.strings = "nan")
window <- read.table("sm_2000.bed", header=FALSE)
names(window) <- c("CHROM","BIN_START", "BIN_STOP")
```
## Merge, reorder and calculate differences
```{r}
PI_merged0 <- merge(PI_MAYUGE, PI_TORORO, by.x = c("CHROM","BIN_START"), by.y = c("CHROM","BIN_START"), all.x = TRUE, all.y = TRUE)
PI_merged <- merge(PI_merged0, windows,by.x = c("CHROM","BIN_START"), by.y = c("CHROM","BIN_START"), all.x = TRUE, all.y = TRUE)
axisdf = PI_merged %>% group_by(CHROM) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
```
## Plot Pi in 2 kb windows for each subpopulation
```{r}
PI_M <- ggplot(subset(PI_merged,N_VARIANTS.x>20), aes(x=ID, y=(PI.x))) +
  geom_point( aes(color=as.factor(CHROM)),size=0.0001) +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center, expand=c(0,0)) +
  labs(y=expression(bold("\U03C0"[Mayuge]))) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.05), breaks=c(0,0.01,0.02,0.03,0.04,0.05)) +
  theme_bw() +
  selection_theme 
PI_T <- ggplot(subset(PI_merged,N_VARIANTS.x>20), aes(x=ID, y=(PI.y))) +
  geom_point( aes(color=as.factor(CHROM)),size=0.0001) +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center, expand=c(0,0)) +
  labs(y=expression(bold("\U03C0"[Tororo]))) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.05), breaks=c(0,0.01,0.02,0.03,0.04,0.05)) +
  theme_bw() +
  selection_theme
```
## Plot the ratio of Pi in 2 kb windows between each subpopulation
```{r}
PI_TM <- ggplot(subset(PI_merged,N_VARIANTS.x>20), aes(x=ID, y=log10(PI.x/PI.y))) +
  geom_point( aes(color=as.factor(CHROM)),size=0.0001) +
  labs(y=expression(bold(log[10]("\U03C0"[Mayuge]*" "*"/"*" "*"\U03C0"[Tororo])))) +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(-3,3), breaks=c(-3,-2,-1,0,1,2,3)) +
  theme_bw() +
  selection_theme
```
## Combine plots
```{r}
plot_grid(PI_M,PI_T,PI_TM,align = "v", rel_heights = c(1,1,1), nrow=3, labels=c("A","B","C"))
```
