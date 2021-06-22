# Supplementary Figure 11: Median Tajima’s D values for each district population.

## Import Tajima's D data for both subpopulations
```{r}
taj_MAYUGE <- read.table("MAYUGE_TAJIMA_D.Tajima.D.2kb.txt", header=TRUE, na.strings = "nan")
taj_TORORO <- read.table("TORORO_TAJIMA_D.Tajima.D.2kb.txt", header=TRUE, na.strings = "nan")
```
## Merge and order data
```{r}
td_merged <- merge(taj_MAYUGE, taj_TORORO, by.x = c("CHROM","BIN_START"), by.y = c("CHROM","BIN_START"), all.x = TRUE, all.y = TRUE)
td_merged<- td_merged[order(td_merged$CHROM, td_merged$BIN_START),]
td_merged$ID <- seq.int(nrow(td_merged))
axisdf = td_merged %>% group_by(CHROM) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
```
## Plot for each population
```{r}
M <- ggplot(subset(td_merged,N_SNPS.x>20), aes(x=ID, y=(TajimaD.x))) +
  geom_point( aes(color=as.factor(CHROM)),size=0.0001) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(-6,6), breaks=c(-6,-4,-2,0,2,4,6)) +
  labs(y=expression(bold("Median Tajima's D"[Mayuge]))) +
  theme_bw() +
  selection_theme

T <- ggplot(subset(td_merged,N_SNPS.x>20), aes(x=ID, y=(TajimaD.y))) +
  geom_point( aes(color=as.factor(CHROM)),size=0.0001) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(-6,6), breaks=c(-6,-4,-2,0,2,4,6)) +
  labs(y=expression(bold("Median Tajima's D"[Tororo]))) +
  theme_bw() +
  selection_theme
```
## Plot the difference between populations
```{r} 
M_TM <- ggplot(subset(td_merged,N_SNPS.x>20), aes(x=ID, y=log10(TajimaD.x/TajimaD.y))) +
  geom_point( aes(color=as.factor(CHROM)),size=0.0001) +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(-4,4), breaks=c(-4,-2,0,2,4)) +
  labs(y=expression(bold("log"[10]("Tajima's D"[Mayuge]*"  /  "*"Tajima's D"[Tororo])))) +
  theme_bw() +
  selection_theme
```
## Combine all plots
```{r}
plot_grid(M,T,M_TM,align = "v", rel_heights = c(1,1,1,0.3), nrow=4)
```
