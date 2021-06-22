# Figure 3: Genome-wide allele frequency patterns. 
## Figure 3A: 1D-SFS for each school subpopulation
### Load and group SFS data
```{r}
sfs <-read.csv("sfs.csv", header=TRUE)
df.summary <- sfs %>%
  group_by(Bin, School) %>%
  subset(Bin>1) %>%
  dplyr::summarise(
    sd = sd((Value), na.rm = TRUE),
    len = mean(Value)
  )
df.summary
```
### Plot 1D-SFS
```{r}
d1sfs <- ggplot(data=subset(df.summary, Bin>1)) + 
  geom_col(aes(x=(Bin-1), y=len/50000, fill=School), position = "dodge") +
  scale_fill_manual(values=c("#56B4E9", "#009e73","#CC79A7","#E69f00","#D55E00","#0072b2"), na.value="grey50") +
  geom_errorbar(aes(ymin=(len/50000)-(sd/50000), ymax=(len/50000)+(sd/50000), x=(Bin-1), group=School), width=.2, position=position_dodge(.9)) +
  scale_y_continuous(limits=c(0,0.8), expand=c(0,0)) + 
  scale_x_continuous(limits=c(0,31), expand=c(0,0), breaks=c(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30)) + 
  theme_bw() +
  theme(axis.ticks.x = element_line(),
        axis.title = element_text(face="bold", size=10)) + 
  xlab("Allele frequency") + 
  ylab("Proportion of alleles") +
  PCA_theme
```
## Figure 3B: Tajima’s D values for each school population 
### Import and order data
```{r}
tajima_d <- read.table("TD.all.txt", header=FALSE)
tajima_d$V5 = factor(tajima_d$V5, levels=c('Bugoto','Bwondha','Musubi','Kocoge'))
```
### Randomly subset data (for plotting)
```{r}
tajima_d_subset <- subset(tajima_d, V3>20)
tajima_d_subset2 <- sample_n(tajima_d, 20000)
```
### Plot boxplot
```{r}
tajima_d_subset_plot <- ggplot(data=tajima_d_subset2, aes(x=V5, y=(V4), fill=V5, color=V5)) + 
  geom_point(position=position_jitterdodge(dodge.width =1,jitter.width = 0.8),alpha=0.3, size=0.0001) +
  geom_boxplot(aes(fill=V5),outlier.alpha = 0.0,notch = TRUE, outlier.colour = "grey35", color="black", alpha=0, width=0.275) +
  theme_bw()+
  scale_fill_manual(values=c("#56B4E9", "#009e73","#E69f00","#CC79A7")) +
  scale_color_manual(values=c("#56B4E9", "#009e73","#E69f00","#CC79A7")) +
  scale_y_continuous(expand=c(0,0), limits=c(-4,4), breaks=c(-4.0,-2.0,0,2.0,4.0)) +
  xlab("A") +
  labs(y=expression(bold("Median Tajima's D"))) +
  PCA_theme + theme(legend.position = "none") +
  theme(legend.text = element_text(face = "bold"), axis.ticks = element_line())
```
### Plot density
```{r}
dens <- ggplot(data=tajima_d_subset2, aes(x=(V4), fill=V5, color=V5)) +
  geom_density(alpha=0.5) + 
  theme_nothing()+
  xlab("") + ylab("") +
  scale_fill_manual(values=c("#56B4E9", "#009e73","#E69f00","#CC79A7")) +
  scale_color_manual(values=c("#56B4E9", "#009e73","#E69f00","#CC79A7"))  +
  scale_x_continuous(expand=c(0,0), limits=c(-4,4), breaks=c(-4.0,-2.0,0,2.0,4.0)) +
  PCA_theme + theme(legend.position = "none") +
  theme(legend.text = element_text(face = "bold"),
        axis.text =element_blank(), panel.border = element_blank()) + coord_flip()
```
### Merge plots
```{r}
td_sum <- plot_grid(tajima_d_subset_plot,dens, nrow=1, align="h", rel_widths = c(1,0.2))
```
