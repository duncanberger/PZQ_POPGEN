# Figure 3: S. mansoni population history
## Setup <a name="setup"></a>
```{r}
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
