# Supplementary Figure 19: The decay of linkage disequilibrium as the squared correlation of allele frequency (r^2) against distance (bp)
```{r}
# Load libraries
library(dplyr)
library(ggplot2)

# Load data
ld_tororo <- read.table("kocoge_median.ld.txt", header=FALSE)
ld_mayuge <- read.table("mayuge_median.ld.txt", header=FALSE)

# Plot for each district
tororo <- ggplot(data=ld_tororo) + 
  geom_smooth(aes(x=V2/1000,y=V3,fill=V1, color=V1)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks=c(0,0.2,0.4,0.6,0.8,1)) +
  scale_x_continuous(expand=c(0,0), limits=c(0,50)) + 
  xlab("Distance (kb)") +
  ylab(expression(bold("Median r"^"2")))+
  theme_bw() + PCA_theme
  
mayuge <- ggplot(data=ld_mayuge) + 
  geom_smooth(aes(x=V2/1000,y=V3,color=V1, fill=V1) ) +
  scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks=c(0,0.2,0.4,0.6,0.8,1)) +
  scale_x_continuous(expand=c(0,0), limits=c(0,50)) + 
  xlab("Distance (kb)") +
  ylab(expression(bold("Median r"^"2")))+
  theme_bw() + PCA_theme

legend <- g_legend(mayuge) 

# Merge plots
plot_grid(mayuge, tororo, legend, nrow=1, align = "h", rel_heights = c(1,1,0.2)) + theme(legend.position = "none")
```
