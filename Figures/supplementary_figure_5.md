# Diversity of host infrapopulations
```{r}
# Load libraries
library(dplyr)
library(ggplot2)
library(stringr)

# Load data
all_infra_pi1 <- read.table("pi.host.txt", header=TRUE)

# Create a new column indicating school population
all_infra_pi1$School <- str_sub(all_infra_pi1$pop, -5,-2)

# Group by host, remove 2 hosts due to consistently low coverage 
all_infra_pi2 <- all_infra_pi1 %>% group_by(pop) %>% sample_n(size = 1500)

# Remove 2 hosts due to consistently low coverage 
all_infra_pi3 <- all_infra_pi2[(all_infra_pi2$pop!="Bu4" & all_infra_pi2$pop!="Bu5"),]

# Order by school
all_infra_pi3$Site = factor(all_infra_pi3$School, levels=c('Bugoto','Bwondha','Musubi','Kocoge'))

# Plot
ggplot(data=all_infra_pi3, aes(x=pop, y=log10(avg_pi), fill=School, color=School)) + 
  geom_point(position=position_jitterdodge(dodge.width =1,jitter.width = 0.8),alpha=0.3, size=0.1) +
  geom_boxplot(aes(fill=School),outlier.alpha = 0.0,notch = TRUE, outlier.colour = "grey35", color="black", alpha=0, width=0.275) +
  theme_bw()+
  scale_fill_manual(values=c("#56B4E9", "#009e73","#E69f00","#CC79A7")) +
  scale_color_manual(values=c("#56B4E9", "#009e73","#E69f00","#CC79A7")) +
  scale_y_continuous(expand=c(0,0), limits=c(-50,50000), breaks=c(-5.0,-4.0,-3.0,-2.0,-1.0,0)) +
  xlab("Host") +
  coord_cartesian(ylim = c(-5, 0)) +
  labs(y=expression(bold(-log[10]*("\U03C0")))) +
  PCA_theme + theme(legend.position = "none") +
  theme(legend.text = element_text(size=6.5, face = "bold"))
```
