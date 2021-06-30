# Figure 1B: Sample summary.
```{r}
# Load packages
library(dplyr)
library(ggplot2)
library(reshape2)

# Load data
pid <- read.table("supplementary_data_9.txt", sep="\t", header=TRUE)

# Create color palette
pca_palette <- c("#56B4E9", "#009e73","#E69f00","#CC79A7")

#Select columns of interest
pid_2 <- pid %>% select('patient_ID','School','mean.post.treatment.posterior.egg.reduction.rate','pre.treatment.miracidia.sequenced.passed.qc','post.treatment.miracidia.sequenced.passed.qc')

# Convert table to ggplot friendly format
pid_2_melt <- melt(pid_2,id.vars = c("child_ID","mean_post.treatment.posterior.egg.reduction.rate","School"), measure.vars=c("pre.treatment.miracidia.sequenced.passed.qc","post.treatment.miracidia.sequenced.passed.qc"))

# Convert pre- and post-treatment numbers to numerical format
pid_2_melt$variable <- ifelse(pid_2_melt$variable == "pre.treatment.miracidia.sequenced.passed.qc", "0", "27")

# Slightly jitter the points so they plot better
pid_2_melt$N_jit <- jitter(pid_2_melt$value, factor=0.5)
pid_2_melt$Days_jit <- jitter(as.numeric(pid_2_melt$variable), factor=0.5)

# Order by school
pid_2_melt$School = factor(pid_2_melt$School, levels=c('Bwondha','Bugoto','Musubi','Kocoge'))

# Plot
ggplot(data=subset(pid_2_melt, value!=0), aes(x=Days_jit, group=child_ID, y=value)) + 
  facet_grid(.~School,scales="free_y", space = "free_y") +
  geom_line(aes(group=child_ID), alpha=0.6) +
  geom_point(aes(fill=School), alpha=0.6,color="black",pch=21, size=4) +
  scale_x_continuous(breaks=c(0,27)) +
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9), limits=c(0,9), expand=c(0,0))+
  scale_fill_manual(values=pca_palette) +
  ylab("Miracidia sampled") +
  xlab("") +
  theme_bw() + theme(axis.title=element_text(face="bold",size=9),
                     panel.background=element_blank(),
                     legend.text = element_text(face="bold"),
                     legend.position = "none",
                     axis.text=element_text(face="bold"),
                     panel.border = element_rect(color="#4c4c4c",fill=NA),
                     panel.grid=element_blank(),
                     legend.title=element_blank())


```
