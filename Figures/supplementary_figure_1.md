# Supplementary Figure 1: Depth of read coverage
```{r}
# Load libraries
library(ggplot2)
library(dplyr)

# Read in data
all.bedtools.cov <- read.table("median.coverage.txt", header=FALSE)

#Reorder windows
all.bedtools.cov<- all.bedtools.cov[order(all.bedtools.cov$V1, all.bedtools.cov$V2),]
all.bedtools.cov$ID <- seq.int(nrow(cov_wind))
axisdf = all.bedtools.cov %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )

# Plot
cov_windows <- ggplot(all.bedtools.cov, aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.01) +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous(label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,75), breaks=c(0,25,50,75))+
  theme_bw() + xlab("") +
  ylab("Median depth") +
  theme( legend.position="none",panel.grid = element_blank(), 
         axis.text.y=element_text(face="bold", color="black"),
         axis.text.x=element_text(face="bold", color="black"),
         axis.title.y = element_text(face="bold", color="black"), axis.ticks.x=element_blank())
```
