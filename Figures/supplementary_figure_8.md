# Supplementary Figure 8: Diversity of infrapopulations.
## Load libraries
```{r}
library(dplyr)
library(ggplot2)
library(stringr)
library("ggpubr")
```
# Load data
```{r}
all_infra_pi1 <- read.table("pi.per_host.txt", header=TRUE)
age_ND <- read.table("ND_AGE.csv", sep=",", header=TRUE) # Not available due to prescence of indirect identifiers

# Create a new column indicating school population
#all_infra_pi1$School <- str_sub(all_infra_pi1$pop, -5,-2)

## Group by host, remove 2 hosts due to consistently low coverage 
all_infra_pi3 <- all_infra_pi1 %>% group_by(pop) %>% sample_n(size = 1500)

## Remove 2 hosts due to consistently low coverage 
#all_infra_pi3 <- all_infra_pi2[(all_infra_pi2$pop!="Bu4" & all_infra_pi2$pop!="Bu5"),]

## Order by school
all_infra_pi3$Site = factor(all_infra_pi3$School, levels=c('Bugoto','Bwondha','Musubi','Kocoge'))

## Plot per-host nucleotide diversity (colored by school)
A <- ggplot(data=all_infra_pi3, aes(x=pop, y=log10(avg_pi), fill=School, color=School)) + 
  geom_point(position=position_jitterdodge(dodge.width =1,jitter.width = 0.8),alpha=0.3, size=0.1) +
  geom_boxplot(aes(fill=School),outlier.alpha = 0.0,notch = TRUE, outlier.colour = "grey35", color="black", 
               alpha=0, width=0.275) +
  theme_bw()+
  scale_fill_manual(values=c("#56B4E9", "#009e73","#E69f00","#CC79A7")) +
  scale_color_manual(values=c("#56B4E9", "#009e73","#E69f00","#CC79A7")) +
  scale_y_continuous(expand=c(0,0), limits=c(-50,50000), breaks=c(-5.0,-4.0,-3.0,-2.0,-1.0,0)) +
  xlab("Child") +
  coord_cartesian(ylim = c(-5, 0)) +
  labs(y=expression(bold(-log[10]*("\U03C0")))) +
  PCA_theme + theme(legend.position = "none",
                    axis.text.x=element_text(face="bold", color="black", size=6),
                    axis.text.y=element_text(face="bold", color="black"),
                    axis.title.y = element_text(face="bold", color="black"),
                    axis.title.x = element_text(face="bold", color="black"))

## Plot per-host nucleotide diversity (colored by sex)        
D <- ggplot(data=all_infra_pi4, aes(x=pop, y=log10(avg_pi), fill=Sex, color=Sex)) + 
  geom_point(position=position_jitterdodge(dodge.width =1,jitter.width = 0.8),alpha=0.2, size=0.1) +
  geom_boxplot(aes(fill=Sex),outlier.alpha = 0.0,notch = TRUE, outlier.colour = "grey35", color="black", 
               alpha=0, width=0.275) +
  theme_bw()+
  scale_fill_manual(values=c("#540d6e","#ffd23f")) +
  scale_color_manual(values=c("#540d6e","#ffd23f")) +
  scale_y_continuous(expand=c(0,0), limits=c(-50,50000), breaks=c(-5.0,-4.0,-3.0,-2.0,-1.0,0)) +
  xlab("Child") +
  coord_cartesian(ylim = c(-5, 0)) +
  labs(y=expression(bold(-log[10]*("\U03C0")))) +
  PCA_theme + theme(legend.position = "none",
                    axis.text.x=element_text(face="bold", color="black", size=6),
                    axis.text.y=element_text(face="bold", color="black"),
                    axis.title.y = element_text(face="bold", color="black"),
                    axis.title.x = element_text(face="bold", color="black"))

## Plot median per-host nucleotide diversity (correlation to host age)
B <- ggscatter(age_ND, x = "Age", y = "Median.Pi", add = "reg.line", size=0) +
  stat_cor(label.y = 0.0035, label.x=10.5, 
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  scale_x_continuous(expand=c(0,0), limits=c(5,13))+
  scale_y_continuous(expand=c(0,0), limits=c(0.000,0.005))+
  xlab("Age (years)") +
  geom_point(aes(x=Age, y=Median.Pi, color=School)) +
  labs(y=expression(bold("Median \U03C0"))) +
  scale_color_manual(values=c("#56B4E9","#009e73","#CC79A7","#E69f00")) +
  theme_bw()+
  PCA_theme +
  theme(legend.position="none",
        panel.grid = element_blank(),
        axis.text.x=element_text(face="bold", color="black"),
        axis.text.y=element_text(face="bold", color="black"),
        axis.title.y = element_text(face="bold", color="black"),
        axis.title.x = element_text(face="bold", color="black"),
        axis.ticks.x=element_blank())

## Plot median per-host nucleotide diversity (correlation to infrapopulation size)
C <- ggscatter(age_ND, x = "Mira", y = "Median.Pi", add = "reg.line", size=0) +
  stat_cor(label.y = 0.001, label.x=5, 
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  scale_x_continuous(expand=c(0,0), limits=c(0,12))+
  scale_y_continuous(expand=c(0,0), limits=c(0.000,0.005))+
  xlab("Total miracidia sequenced") +
  geom_point(aes(x=Mira, y=Median.Pi, color=School)) +
  labs(y=expression(bold("Median \U03C0"))) +
  scale_color_manual(values=c("#56B4E9","#009e73","#CC79A7","#E69f00")) +
  theme_bw()+
  PCA_theme+
  theme(legend.position="none",
        panel.grid = element_blank(),
        axis.text.x=element_text(face="bold", color="black"),
        axis.text.y=element_text(face="bold", color="black"),
        axis.title.y = element_text(face="bold", color="black"),
        axis.title.x = element_text(face="bold", color="black"),
        axis.ticks.x=element_blank())

## Merge plots
E <- plot_grid(B,C,nrow=1, labels=c("C","D"))
plot_grid(A,D,E,nrow=3, labels=c("A","B",""))
```
## Quick t-test for significance (host sex vs nucleotide diversity)
```{r}
t.test(Median.Pi ~ Sex, data=subset(all_infra_pi4, pop!="Bu4" | pop!="Bu5"))
t.test(Median.Pi ~ Sex, data=subset(age_ND))
```
