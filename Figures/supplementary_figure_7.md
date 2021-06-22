# Supplementary Figure 7: Inference of demographic history using SMC++ for *Schistosoma mansoni* subpopulations. 
# Setup <a name="setup"></a>
```{r}
# Load required packages
library("ggplot2")
library("scales")
```
# Load data
```
smc <- read.table("smcpp.csv", header=TRUE, sep="\t", check.names = FALSE, comment.char = "")
```
# Plot data
```
ggplot(data=smc) + 
  geom_line(aes(x=(x),y=(y), color=label,linetype=label), size=1) +
  scale_x_log10(expand=c(0,0), breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(expand=c(0,0),breaks = trans_breaks("log10", function(x) 10^x), 
              labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values=c("#56B4E9","#0072b2","#009e73","grey","grey","#CC79A7","#E69f00","grey","#D55E00","grey","grey")) +
  scale_linetype_manual(values=c("solid", "solid","solid","twodash","dotted","solid","solid","solid","solid"))+
  xlab("Years before present") + 
  labs(y=expression(bold("Effective population size "*bolditalic((N[e]))))) +
  coord_cartesian(xlim = c(10^2, 10^5), ylim= c(10^2, 10^6)) +
  theme_bw() + 
  PCA_theme + 
  theme(legend.position="none", legend.key.size = unit(0.5,"cm"), legend.title =element_blank())  +
  annotation_logticks(sides="lb",outside = FALSE,
                      short = unit(0.05, "cm"),
                      mid = unit(0.1, "cm"),
                      long = unit(0.1,"cm"))
 ```
