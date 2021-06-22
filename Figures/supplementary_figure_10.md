# Supplementary Figure 10: Genome-wide allele frequency patterns.
## Load and group SFS data
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
## Plot 1D-SFS
```{r}
d1sfs <- ggplot(data=subset(sfs, Bin>1)) + 
  geom_col(data=df.summary,aes(x=(Bin-1), y=(len/50000), fill=School), position = "dodge") +
  geom_point(aes(x=(Bin-1), y=(Value/50000), color=School),
             position=position_jitterdodge(dodge.width =1,jitter.width = 0.2), size=0.3) +
  geom_errorbar(data=df.summary,aes(ymin=(len/50000)-(sd/50000), ymax=(len/50000)+(sd/50000), x=(Bin-1), group=School), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values=c("#56B4E9", "#009e73","#CC79A7","#E69f00","#D55E00","#0072b2"), na.value="grey50") +
  scale_color_manual(values=c("grey50", "grey50","grey50","grey50","grey50","grey50"), na.value="grey50") +
  scale_y_continuous(limits=c(0,0.7), expand=c(0,0),breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)) + 
  scale_x_continuous(limits=c(0,31), expand=c(0,0), breaks=c(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30)) + 
  theme_bw() +
  theme(axis.ticks.x = element_line(),
        axis.title = element_text(face="bold", size=7),
        axis.text = element_text(face="bold", size=5)) + 
  xlab("Allele frequency") + 
  ylab("Proportion of alleles") +
  PCA_theme + theme(axis.text = element_text(color="black"))
```
## Load and group 1D-SFS residuals
```{r}
sfs_res <-read.csv("sfs_res.csv", header=FALSE)
df.summary_res <- sfs_res %>%
  group_by(V2, V3) %>%
  subset(V2>1) %>%
  dplyr::summarise(
    sd = sd((V1), na.rm = TRUE),
    len = median(V1)
  )
df.summary_res
```
## Plot 1D-SFS residuals
```{r}
d1sfs_res <- ggplot(data=subset(df.summary_res, V2>1)) + 
  geom_col(aes(x=(V2-1), y=len, fill=V3), position = "dodge") +
  geom_point(data=subset(sfs_res, V2>1),aes(x=(V2-1), y=(V1), color=V3),
             position=position_jitterdodge(dodge.width =1,jitter.width = 0.2), size=0.3) +
  scale_fill_manual(values=c("#56B4E9", "#009e73","#CC79A7","#E69f00","#D55E00","#0072b2"), na.value="grey50") +
  scale_color_manual(values=c("grey50", "grey50","grey50","grey50","grey50","grey50"), na.value="grey50") +
  scale_x_continuous(expand=c(0,0),limits=c(0,31), breaks=c(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30)) + 
  scale_y_continuous(limits=c(0,2), expand=c(0,0)) + 
  theme_bw() +
  theme(axis.ticks.x = element_line(),
        axis.title = element_text(face="bold", size=7, color="black"),
        axis.text = element_text(face="bold", size=5, color="black")) + 
  xlab("Allele frequency") + 
  ylab("Residuals") +
  PCA_theme
  ```
