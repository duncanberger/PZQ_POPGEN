# Supplementary Figure 9: Identification of sex of all samples using differential depth of coverage.
## Import and reorder data
```{r}
sex <- read.table("sexing.csv", header=TRUE, sep=',')
sex <- sex[order(sex$ZSR.PAR),]
sex$ID <- seq.int(nrow(sex))
```
## Plot data
```{r}
ggplot(data=sex) +
  geom_col(aes(x=ID, y=as.numeric(ZSR.PAR), fill=Sex), color="white") +
  geom_hline(yintercept = 0.75, linetype="dashed", color="grey75") +
 scale_y_continuous(expand=c(0,0), limits=c (0,1.5))+
  scale_x_discrete(expand=c(0,0))+
  scale_color_manual(values=c("#540d6e","#ee4266","#ffd23f","#0ead69")) +
  scale_fill_manual(values=c("#540d6e","#ee4266","#ffd23f","#0ead69")) +
  xlab("Sample") + ylab("Depth of read coverage ratio (ZSR/PAR)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.y=element_text(face="bold", color="black"),
        axis.title.y = element_text(face="bold", color="black"),
        axis.title.x = element_text(face="bold", color="black"),
        axis.text.x=element_text(face="bold", color="black"))
```
