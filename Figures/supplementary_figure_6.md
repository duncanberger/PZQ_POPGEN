# Supplementary Figure 6: Uncertainty in inference of demographic history using SMC++ for *Schistosoma mansoni* populations. 

## Import data (SMC++ estimates of demographic history)
```{r}
KOCOGE_SMC <- read.csv("KOCOGE_SE.csv", header=TRUE)
MUSUBI_SMC <- read.csv("MUSUBI_SE.csv", header=TRUE)
BULOOSI_SMC <- read.csv("BULOOSI_SE.csv", header=TRUE)
WALUKUBA_SMC <- read.csv("WALUKUBA_SE.csv", header=TRUE)
BWONDHA_SMC <- read.csv("BWONDHA_SE.csv", header=TRUE)
KENYA_SMC <- read.csv("KENYA_SE.csv", header=TRUE)
SENEGAL_SMC <- read.csv("SENEGAL_SE.csv", header=TRUE)
CAMEROON_SMC <- read.csv("CAMEROON_SE.csv", header=TRUE)
BUGOTO_SMC <- read.csv("BUGOTO_SE.csv", header=TRUE)
```
## Plot data for each subpopulation
```{r}
KO_SE_SMC <- ggplot(data=subset(KOCOGE_SMC)) + 
  geom_line(aes(x=(x),y=(y), color=label,linetype=label), size=0.5, alpha=0.7) +
  scale_x_log10(expand=c(0,0), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(expand=c(0,0),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) + scale_color_manual(values=c("grey","#CC79A7")) +
  scale_linetype_manual(values=c("solid", "solid","solid","twodash","dotted","solid","solid","solid","solid"))+xlab("Years before sampling") + labs(y=expression(bolditalic(N[e]))) +
  coord_cartesian(xlim = c(10^2, 10^5), ylim= c(10^2, 10^6)) + 
  theme_bw() + PCA_theme + theme(legend.position="none")  +
  annotation_logticks(sides="lb",outside = FALSE,short = unit(0.05, "cm"),mid = unit(0.1, "cm"),long = unit(0.1,"cm"))

MU_SE_SMC <- ggplot(data=subset(MUSUBI_SMC)) + 
  geom_line(aes(x=(x),y=(y), color=label,linetype=label), size=0.5, alpha=0.7) +
  scale_x_log10(expand=c(0,0), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(expand=c(0,0),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) + scale_color_manual(values=c("grey","#E69f00")) +
  scale_linetype_manual(values=c("solid", "solid","solid","twodash","dotted","solid","solid","solid","solid"))+xlab("Years before sampling") + labs(y=expression(bolditalic(N[e]))) +
  coord_cartesian(xlim = c(10^2, 10^5), ylim= c(10^2, 10^6)) + 
  theme_bw() + PCA_theme + theme(legend.position="none")  +
  annotation_logticks(sides="lb",outside = FALSE,short = unit(0.05, "cm"),mid = unit(0.1, "cm"),long = unit(0.1,"cm"))

BUL_SE_SMC <- ggplot(data=subset(BULOOSI_SMC)) + 
  geom_line(aes(x=(x),y=(y), color=label,linetype=label), size=0.5, alpha=0.7) +
  scale_x_log10(expand=c(0,0), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(expand=c(0,0),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) + scale_color_manual(values=c("grey","#0072b2")) +
  scale_linetype_manual(values=c("solid", "solid","solid","twodash","dotted","solid","solid","solid","solid"))+xlab("Years before sampling") + labs(y=expression(bolditalic(N[e]))) +
  coord_cartesian(xlim = c(10^2, 10^5), ylim= c(10^2, 10^6)) + 
  theme_bw() + PCA_theme + theme(legend.position="none")  +
  annotation_logticks(sides="lb",outside = FALSE,short = unit(0.05, "cm"),mid = unit(0.1, "cm"),long = unit(0.1,"cm"))

WA_SE_SMC <- ggplot(data=subset(WALUKUBA_SMC)) + 
  geom_line(aes(x=(x),y=(y), color=label,linetype=label), size=0.5, alpha=0.7) +
  scale_x_log10(expand=c(0,0), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(expand=c(0,0),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) + scale_color_manual(values=c("grey","#D55E00")) +
  scale_linetype_manual(values=c("solid", "solid","solid","twodash","dotted","solid","solid","solid","solid"))+xlab("Years before sampling") + labs(y=expression(bolditalic(N[e]))) +
  coord_cartesian(xlim = c(10^2, 10^5), ylim= c(10^2, 10^6)) + 
  theme_bw() + PCA_theme + theme(legend.position="none")  +
  annotation_logticks(sides="lb",outside = FALSE,short = unit(0.05, "cm"),mid = unit(0.1, "cm"),long = unit(0.1,"cm"))

BW_SE_SMC <- ggplot(data=subset(BWONDHA_SMC)) + 
  geom_line(aes(x=(x),y=(y), color=label,linetype=label), size=0.5, alpha=0.7) +
  scale_x_log10(expand=c(0,0), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(expand=c(0,0),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) + scale_color_manual(values=c("grey","#009e73")) +
  scale_linetype_manual(values=c("solid", "solid","solid","twodash","dotted","solid","solid","solid","solid"))+xlab("Years before sampling") + labs(y=expression(bolditalic(N[e]))) +
  coord_cartesian(xlim = c(10^2, 10^5), ylim= c(10^2, 10^6)) + 
  theme_bw() + PCA_theme + theme(legend.position="none")  +
  annotation_logticks(sides="lb",outside = FALSE,short = unit(0.05, "cm"),mid = unit(0.1, "cm"),long = unit(0.1,"cm"))

KE_SE_SMC <- ggplot(data=subset(KENYA_SMC)) + 
  geom_line(aes(x=(x),y=(y), color=label,linetype=label), size=0.5, alpha=0.7) +
  scale_x_log10(expand=c(0,0), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(expand=c(0,0),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) + scale_color_manual(values=c("grey","black")) +
  scale_linetype_manual(values=c("solid", "solid","solid","twodash","dotted","solid","solid","solid","solid"))+xlab("Years before sampling") + labs(y=expression(bolditalic(N[e]))) +
  coord_cartesian(xlim = c(10^2, 10^5), ylim= c(10^2, 10^6)) + 
  theme_bw() + PCA_theme + theme(legend.position="none")  +
  annotation_logticks(sides="lb",outside = FALSE,short = unit(0.05, "cm"),mid = unit(0.1, "cm"),long = unit(0.1,"cm"))

SEN_SE_SMC <- ggplot(data=subset(SENEGAL_SMC)) + 
  geom_line(aes(x=(x),y=(y), color=label,linetype=label), size=0.5, alpha=0.7) +
  scale_x_log10(expand=c(0,0), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(expand=c(0,0),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) + scale_color_manual(values=c("grey","black")) +
  scale_linetype_manual(values=c("solid", "solid","solid","twodash","dotted","solid","solid","solid","solid"))+xlab("Years before sampling") + labs(y=expression(bolditalic(N[e]))) +
  coord_cartesian(xlim = c(10^2, 10^5), ylim= c(10^2, 10^6)) + 
  theme_bw() + PCA_theme + theme(legend.position="none")  +
  annotation_logticks(sides="lb",outside = FALSE,short = unit(0.05, "cm"),mid = unit(0.1, "cm"),long = unit(0.1,"cm"))

CAM_SE_SMC <- ggplot(data=subset(CAMEROON_SMC)) + 
  geom_line(aes(x=(x),y=(y), color=label,linetype=label), size=0.5, alpha=0.7) +
  scale_x_log10(expand=c(0,0), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(expand=c(0,0),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) + scale_color_manual(values=c("grey","black")) +
  scale_linetype_manual(values=c("solid", "solid","solid","twodash","dotted","solid","solid","solid","solid"))+xlab("Years before sampling") + labs(y=expression(bolditalic(N[e]))) +
  coord_cartesian(xlim = c(10^2, 10^5), ylim= c(10^2, 10^6)) + 
  theme_bw() + PCA_theme + theme(legend.position="none")  +
  annotation_logticks(sides="lb",outside = FALSE,short = unit(0.05, "cm"),mid = unit(0.1, "cm"),long = unit(0.1,"cm"))

BUG_SE_SMC <- ggplot(data=subset(BUGOTO_SMC)) + 
  geom_line(aes(x=(x),y=(y), color=label,linetype=label), size=0.5, alpha=0.7) +
  scale_x_log10(expand=c(0,0), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(expand=c(0,0),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) + scale_color_manual(values=c("grey","#56B4E9")) +
  scale_linetype_manual(values=c("solid", "solid","solid","twodash","dotted","solid","solid","solid","solid"))+xlab("Years before sampling") + labs(y=expression(bolditalic(N[e]))) +
  coord_cartesian(xlim = c(10^2, 10^5), ylim= c(10^2, 10^6)) + 
  theme_bw() + PCA_theme + theme(legend.position="none")  +
  annotation_logticks(sides="lb",outside = FALSE,short = unit(0.05, "cm"),mid = unit(0.1, "cm"),long = unit(0.1,"cm"))
```
## Combine plots
```{r}
plot <- plot_grid(BUG_SE_SMC,BW_SE_SMC,KO_SE_SMC,MU_SE_SMC,BUL_SE_SMC,WA_SE_SMC,
                  KE_SE_SMC,SEN_SE_SMC,CAM_SE_SMC, nrow=5,ncol=2, align="v",
                  labels = "AUTO")
```
