# Supplementary Figure 4: Population structure in *Schistosoma mansoni* samples from Uganda overlaid with host or drug treatment-phenotype details.
## Setup <a name="setup"></a>
```{r}
# Load required packages
library("ggplot2")
library("reshape2")

# Load metadata
key <- read.table("supplementary_table_10.txt", header=TRUE, sep="\t", check.names = FALSE, comment.char = "")
```
## Figure 2A: Principal component analysis
```{r}
# Create a theme and palette
PCA_theme <- theme(axis.title=element_text(face="bold",size=9),
                   panel.background=element_blank(),
                   legend.text = element_text(face="bold"),
                   legend.position = "none",
                   axis.text=element_text(face="bold"),
                   panel.border = element_rect(color="#4c4c4c",fill=NA),
                   panel.grid=element_blank(),
                   legend.title=element_blank())

# Load data (same input as supplementary figure 3)
eigenvec <- read.delim("prunedData.eigenvec", header=TRUE, sep="\t")
eigenval <- read.delim("prunedData.eigenval", sep="\t", header=FALSE)
eigenvec_merged <- (merge(key, eigenvec, all=TRUE, by.y = "IID", by.x='sample_ID'))

# Calculate contribution of each eigenvalue to total variance
sum_eigenval<-sum(eigenval$V1)
sum_eigenval<-lapply(eigenval$V1,function(x){
  rt<-(x/sum_eigenval)*100
  rt<-round(rt)
  return(rt) 
})
```
### Plot supplementary figure 4A
```{r}
# Plot first two principal components, color points by host
sfig_4_A1 <- ggplot(eigenvec_merged, aes((PC1),(PC2))) + 
  geom_point(size=0.75, aes(color=sample_ID, fill=sample_ID)) +
  xlab(paste0("PC1 (",sum_eigenval[[1]],"%)")) + 
  ylab(paste0("PC2 (",sum_eigenval[[2]],"%)")) + 
  scale_shape_manual(values=c(21,23,22,24)) +
  scale_x_continuous(limits=c(-0.100000001,0.400000001), expand=c(0,0)) + 
  scale_y_continuous(limits=c(-0.200000001,0.200000001), expand=c(0,0)) +
  theme_bw() + PCA_theme 

# Plot first third and fourth principal components, color points by host
sfig_4_A2 <- ggplot(eigenvec_merged, aes((PC3),(PC4))) + 
  geom_point(size=0.75, aes(color=sample_ID, fill=sample_ID)) +
  xlab(paste0("PC3 (",sum_eigenval[[3]],"%)")) + 
  ylab(paste0("PC4 (",sum_eigenval[[4]],"%)")) + 
  scale_color_manual(values=pca_palette)+
  scale_fill_manual(values=pca_palette)+
  scale_shape_manual(values=c(21,23,22,24)) +
  scale_x_continuous(limits=c(-0.400000001,0.600000001), expand=c(0,0)) + 
  scale_y_continuous(limits=c(-0.2000000001,0.800000001), expand=c(0,0)) +
  theme_bw() + PCA_theme
```
### Plot supplementary figure 4B
```{r}
# Plot first two principal components, color points by sampling point
sfig_4_B1 <- ggplot(eigenvec_merged, aes((PC1),(PC2))) + 
  geom_point(size=0.75, aes(color=sample_point, fill=sample_point)) +
  xlab(paste0("PC1 (",sum_eigenval[[1]],"%)")) + 
  ylab(paste0("PC2 (",sum_eigenval[[2]],"%)")) + 
  scale_color_manual(values=pca_palette)+
  scale_fill_manual(values=pca_palette)+
  scale_shape_manual(values=c(21,23,22,24)) +
  scale_x_continuous(limits=c(-0.100000001,0.400000001), expand=c(0,0)) + 
  scale_y_continuous(limits=c(-0.200000001,0.200000001), expand=c(0,0)) +
  theme_bw() + PCA_theme 

# Plot first third and fourth principal components, color points by sampling point
sfig_4_B2 <- ggplot(eigenvec_merged, aes((PC3),(PC4))) + 
  geom_point(size=0.75, aes(color=sample_point, fill=sample_point)) +
  xlab(paste0("PC3 (",sum_eigenval[[3]],"%)")) + 
  ylab(paste0("PC4 (",sum_eigenval[[4]],"%)")) + 
  scale_color_manual(values=pca_palette)+
  scale_fill_manual(values=pca_palette)+
  scale_shape_manual(values=c(21,23,22,24)) +
  scale_x_continuous(limits=c(-0.400000001,0.600000001), expand=c(0,0)) + 
  scale_y_continuous(limits=c(-0.2000000001,0.800000001), expand=c(0,0)) +
  theme_bw() + PCA_theme
```
### Merge
```{r}
# Merge figures
plot_grid(sfig_4_A1,sfig_4_A2,sfig_4_B1,sfig_4_B2, ncol=2, nrow=2, rel_widths=c(1,1), rel_heights = c(1,0.6)) + 
```
