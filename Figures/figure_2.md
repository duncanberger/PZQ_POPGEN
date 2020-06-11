# Figure 2: S. mansoni population structure in Uganda

0. [Project setup](#setup)
1. [Figure 2A](#figure2a)
2. [Figure 2B](#figure2b)
3. [Figure 2C](#figure2c)
4. [Figure 2D](#figure2d)
5. [Figure 2E](#figure2e)
## Project setup <a name="setup"></a>
```{r}
# Load required packages
library("ggplot2")
library("reshape2")
library("lattice")
library("cowplot")
library("SNPRelate")
library("ape")
library("ggtree")
library("phangorn")
# Load metadata
key <- read.table("supplementary_table_2.txt", header=TRUE, sep="\t", check.names = FALSE, comment.char = "")
```

## Figure 2A: Principal component analysis <a name="figure2a"></a>
```R
# Create theme
PCA_theme <- theme(axis.title=element_text(face="bold",size=9),
                   panel.background=element_blank(),
                   legend.text = element_text(face="bold"),
                   legend.position = "none",
                   axis.text=element_text(face="bold"),
                   panel.border = element_rect(color="#4c4c4c",fill=NA),
                   panel.grid=element_blank(),
                   legend.title=element_blank())

# Load data, produce in STEP X, and merge with metadata
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

# Plot first two principal components, color points by school
pc1_pc2 <- ggplot(eigenvec_merged, aes((PC1),(PC2))) + 
  geom_point(size=0.75, aes(color=school, fill=school)) +
  xlab(paste0("PC1 (",sum_eigenval[[1]],"%)")) + 
  ylab(paste0("PC2 (",sum_eigenval[[2]],"%)")) + 
  scale_color_manual(values=pca_palette)+
  scale_fill_manual(values=pca_palette)+
  scale_shape_manual(values=c(21,23,22,24)) +
  scale_x_continuous(limits=c(-0.100000001,0.400000001), expand=c(0,0)) + 
  scale_y_continuous(limits=c(-0.200000001,0.200000001), expand=c(0,0)) +
  theme_bw() + PCA_theme 

# Plot first third and fourth principal components, color points by school
pc3_pc4 <- ggplot(eigenvec_merged, aes((PC3),(PC4))) + 
  geom_point(size=0.75, aes(color=school, fill=school)) +
  xlab(paste0("PC3 (",sum_eigenval[[3]],"%)")) + 
  ylab(paste0("PC4 (",sum_eigenval[[4]],"%)")) + 
  scale_color_manual(values=c( "#56B4E9", "#009e73","#CC79A7","#E69f00"))+
  scale_fill_manual(values=c( "#56B4E9", "#009e73","#CC79A7","#E69f00"))+
  scale_shape_manual(values=c(21,23,22,24)) +
  scale_x_continuous(limits=c(-0.400000001,0.600000001), expand=c(0,0)) + 
  scale_y_continuous(limits=c(-0.2000000001,0.800000001), expand=c(0,0)) +
  theme_bw() + PCA_theme

                   
                   
                   
                   
                   
                   
                   
                   
                   