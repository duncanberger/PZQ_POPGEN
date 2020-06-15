# Coefficient of variation from ADMIXTURE
```{r}
# Load libraries
library(dplyr)
library(ggplot2)

# Read in data
cv_scores <- read.table("cv_scores.txt", header=FALSE)

# Calculate, mean, median, and standard deviations
cv_scores.summary <- cv_scores %>%
  group_by(V1) %>%
  summarise(
    med = median(V2),
    sd = sd(V2, na.rm = TRUE),
    len = mean(V2)
  )

# Plot
ggplot(data=cv_scores.summary) + 
  geom_point(aes(x=V1, y=med)) +
  geom_errorbar(aes(x = V1, y = len, ymin = med-sd, ymax = med+sd)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,1)) +
  scale_x_continuous(expand=c(0,0), limits=c(0,21), breaks=c(0,2,4,6,8,10,12,14,16,18,20)) +
  xlab("Number of clusters (K)") + 
  ylab("Coefficient of Variation Error") +
  theme( legend.position="none",panel.grid = element_blank(), 
         axis.text.y=element_text(face="bold", color="black", size=8),
         axis.title.y = element_text(face="bold", color="black", size=10),
         axis.text.x=element_text(face="bold", color="black", size=8),
         axis.title.x = element_text(face="bold", color="black", size=10),  
         panel.border = element_rect(color="black",fill = NA),
         panel.background = element_blank())
```
