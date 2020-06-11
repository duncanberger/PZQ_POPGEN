library(dplyr)
reshape2 ggplot

pid <- read.table("p_id.csv", sep=",", header=TRUE)
pid <- read.table("supplementary_table_1.txt", sep="\t", header=TRUE)

pid_2 <- pid %>%
  select('patient_ID','School','mean_post.treatment_posterior_egg_reduction_rate','pre.treatment_miracidia_sequenced_passed_qc','post.treatment_miracidia_sequenced_passed_qc')

pid_2_melt <- melt(pid_2,id.vars = c("patient_ID","mean_post.treatment_posterior_egg_reduction_rate","School"), measure.vars=c("pre.treatment_miracidia_sequenced_passed_qc","post.treatment_miracidia_sequenced_passed_qc"))
pid_2_melt$variable <- ifelse(pid_2_melt$variable == "pre.treatment_miracidia_sequenced_passed_qc", "0", "27")

pid_2_melt$N_jit <- jitter(pid_2_melt$value, factor=0.5)
pid_2_melt$Days_jit <- jitter(as.numeric(pid_2_melt$variable), factor=0.5)

pca_palette <- c("#56B4E9", "#009e73","#E69f00","#CC79A7")
pid$School = factor(pid$School, levels=c('Bugoto','Bwondha','Musubi','Kocoge'))


ggplot(data=subset(pid_2_melt, value!=0), aes(x=Days_jit, group=patient_ID, y=value)) + 
  geom_line(aes(group=patient_ID), alpha=0.6) +
  geom_point(aes(size=value, fill=School), alpha=0.6,color="black",pch=21) +
  facet_grid(cols=vars(School),scales="free_y", space = "free_y") +
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
