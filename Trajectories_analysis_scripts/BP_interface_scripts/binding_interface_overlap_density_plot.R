# Libraries
library(ggplot2)
library(dplyr)
library(forcats)
png(file="./DNA_overlap_ratio_5_density_plot_displacement.png",width=10,height=6,units="in",res=300)
# Create data
df <-read.table("./Nucleosome_contacts/overlap_ratio_5_displacement.csv",header=FALSE, sep=",", col.names = c("name","ratio"))
#df["PDB"] = rep("PDB", 43)
#data$ratio %>%
#  mutate(name = fct_reorder(name, ratio)) %>%
# Horizontal version
  ggplot(data =df, aes(x=ratio ))   +
#  geom_segment( aes(xend=name, yend=0),lwd=0,color="Blue") +
  geom_density(fill="#69b3a2", bw = "nrd0") + 
  theme_bw(base_line_size =0 ) + 
  xlab("DNA binding interface overlap ") +
  ylab("Density of structure number") +
  xlim(c(0,1)) +
#  ggtitle("Title") 
  theme(axis.title.x=element_blank()) + theme(panel.border = element_rect(size =1.5, fill = NA)) +
  theme(axis.title.y=element_text(angle=90, vjust=4)) +
  theme(plot.title=element_text(size=16, vjust=3)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(
    axis.title.x = element_text(size = 28,face="bold"),
    axis.text.x = element_text(size = 28,face="bold"),
    axis.title.y = element_text(size = 28,face="bold"),
    axis.text.y = element_text(size = 28,face="bold"))


dev.off()