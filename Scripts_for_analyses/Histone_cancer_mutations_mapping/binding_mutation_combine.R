# Libraries
library(ggplot2)
library(dplyr)
library(forcats)
png(file="./H4_binding_mutation_sites_msa.png",width=10,height=4,units="in",res=300)
# Create data
binding <-read.table("./H4_binding_sites_msa.txt",header=FALSE, sep=",", col.names = c("residue","number"))
mutation <-read.table("./H4_mutation_sites_msa.txt",header=FALSE, sep=",", col.names = c("residue","number"))

  ggplot(binding,aes(x=residue, y=number)) +
    geom_hline(yintercept = 0, size=1) +
  geom_segment(data = binding, aes(x=residue,xend=residue, y=0,yend=number),lwd=0.5,color="grey") +
  geom_point(data = binding, aes(x=residue, y=number), size=2, color="dark grey", 
             fill=alpha("dark turquoise", 0.3), alpha=0.7, shape=21, stroke=0.5) +

  geom_segment(data = mutation, aes(x=residue,xend=residue, y=0,yend=-number*2),lwd=0.5,color="grey") +
  geom_point(data = mutation, aes(x=residue, y=-number*2), size=2, color="dark grey",
             fill=alpha("hot pink", 0.3), alpha=0.7, shape=21, stroke=0.5) +  
  theme_classic() + 
    scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Number of mutations",
                                           breaks=c(0,-10,-20,-30), 
                                           labels=c(0,10,20,30)),
                       breaks=c(0,10,20,30,40,50), labels=c(0,10,20,30,40,50), limits = c(-30,30))  + 
  xlab("Residue Number") +
  ylab("Number of proteins") +

  theme(axis.title.x=element_blank()) + theme(panel.border = element_rect(size =1.5, fill = NA)) +
  theme(axis.title.y.right =element_text(angle=-90, vjust=3)) +
  theme(axis.title.y.left  =element_text(angle=90, vjust=3)) +
  theme(axis.title.x.bottom =element_text(angle=0, vjust=-2)) +
  theme(plot.title=element_text(size=18, vjust=4)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 20,face="bold"),
    axis.title.y.left = element_text(size = 24,face="bold",colour = "dark turquoise"),
    axis.title.y.right = element_text(size = 24,face="bold",colour = "hot pink"),
    axis.title.x.bottom = element_text(size = 23,face="bold",colour = "black"),
    axis.text.y = element_text(size = 20,face="bold"))


dev.off()