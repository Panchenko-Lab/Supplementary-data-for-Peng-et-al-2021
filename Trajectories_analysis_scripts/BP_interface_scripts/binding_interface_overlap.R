# Libraries
library(ggplot2)
library(dplyr)
library(forcats)
png(file="./DNA_overlap_ratio_3.png",width=10,height=6,units="in",res=300)
# Create data
data<-read.table("./overlap_ratio_3.csv",header=FALSE, sep=",", col.names = c("name","ratio"))
data %>%
  mutate(name = fct_reorder(name, ratio)) %>%
# Horizontal version
  ggplot(aes(x=name, y=ratio))   +
  geom_segment( aes(xend=name, yend=0),lwd=0,color="Blue") +
  geom_point( size=3, color="red") + coord_flip() + geom_line(size=0)+
  theme_bw(base_line_size =0 ) + 
  xlab("Nucleosome Binding Partner") +
  ylab("Binding Interface Overlaping Ratio ") +
  ylim(c(0,1)) +
#  ggtitle("Title") 
  theme(axis.title.x=element_blank()) + theme(panel.border = element_rect(size =1.5, fill = NA)) +
  theme(axis.title.y=element_text(angle=90, vjust=4)) +
  theme(plot.title=element_text(size=16, vjust=3)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(
    axis.title.x = element_text(size = 20,face="bold"),
    axis.text.x = element_text(size = 15,face="bold"),
    axis.title.y = element_text(size = 20,face="bold"),
    axis.text.y = element_blank())


dev.off()