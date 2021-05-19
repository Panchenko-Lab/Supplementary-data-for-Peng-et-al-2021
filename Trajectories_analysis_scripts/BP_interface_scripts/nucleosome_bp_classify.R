library(ggplot2)
library(scales)


blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=20, face="bold"))

aa_wt <- data.frame(
  type = c("Remodelling", "Transcription", "Histone Modification",
           "Chromosome Segregation", "Ubiquitination",
           "DNA Repair","Others"),
  count = c(44, 41, 39, 20, 7,7,16))
aa_wt$type <- factor(aa_wt$type,levels = c("Remodelling", "Transcription","Histone Modification",
                                           "Chromosome Segregation", "Ubiquitination",
                                           "DNA Repair","Others"))
pie1 <- ggplot(aa_wt, aes(x="", y=count, fill=type))+
  geom_bar(width = 5, stat = "identity",color='black',lwd=1.5) + coord_polar("y")  + blank_theme + 
  guides(fill=guide_legend(override.aes=list(colour=NA)))+ #remove line in legend 
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#FF00FF", "#00FF00","#FFFACD","#808000","#C0C0C0")) +
  theme(axis.text = element_text(face="bold", color="#993333", 
                                 size=12, angle=0),
        axis.title = element_text(face="bold", color="#993333", 
                                  size=12, angle=0),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold", color="black", 
                                   size=16, angle=0),
        legend.position = "top") + guides(fill=guide_legend(nrow=4,byrow=TRUE))
#  geom_text(aes(y = value, 
#                label = percent(value/100)), size=5) +


ggsave("nucleosme_bp_classify.png",plot = pie1,dpi = 300, height = 6, width = 10)
