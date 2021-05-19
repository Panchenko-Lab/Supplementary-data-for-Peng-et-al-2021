library("ggpubr")
library(matrixStats)
df1<-read.table("./amber/sasa_4500ns/dsasa_mean.csv", skip = 1, header = FALSE, col.names = "sasa")
df2<-read.table("./amber/DNA_contacts_4500ns/contact_mean_half.csv", skip = 1,
                header = FALSE, col.names = "contacts")

my_data <- data.frame(sasa = df1, contacts = df2)
p1 <-ggscatter(my_data, x = "sasa", y = "contacts", 
               #          add = "reg.line", conf.int = FALSE, 
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "Mean SASA change", 
               ylab = "Mean contacts number") + 
#  geom_abline(intercept = 0, slope = 1,linetype = "dashed") #+ xlim(0,25) + ylim(0,25) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
ggsave(plot = p1,dpi=300,width = 3.2,height = 2.8,units = c("in"),
       filename ="./amber/sasa_4500ns/SASA_correlation.png" )