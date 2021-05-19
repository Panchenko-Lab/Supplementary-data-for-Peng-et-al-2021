library("ggplot2")

my_data <- data.frame(Kd1 = log(c(0.002609355,0.147248739,0.032770606,0.005971343,0.003641417)), 
                      Kd2 = log(c(4.50219E-23,4.88211E-07,1.73313E-21,3.97414E-48,1.57616E-18)))

p1 <-ggscatter(my_data, x = "Kd1", y = "Kd2", 
                         add = "reg.line", conf.int = FALSE, 
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "Ln(Kd)_conformation", 
               ylab = "Ln(Kd)_MM/GBSA") + 
  geom_abline(intercept = 0, slope = 1,linetype = "dashed") + #xlim(0,20) + ylim(0,20) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
ggsave(plot = p1,dpi=300,width = 4,height = 4,units = c("in"),
       filename ="./Kd_correlation.png" )