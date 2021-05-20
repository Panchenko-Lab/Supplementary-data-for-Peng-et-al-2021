library(ggplot2)
library(matrixStats)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(forcats)

h2a_residence_time <-read.table("h2a_residue_residence_time.csv",sep = ",", header = FALSE)
h2b_residence_time<-read.table("h2b_residue_residence_time.csv",sep = ",", header = FALSE)
h3_residence_time<-read.table("h3_residue_residence_time.csv",sep = ",", header = FALSE)
h4_residence_time<-read.table("h4_residue_residence_time.csv",sep = ",", header = FALSE)

be11<-read.table("MMPBSA_ddG_modelA_run1.csv",sep = ",", skip = 4, header=FALSE)
be12<-read.table("MMPBSA_ddG_modelA_run2.csv",sep = ",", skip = 4, header=FALSE)
be13<-read.table("MMPBSA_ddG_modelA_run3.csv",sep = ",", skip = 4, header=FALSE)
be14<-read.table("MMPBSA_ddG_modelA_run4.csv",sep = ",", skip = 4, header=FALSE)
be15<-read.table("MMPBSA_ddG_modelA_run5.csv",sep = ",", skip = 4, header=FALSE)

be21<-read.table("MMPBSA_ddG_modelB_run1.csv",sep = ",", skip = 4, header=FALSE)
be22<-read.table("MMPBSA_ddG_modelB_run2.csv",sep = ",", skip = 4, header=FALSE)
be23<-read.table("MMPBSA_ddG_modelB_run3.csv",sep = ",", skip = 4, header=FALSE)
be24<-read.table("MMPBSA_ddG_modelB_run4.csv",sep = ",", skip = 4, header=FALSE)
be25<-read.table("MMPBSA_ddG_modelB_run5.csv",sep = ",", skip = 4, header=FALSE)

be31<-read.table("MMPBSA_ddG_modelC_run1.csv",sep = ",", skip = 4, header=FALSE)
be32<-read.table("MMPBSA_ddG_modelC_run2.csv",sep = ",", skip = 4, header=FALSE)
be33<-read.table("MMPBSA_ddG_modelC_run3.csv",sep = ",", skip = 4, header=FALSE)
be34<-read.table("MMPBSA_ddG_modelC_run4.csv",sep = ",", skip = 4, header=FALSE)
be35<-read.table("MMPBSA_ddG_modelC_run5.csv",sep = ",", skip = 4, header=FALSE)

be41<-read.table("MMPBSA_ddG_modelD_run1.csv",sep = ",", skip = 4, header=FALSE)
be42<-read.table("MMPBSA_ddG_modelD_run2.csv",sep = ",", skip = 4, header=FALSE)
be43<-read.table("MMPBSA_ddG_modelD_run3.csv",sep = ",", skip = 4, header=FALSE)
be44<-read.table("MMPBSA_ddG_modelD_run4.csv",sep = ",", skip = 4, header=FALSE)
be45<-read.table("MMPBSA_ddG_modelD_run5.csv",sep = ",", skip = 4, header=FALSE)

be51<-read.table("MMPBSA_ddG_ModelD_gromacs_run1.csv",sep = ",", skip = 4, header=FALSE)
be52<-read.table("MMPBSA_ddG_ModelD_gromacs_run2.csv",sep = ",", skip = 4, header=FALSE)


## Mean binding energy and standard error per tail type
mean_be_h3 <- (be11$V18[984:1019] + be11$V18[625:660] + be12$V18[984:1019] + be12$V18[625:660] + be13$V18[984:1019] + be13$V18[625:660] +
            be14$V18[984:1019] + be14$V18[625:660] + be15$V18[984:1019] + be15$V18[625:660] +
          be21$V18[375:410] + be21$V18[1092:1127] + be22$V18[375:410] + be22$V18[1092:1127] + be23$V18[375:410] + be23$V18[1092:1127] +
            be24$V18[375:410] + be24$V18[1092:1127] + be25$V18[375:410] + be25$V18[1092:1127] +
          be31$V18[1112:1147] + be31$V18[605:640] + be32$V18[1112:1147] + be32$V18[605:640] + be33$V18[1112:1147] + be33$V18[605:640] +
            be34$V18[1112:1147] + be34$V18[605:640] + be35$V18[1112:1147] + be35$V18[605:640] +
          be41$V18[257:292] + be41$V18[738:773] + be42$V18[257:292] + be42$V18[738:773] +  be43$V18[257:292] + be43$V18[738:773] +
            be44$V18[257:292] + be44$V18[738:773]  + be45$V18[257:292] + be45$V18[738:773] + 
          be51$V18[605:640] + be51$V18[1112:1147]  + be52$V18[605:640] + be52$V18[1112:1147] 
          ) /44

se_h3 <- colSds(rbind(be11$V18[984:1019] , be11$V18[625:660] , be12$V18[984:1019] , be12$V18[625:660] , be13$V18[984:1019] , be13$V18[625:660],
                        be14$V18[984:1019] , be14$V18[625:660] , be15$V18[984:1019] , be15$V18[625:660] ,
                        be21$V18[375:410] , be21$V18[1092:1127] , be22$V18[375:410] , be22$V18[1092:1127] , be23$V18[375:410] , be23$V18[1092:1127],
                        be24$V18[375:410] , be24$V18[1092:1127] , be25$V18[375:410] , be25$V18[1092:1127] ,
                        be31$V18[1112:1147] , be31$V18[605:640] , be32$V18[1112:1147] , be32$V18[605:640] , be33$V18[1112:1147] ,be33$V18[605:640],
                        be34$V18[1112:1147] , be34$V18[605:640] , be35$V18[1112:1147] , be35$V18[605:640] ,
                        be41$V18[257:292] , be41$V18[738:773] , be42$V18[257:292] , be42$V18[738:773] ,  be43$V18[257:292], be43$V18[738:773],
                        be44$V18[257:292] , be44$V18[738:773] , be45$V18[257:292] , be45$V18[738:773] , 
                        be51$V18[605:640] , be51$V18[1112:1147] , be52$V18[605:640] , be52$V18[1112:1147] ))/sqrt(44)

mean_be_h4 <- (be11$V18[1119:1138] + be11$V18[760:779] + be12$V18[1119:1138] + be12$V18[760:779] +  be13$V18[1119:1138] + be13$V18[760:779] + 
            be14$V18[1119:1138] + be14$V18[760:779] +  be15$V18[1119:1138] + be15$V18[760:779] + 
          be21$V18[510:529] + be21$V18[740:759] + be22$V18[510:529] + be22$V18[740:759] + be23$V18[510:529] + be23$V18[740:759] +
            be24$V18[510:529] + be24$V18[740:759] + be25$V18[510:529] + be25$V18[740:759] +
          be31$V18[503:522] + be31$V18[1247:1266] + be32$V18[503:522] + be32$V18[1247:1266] + be33$V18[503:522] + be33$V18[1247:1266] +
            be34$V18[503:522] + be34$V18[1247:1266] + be35$V18[503:522] + be35$V18[1247:1266] +
          be41$V18[514:533] + be41$V18[873:892] + be42$V18[514:533] + be42$V18[873:892] + be43$V18[514:533] + be43$V18[873:892] +
            be44$V18[514:533] + be44$V18[873:892] + be45$V18[514:533] + be45$V18[873:892] + 
          be51$V18[503:522] + be51$V18[1247:1266]  + be52$V18[503:522] + be52$V18[1247:1266]  
          ) /44

se_h4 <- colSds(rbind(be11$V18[1119:1138] , be11$V18[760:779] , be12$V118[1119:1138] , be12$V18[760:779] ,  be13$V18[1119:1138] , be13$V18[760:779], 
                        be14$V18[1119:1138] , be14$V18[760:779] ,  be15$V18[1119:1138] , be15$V18[760:779] , 
                        be21$V18[510:529] , be21$V18[740:759] , be22$V18[510:529] , be22$V18[740:759] , be23$V18[510:529] , be23$V18[740:759],
                        be24$V18[510:529] , be24$V18[740:759] , be25$V18[510:529] , be25$V18[740:759] ,
                        be31$V18[503:522] , be31$V18[1247:1266] , be32$V18[503:522] , be32$V18[1247:1266] , be33$V18[503:522] , be33$V18[1247:1266],
                        be34$V18[503:522] , be34$V18[1247:1266] , be35$V18[503:522] , be35$V18[1247:1266] ,
                        be41$V18[514:533] , be41$V18[873:892] , be42$V18[514:533] , be42$V18[873:892] , be43$V18[514:533] , be43$V18[873:892],
                        be44$V18[514:533] , be44$V18[873:892] , be45$V18[514:533] , be45$V18[873:892] , 
                        be51$V18[503:522] , be51$V18[1247:1266]  , be52$V18[503:522] , be52$V18[1247:1266] ))/sqrt(44)


mean_be_h2a_N <- (be11$V18[375:387] + be11$V18[1221:1233] + be12$V18[375:387] + be12$V18[1221:1233] +  be13$V18[375:387] + be13$V18[1221:1233] +  
               be14$V18[375:387] + be14$V18[1221:1233] +  be15$V18[375:387] + be15$V18[1221:1233] + 
             be21$V18[612:624] + be21$V18[842:854] +  be22$V18[612:624] + be22$V18[842:854] +  be23$V18[612:624] + be23$V18[842:854] +
               be24$V18[612:624] + be24$V18[842:854] +  be25$V18[612:624] + be25$V18[842:854] +
             be31$V18[984:996] + be31$V18[375:387] + be32$V18[984:996] + be32$V18[375:387] + be33$V18[984:996] + be33$V18[375:387] +
               be34$V18[984:996] + be34$V18[375:387] + be35$V18[984:996] + be35$V18[375:387] +
             be41$V18[129:141] + be41$V18[1:13] + be42$V18[129:141] + be42$V18[1:13] +  be43$V18[129:141] + be43$V18[1:13] +  
               be44$V18[129:141] + be44$V18[1:13] +  be45$V18[129:141] + be45$V18[1:13] +
             be51$V18[375:387] + be51$V18[984:996] + be52$V18[375:387] + be52$V18[984:996]
             )/44

se_h2a_N <- colSds(rbind(be11$V18[375:387] , be11$V18[1221:1233] , be12$V18[375:387] , be12$V18[1221:1233] ,  be13$V18[375:387] , be13$V18[1221:1233] ,  
                           be14$V18[375:387] , be14$V18[1221:1233] ,  be15$V18[375:387] , be15$V18[1221:1233] , 
                           be21$V18[612:624] , be21$V18[842:854] ,  be22$V18[612:624] , be22$V18[842:854] ,  be23$V18[612:624] , be23$V18[842:854] ,
                           be24$V18[612:624] , be24$V18[842:854] ,  be25$V18[612:624] , be25$V18[842:854] ,
                           be31$V18[984:996] , be31$V18[375:387] , be32$V18[984:996] , be32$V18[375:387] , be33$V18[984:996] , be33$V18[375:387] ,
                           be34$V18[984:996] , be34$V18[375:387] , be35$V18[984:996] , be35$V18[375:387] ,
                           be41$V18[129:141] , be41$V18[1:13] , be42$V18[129:141] , be42$V18[1:13] ,  be43$V18[129:141] , be43$V18[1:13] +  
                           be44$V18[129:141] , be44$V18[1:13] ,  be45$V18[129:141] , be45$V18[1:13] ,
                           be51$V18[375:387] , be51$V18[984:996] , be52$V18[375:387] , be52$V18[984:996]))/sqrt(44)


mean_be_h2a_C <- (be11$V18[493:502] + be11$V18[1339:1348] +  be12$V18[493:502] + be12$V18[1339:1348] + be13$V18[493:502] + be13$V18[1339:1348] +
               be14$V18[493:502] + be14$V18[1339:1348] + be15$V18[493:502] + be15$V18[1339:1348] + 
             be21$V18[730:739] + be21$V18[960:969] +  be22$V18[730:739] + be22$V18[960:969] +   be23$V18[730:739] + be23$V18[960:969] + 
               be24$V18[730:739] + be24$V18[960:969] +  be25$V18[730:739] + be25$V18[960:969] + 
             be31$V18[1102:1111] + be31$V18[493:502] + be32$V18[1102:1111] + be32$V18[493:502] + be33$V18[1102:1111] + be33$V18[493:502] +
               be34$V18[1102:1111] + be34$V18[493:502] + be35$V18[1102:1111] + be35$V18[493:502] +
             be41$V18[247:256] + be41$V18[119:128] + be42$V18[247:256] + be42$V18[119:128] +  be43$V18[247:256] + be43$V18[119:128] + 
               be44$V18[247:256] + be44$V18[119:128] +  be45$V18[247:256] + be45$V18[119:128] +
             be51$V18[493:502] + be51$V18[1102:1111] + be52$V18[493:502] + be52$V18[1102:1111]
             ) /44

se_h2a_C <- colSds(rbind(be11$V18[493:502] , be11$V18[1339:1348] ,  be12$V18[493:502] , be12$V18[1339:1348] , be13$V18[493:502] , be13$V18[1339:1348],
                           be14$V18[493:502] , be14$V18[1339:1348] , be15$V18[493:502] , be15$V18[1339:1348] , 
                           be21$V18[730:739] , be21$V18[960:969] ,  be22$V18[730:739] , be22$V18[960:969] ,   be23$V18[730:739] , be23$V18[960:969],
                           be24$V18[730:739] , be24$V18[960:969] ,  be25$V18[730:739] , be25$V18[960:969] , 
                           be31$V18[1102:1111] , be31$V18[493:502] , be32$V18[1102:1111] , be32$V18[493:502] , be33$V18[1102:1111] , be33$V18[493:502],
                           be34$V18[1102:1111] , be34$V18[493:502] , be35$V18[1102:1111] , be35$V18[493:502] ,
                           be41$V18[247:256] , be41$V18[119:128] , be42$V18[247:256] , be42$V18[119:128] ,  be43$V18[247:256] , be43$V18[119:128],
                           be44$V18[247:256] , be44$V18[119:128] ,  be45$V18[247:256] , be45$V18[119:128] ,
                           be51$V18[493:502] , be51$V18[1102:1111] , be52$V18[493:502] , be52$V18[1102:1111]))/sqrt(44)


mean_be_h2b <- (be11$V18[503:525] + be11$V18[862:884] + be12$V18[503:525] + be12$V18[862:884] +  be13$V18[503:525] + be13$V18[862:884] + 
             be14$V18[503:525] + be14$V18[862:884] + be15$V18[503:525] + be15$V18[862:884] + 
           be21$V18[1227:1249] + be21$V18[970:992] + be22$V18[1227:1249] + be22$V18[970:992] + be23$V18[1227:1249] + be23$V18[970:992] +
             be24$V18[1227:1249] + be24$V18[970:992] + be25$V18[1227:1249] + be25$V18[970:992] +
           be31$V18[862:884] + be31$V18[740:762] + be32$V18[862:884] + be32$V18[740:762] +  be33$V18[862:884] + be33$V18[740:762] + 
             be34$V18[862:884] + be34$V18[740:762] + be35$V18[862:884] + be35$V18[740:762] + 
           be41$V18[392:414] + be41$V18[616:638] + be42$V18[392:414] + be42$V18[616:638] + be43$V18[392:414] + be43$V18[616:638] +
             be44$V18[392:414] + be44$V18[616:638] + be45$V18[392:414] + be45$V18[616:638] +
           be51$V18[740:762] + be51$V18[862:884] + be52$V18[740:762] + be52$V18[862:884]
           ) /44


se_h2b <- colSds(rbind(be11$V18[503:525] , be11$V18[862:884] , be12$V18[503:525] , be12$V18[862:884] ,  be13$V18[503:525] , be13$V18[862:884] ,
                         be14$V18[503:525] , be14$V18[862:884] , be15$V18[503:525] , be15$V18[862:884] , 
                         be21$V18[1227:1249] , be21$V18[970:992] , be22$V18[1227:1249] , be22$V18[970:992] , be23$V18[1227:1249] , be23$V18[970:992],
                         be24$V18[1227:1249] , be24$V18[970:992] , be25$V18[1227:1249] , be25$V18[970:992],
                         be31$V18[862:884] , be31$V18[740:762] , be32$V18[862:884] , be32$V18[740:762] ,  be33$V18[862:884] , be33$V18[740:762],
                         be34$V18[862:884] , be34$V18[740:762] , be35$V18[862:884] , be35$V18[740:762] , 
                         be41$V18[392:414] , be41$V18[616:638] , be42$V18[392:414] , be42$V18[616:638] , be43$V18[392:414] , be43$V18[616:638],
                         be44$V18[392:414] , be44$V18[616:638] , be45$V18[392:414] , be45$V18[616:638] ,
                         be51$V18[740:762] , be51$V18[862:884] , be52$V18[740:762] , be52$V18[862:884]))/sqrt(44)

## Residue name
h2a_Ntail = c("S1","G2","R3","G4","K5","Q6","G7","G8","K9","T10","R11","A12","K13")

h2a_Ctail = c("K119","T120","E121","S122","S123","K124","S125","K126","S127","K128" )

h2b_tail = c("A1","K2","S3","A4","P5","A6","P7","K8","K9","G10","S11","K12","K13","A14","V15",
             "T16","K17","T18","Q19","K20","K21","D22","G23")

h3_tail <- c("A1","R2","T3","K4","Q5","T6","A7","R8","K9","S10","T11","G12","G13",
             "K14","A15","P16","R17","K18","Q19","L20","A21","T22","K23","A24","A25","R26",
             "K27","S28","A29","P30","A31","T32","G33","G34","V35","K36")

h4_tail = c("S1","G2","R3","G4","K5","G6","G7","K8","G9","L10","G11","K12","G13","G14","A15","K16","R17","H18","R19","K20")

RT <- data.frame( number = h3_residence_time[,1], contacts = h3_residence_time[,2],
                  tail = h3_tail,RT_se=h3_residence_time[,3])
BE <- data.frame(label = h3_tail ,be = mean_be_h3, be_se = se_h3)

residue_number =seq(1,36)
ggplot() +
  geom_line(data=RT,aes(x=number, y=RT$contacts),stat="identity", colour="#0000FF") +
      labs(x = "H2B Tail Residues", y ="Mean Residence Time (ns)") +
  geom_errorbar(data=RT, aes(x=number, y=RT$contacts,ymin=RT$contacts-RT$RT_se/2, ymax=RT$contacts+RT$RT_se/2), width=.2, colour="#0000FF",
                position=position_dodge(0.9))+
  geom_hline(yintercept=0, linetype="dashed") +
  geom_line(data=BE,aes(x=number, y=BE$be),stat="identity", colour="#FF0000")+ 
  scale_x_continuous(breaks=number,labels = BE$label) +
  geom_errorbar(data=BE, aes(x=number, y=BE$be,ymin=BE$be-BE$be_se/2, ymax=BE$be+BE$be_se/2), width=.2, colour="#FF0000",
                position=position_dodge(0.9))+ theme_light() +  
  scale_y_continuous(breaks=c(-20,-15,-10,-5,0,5,10,15,20,25),labels = c(-20,-15,-10,-5,0,5,10,15,20,25),
                     sec.axis = sec_axis(~. ,breaks=c(-20,-15,-10,-5,0,5,10,15,20,25),labels = c(-20,-15,-10,-5,0,5,10,15,20,25), name = "Binding Free Energy(kcal/mol)"))  + #limits=c(-4,3) 
  theme(axis.title.y.right = element_text(colour = "#FF0000", face="bold",size = 13),
        panel.background = element_rect(fill = NA,colour = "black"),panel.grid = element_blank(),
        axis.title.y.left = element_text(colour = "#0000FF",face="bold",size = 13),
        axis.text.y.right = element_text(colour = "#FF0000",face="bold",size = 16),
        axis.text.y.left = element_text(colour = "#0000FF",face="bold",size = 16),
        axis.title.x = element_text( face="bold",size = 14,colour="black"),
        axis.text.x = element_text(face="bold",size = 14,colour="black",angle = 90, vjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.margin = unit(c(1,0.2,0.2,0.2), "cm"))

ggsave("Residence_time_BE_H3.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 8, height = 3, units = c("in"),
       dpi = 300, limitsize = TRUE)

