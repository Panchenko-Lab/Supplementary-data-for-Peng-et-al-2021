library(ggplot2)
library(matrixStats)

tail1 = "h2a_C"
tail2 = "h2a_G"
df1<-read.table(paste("../tail_contacts/tail_contacts_all_1aoi_",tail1,".dat",sep=""),header=TRUE)
df2<-read.table(paste("../tail_contacts/tail_contacts_all_1eqz_",tail1,".dat",sep=""),header=TRUE)
df3<-read.table(paste("../tail_contacts/tail_contacts_all_ext2_",tail1,".dat",sep=""),header=TRUE)
df4<-read.table(paste("../tail_contacts/tail_contacts_all_sym_",tail1,".dat",sep=""),header=TRUE)
df5<-read.table(paste("../tail_contacts/tail_contacts_all_1aoi_",tail2,".dat",sep=""),header=TRUE)
df6<-read.table(paste("../tail_contacts/tail_contacts_all_1eqz_",tail2,".dat",sep=""),header=TRUE)
df7<-read.table(paste("../tail_contacts/tail_contacts_all_ext2_",tail2,".dat",sep=""),header=TRUE)
df8<-read.table(paste("../tail_contacts/tail_contacts_all_sym_",tail2,".dat",sep=""),header=TRUE)

be1<-read.table("./MMPBSA_ddG_1aoi.csv",sep = ",", skip = 4, header=FALSE)
be2<-read.table("./MMPBSA_ddG_1eqz.csv",sep = ",", skip = 4, header=FALSE)
be3<-read.table("./MMPBSA_ddG_ext2.csv",sep = ",", skip = 4, header=FALSE)
be4<-read.table("./MMPBSA_ddG_sym.csv",sep = ",", skip = 4, header=FALSE)

be_h3 <- (be1$V18[984:1019] + be1$V18[625:660] + be2$V18[375:410] + be2$V18[1092:1127]
    + be3$V18[1112:1147] + be3$V18[605:640] + be4$V18[257:292] + be4$V18[738:773]) /6
sd_h3 <- (be1$V19[984:1019] + be1$V19[625:660] + be2$V19[375:410] + be2$V19[1092:1127]
          + be3$V19[1112:1147] + be3$V19[605:640] + be4$V19[257:292] + be4$V19[738:773]) /6
se_h3 <- (be1$V20[984:1019] + be1$V20[625:660] + be2$V20[375:410] + be2$V20[1092:1127]
          + be3$V20[1112:1147] + be3$V20[605:640] + be4$V20[257:292] + be4$V20[738:773]) /6

be_h4 <- (be1$V18[1119:1133] + be1$V18[760:774] + be2$V18[510:524] + be2$V18[740:754]
          + be3$V18[503:517] + be3$V18[1247:1261] + be4$V18[514:528] + be4$V18[873:887]) /6
sd_h4 <- (be1$V19[1119:1133] + be1$V19[760:774] + be2$V19[510:524] + be2$V19[740:754]
          + be3$V19[503:517] + be3$V19[1247:1261] + be4$V19[514:528] + be4$V19[873:887]) /6
se_h4 <- (be1$V20[1119:1133] + be1$V20[760:774] + be2$V20[510:524] + be2$V20[740:754]
          + be3$V20[503:517] + be3$V20[1247:1261] + be4$V20[514:528] + be4$V20[873:887]) /6

be_h2a_N <- (be1$V18[375:385] + be1$V18[1221:1231] + be2$V18[612:622] + be2$V18[842:852]
          + be3$V18[984:994] + be3$V18[375:385] + be4$V18[129:139] + be4$V18[1:11]) /6
sd_h2a_N <- (be1$V19[375:385] + be1$V19[1221:1231] + be2$V19[612:622] + be2$V19[842:852]
             + be3$V19[984:994] + be3$V19[375:385] + be4$V19[129:139] + be4$V19[1:11]) /6
se_h2a_N <- (be1$V20[375:385] + be1$V20[1221:1231] + be2$V20[612:622] + be2$V20[842:852]
             + be3$V20[984:994] + be3$V20[375:385] + be4$V20[129:139] + be4$V20[1:11]) /6

be_h2a_C <- (be1$V18[493:502] + be1$V18[1339:1348] + be2$V18[730:739] + be2$V18[960:969]
             + be3$V18[1102:1111] + be3$V18[493:502] + be4$V18[247:256] + be4$V18[119:128]) /6
sd_h2a_C <- (be1$V19[493:502] + be1$V19[1339:1348] + be2$V19[730:739] + be2$V19[960:969]
             + be3$V19[1102:1111] + be3$V19[493:502] + be4$V19[247:256] + be4$V19[119:128]) /6
se_h2a_C <- (be1$V20[493:502] + be1$V20[1339:1348] + be2$V20[730:739] + be2$V20[960:969]
             + be3$V20[1102:1111] + be3$V20[493:502] + be4$V20[247:256] + be4$V20[119:128]) /6


be_h2b <- (be1$V18[503:525] + be1$V18[862:884] + be2$V18[1227:1249] + be2$V18[970:992]
          + be3$V18[862:884] + be3$V18[740:762] + be4$V18[392:414] + be4$V18[616:638]) /6
sd_h2b <- (be1$V19[503:525] + be1$V19[862:884] + be2$V19[1227:1249] + be2$V19[970:992]
           + be3$V19[862:884] + be3$V19[740:762] + be4$V19[392:414] + be4$V19[616:638]) /6
se_h2b <- (be1$V20[503:525] + be1$V20[862:884] + be2$V20[1227:1249] + be2$V20[970:992]
           + be3$V20[862:884] + be3$V20[740:762] + be4$V20[392:414] + be4$V20[616:638]) /6

# find the mean contacts and standard error of mean
df <-rbind(df1,df2,df3,df4,df5,df6,df7,df8)
df <-as.matrix(df[,-1])

mean <-colMeans(df)
sd <- colSds(df)
se <- colSds(df)/sqrt(nrow(df))

h2a_Ntail = c("S1","G2","R3","G4","K5","Q6","G7","G8","K9","T10","R11")

h2a_Ctail = c("K119","T120","E121","S122","S123","K124","S125","K126","S127","K128" )

h2b_tail = c("A1","K2","S3","A4","P5","A6","P7","K8","K9","G10","S11","K12","K13","A14","V15",
                                         "T16","K17","T18","Q19","K20","K21","D22","G23")

h3_tail <- c("A1","R2","T3","K4","Q5","T6","A7","R8","K9","S10","T11","G12","G13",
             "K14","A15","P16","R17","K18","Q19","L20","A21","T22","K23","A24","A25","R26",
             "K27","S28","A29","P30","A31","T32","G33","G34","V35","K36")

h4_tail = c("S1","G2","R3","G4","K5","G6","G7","K8","G9","L10","G11","K12","G13","G14","A15")

df <- data.frame(tail = h2a_Ctail, contacts = mean[12:21], contacts_sd = sd[12:21], be = be_h2a_C, be_sd = sd_h2a_C)
#df <- data.frame(tail = h4_tail, contacts = mean, contacts_sd = sd, be = be_h4, be_sd = sd_h4)
ggplot(df)  + 
  geom_bar(aes(x=df$tail, y=df$contacts),stat="identity", fill="tan1", colour="sienna3")+ 
  geom_errorbar(aes(x=df$tail, ymin=df$contacts - df$contacts_sd/4, ymax=df$contacts + df$contacts_sd/4), width=.2,colour="sienna3",
                position=position_dodge(.9)) + 
  labs(x = "H2A C-terminal tail residues", y ="Mean number of contacts with DNA") + 
  geom_line(aes(x=seq(1,10), y=df$be),stat="identity", colour="blue")+ scale_x_discrete(limits = df$tail) +
  geom_errorbar(aes(x=seq(1,10), ymin=df$be-df$be_sd/2, ymax=df$be+df$be_sd/2), width=.2, colour="blue",
                position=position_dodge(0.9))+ theme_light() +   
  scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Binding free energy(kcal/mol)")) +
  theme(axis.title.y.right = element_text(colour = "blue", face="bold"),
        axis.title.y.left = element_text(colour = "sienna3",face="bold"),
        axis.text.y.right = element_text(colour = "blue",face="bold"),
        axis.text.y.left = element_text(colour = "sienna3",face="bold"),
        axis.title.x = element_text( face="bold"),
        axis.text.x = element_text(face="bold"))  + expand_limits(y=c(-15, 20))

ggsave("./contacts_BE_H2A_C.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 10, height = 6, units = c("in"),
       dpi = 300, limitsize = TRUE)

