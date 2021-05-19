#R-script to diplay DNA hbonds
library(matrixStats)

df11<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_1aoi_h2a_LYS.dat",sep=""),header=TRUE)
df12<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_1aoi_h2b_LYS.dat",sep=""),header=TRUE)
df13<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_1aoi_h3_LYS.dat",sep=""),header=TRUE)
df14<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_1aoi_h4_LYS.dat",sep=""),header=TRUE)

df21<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_1eqz_h2a_LYS.dat",sep=""),header=TRUE)
df22<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_1eqz_h2b_LYS.dat",sep=""),header=TRUE)
df23<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_1eqz_h3_LYS.dat",sep=""),header=TRUE)
df24<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_1eqz_h4_LYS.dat",sep=""),header=TRUE)

df31<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_ext2_h2a_LYS.dat",sep=""),header=TRUE)
df32<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_ext2_h2b_LYS.dat",sep=""),header=TRUE)
df33<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_ext2_h3_LYS.dat",sep=""),header=TRUE)
df34<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_ext2_h4_LYS.dat",sep=""),header=TRUE)

df41<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_sym_h2a_LYS.dat",sep=""),header=TRUE)
df42<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_sym_h2b_LYS.dat",sep=""),header=TRUE)
df43<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_sym_h3_LYS.dat",sep=""),header=TRUE)
df44<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_sym_h4_LYS.dat",sep=""),header=TRUE)

df51<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_shuxiang_ext_sym_amber_run1_h2a_LYS.dat",sep=""),header=TRUE)
df52<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_shuxiang_ext_sym_amber_run1_h2b_LYS.dat",sep=""),header=TRUE)
df53<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_shuxiang_ext_sym_amber_run1_h3_LYS.dat",sep=""),header=TRUE)
df54<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_shuxiang_ext_sym_amber_run1_h4_LYS.dat",sep=""),header=TRUE)

df61<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_shuxiang_ext_sym_amber_run2_h2a_LYS.dat",sep=""),header=TRUE)
df62<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_shuxiang_ext_sym_amber_run2_h2b_LYS.dat",sep=""),header=TRUE)
df63<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_shuxiang_ext_sym_amber_run2_h3_LYS.dat",sep=""),header=TRUE)
df64<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_shuxiang_ext_sym_amber_run2_h4_LYS.dat",sep=""),header=TRUE)

df11<-as.matrix(df11[,-1])
df12<-as.matrix(df12[,-1])
df13<-as.matrix(df13[,-1])
df14<-as.matrix(df14[,-1])

df21<-as.matrix(df21[,-1])
df22<-as.matrix(df22[,-1])
df23<-as.matrix(df23[,-1])
df24<-as.matrix(df24[,-1])

df31<-as.matrix(df31[,-1])
df32<-as.matrix(df32[,-1])
df33<-as.matrix(df33[,-1])
df34<-as.matrix(df34[,-1])

df41<-as.matrix(df41[,-1])
df42<-as.matrix(df42[,-1])
df43<-as.matrix(df43[,-1])
df44<-as.matrix(df44[,-1])

df51<-as.matrix(df51[,-1])
df52<-as.matrix(df52[,-1])
df53<-as.matrix(df53[,-1])
df54<-as.matrix(df54[,-1])

df61<-as.matrix(df61[,-1])
df62<-as.matrix(df62[,-1])
df63<-as.matrix(df63[,-1])
df64<-as.matrix(df64[,-1])

df1<-rbind(df11,df21,df31,df41,df51,df61)
df2<-rbind(df12,df22,df32,df42,df52,df62)
df3<-rbind(df13,df23,df33,df43,df53,df63)
df4<-rbind(df14,df24,df34,df44,df54,df64)

df11 <- df11[seq(1, nrow(df11), 5), ] ## 1ns interval of frames
df12 <- df12[seq(1, nrow(df12), 5), ]
df13 <- df13[seq(1, nrow(df13), 5), ]
df14 <- df14[seq(1, nrow(df14), 5), ]

df21 <- df21[seq(1, nrow(df21), 5), ] ## 1ns interval of frames
df22 <- df22[seq(1, nrow(df22), 5), ]
df23 <- df23[seq(1, nrow(df23), 5), ]
df24 <- df24[seq(1, nrow(df24), 5), ]

df31 <- df31[seq(1, nrow(df31), 5), ] ## 1ns interval of frames
df32 <- df32[seq(1, nrow(df32), 5), ]
df33 <- df33[seq(1, nrow(df33), 5), ]
df34 <- df34[seq(1, nrow(df34), 5), ]

df41 <- df41[seq(1, nrow(df41), 5), ] ## 1ns interval of frames
df42 <- df42[seq(1, nrow(df42), 5), ]
df43 <- df43[seq(1, nrow(df43), 5), ]
df44 <- df44[seq(1, nrow(df44), 5), ]

df51 <- df51[seq(1, nrow(df51), 5), ] ## 1ns interval of frames
df52 <- df52[seq(1, nrow(df52), 5), ]
df53 <- df53[seq(1, nrow(df53), 5), ]
df54 <- df54[seq(1, nrow(df54), 5), ]

df61 <- df61[seq(1, nrow(df61), 5), ] ## 1ns interval of frames
df62 <- df62[seq(1, nrow(df62), 5), ]
df63 <- df63[seq(1, nrow(df63), 5), ]
df64 <- df64[seq(1, nrow(df64), 5), ]

mean1<-colMeans(df1)
mean2<-colMeans(df2)
mean3<-colMeans(df3)
mean4<-colMeans(df4)


## Means for each run
#H2A
mean111 <- colMeans(df11[1:4300,])
mean112 <- colMeans(df11[4300:(4300+750),])
mean113 <- colMeans(df11[(4300+750):(4300+750*2),])
mean114 <- colMeans(df11[(4300+750*2):(4300+750*3),])
mean115 <- colMeans(df11[(4300+750*3):(4300+750*4),])

mean121 <- colMeans(df21[1:4800,])
mean122 <- colMeans(df21[4800:(4800+750),])
mean123 <- colMeans(df21[(4800+750):(4800+750*2),])
mean124 <- colMeans(df21[(4800+750*2):(4800+750*3),])
mean125 <- colMeans(df21[(4800+750*3):(4800+750*4),])

mean131 <- colMeans(df31[1:4300,])
mean132 <- colMeans(df31[4300:(4300+750),])
mean133 <- colMeans(df31[(4300+750):(4300+750*2),])
mean134 <- colMeans(df31[(4300+750*2):(4300+750*3),])
mean135 <- colMeans(df31[(4300+750*3):(4300+750*4),])

mean141 <- colMeans(df41[1:3800,])
mean142 <- colMeans(df41[3800:(3800+750),])
mean143 <- colMeans(df41[(3800+750):(3800+750*2),])
mean144 <- colMeans(df41[(3800+750*2):(3800+750*3),])
mean145 <- colMeans(df41[(3800+750*3):(3800+750*4),])

mean151 <- colMeans(df51[1:4800,])
mean161 <- colMeans(df61[1:4800,])

#H2B
mean211 <- colMeans(df12[1:4300,])
mean212 <- colMeans(df12[4300:(4300+750),])
mean213 <- colMeans(df12[(4300+750):(4300+750*2),])
mean214 <- colMeans(df12[(4300+750*2):(4300+750*3),])
mean215 <- colMeans(df12[(4300+750*3):(4300+750*4),])

mean221 <- colMeans(df22[1:4800,])
mean222 <- colMeans(df22[4800:(4800+750),])
mean223 <- colMeans(df22[(4800+750):(4800+750*2),])
mean224 <- colMeans(df22[(4800+750*2):(4800+750*3),])
mean225 <- colMeans(df22[(4800+750*3):(4800+750*4),])

mean231 <- colMeans(df32[1:4300,])
mean232 <- colMeans(df32[4300:(4300+750),])
mean233 <- colMeans(df32[(4300+750):(4300+750*2),])
mean234 <- colMeans(df32[(4300+750*2):(4300+750*3),])
mean235 <- colMeans(df32[(4300+750*3):(4300+750*4),])

mean241 <- colMeans(df42[1:3800,])
mean242 <- colMeans(df42[3800:(3800+750),])
mean243 <- colMeans(df42[(3800+750):(3800+750*2),])
mean244 <- colMeans(df42[(3800+750*2):(3800+750*3),])
mean245 <- colMeans(df42[(3800+750*3):(3800+750*4),])

mean251 <- colMeans(df52[1:4800,])
mean261 <- colMeans(df62[1:4800,])

#H3
mean311 <- colMeans(df13[1:4300,])
mean312 <- colMeans(df13[4300:(4300+750),])
mean313 <- colMeans(df13[(4300+750):(4300+750*2),])
mean314 <- colMeans(df13[(4300+750*2):(4300+750*3),])
mean315 <- colMeans(df13[(4300+750*3):(4300+750*4),])

mean321 <- colMeans(df23[1:4800,])
mean322 <- colMeans(df23[4800:(4800+750),])
mean323 <- colMeans(df23[(4800+750):(4800+750*2),])
mean324 <- colMeans(df23[(4800+750*2):(4800+750*3),])
mean325 <- colMeans(df23[(4800+750*3):(4800+750*4),])

mean331 <- colMeans(df33[1:4300,])
mean332 <- colMeans(df33[4300:(4300+750),])
mean333 <- colMeans(df33[(4300+750):(4300+750*2),])
mean334 <- colMeans(df33[(4300+750*2):(4300+750*3),])
mean335 <- colMeans(df33[(4300+750*3):(4300+750*4),])

mean341 <- colMeans(df43[1:3800,])
mean342 <- colMeans(df43[3800:(3800+750),])
mean343 <- colMeans(df43[(3800+750):(3800+750*2),])
mean344 <- colMeans(df43[(3800+750*2):(3800+750*3),])
mean345 <- colMeans(df43[(3800+750*3):(3800+750*4),])

mean351 <- colMeans(df53[1:4800,])
mean361 <- colMeans(df63[1:4800,])

#H4
mean411 <- colMeans(df14[1:4300,])
mean412 <- colMeans(df14[4300:(4300+750),])
mean413 <- colMeans(df14[(4300+750):(4300+750*2),])
mean414 <- colMeans(df14[(4300+750*2):(4300+750*3),])
mean415 <- colMeans(df14[(4300+750*3):(4300+750*4),])

mean421 <- colMeans(df24[1:4800,])
mean422 <- colMeans(df24[4800:(4800+750),])
mean423 <- colMeans(df24[(4800+750):(4800+750*2),])
mean424 <- colMeans(df24[(4800+750*2):(4800+750*3),])
mean425 <- colMeans(df24[(4800+750*3):(4800+750*4),])

mean431 <- colMeans(df34[1:4300,])
mean432 <- colMeans(df34[4300:(4300+750),])
mean433 <- colMeans(df34[(4300+750):(4300+750*2),])
mean434 <- colMeans(df34[(4300+750*2):(4300+750*3),])
mean435 <- colMeans(df34[(4300+750*3):(4300+750*4),])

mean441 <- colMeans(df44[1:3800,])
mean442 <- colMeans(df44[3800:(3800+750),])
mean443 <- colMeans(df44[(3800+750):(3800+750*2),])
mean444 <- colMeans(df44[(3800+750*2):(3800+750*3),])
mean445 <- colMeans(df44[(3800+750*3):(3800+750*4),])

mean451 <- colMeans(df54[1:4800,])
mean461 <- colMeans(df64[1:4800,])


sderr1 <- colSds(rbind(mean111,mean112,mean113,mean114,mean115,
                       mean121,mean122,mean123,mean124,mean125,
                       mean131,mean132,mean133,mean134,mean135,
                       mean141,mean142,mean143,mean144,mean145,
                       mean151,mean161))/sqrt(22)

sderr2 <- colSds(rbind(mean211,mean212,mean213,mean214,mean215,
                       mean221,mean222,mean223,mean224,mean225,
                       mean231,mean232,mean233,mean234,mean235,
                       mean241,mean242,mean243,mean244,mean245,
                       mean251,mean261))/sqrt(22)

sderr3 <- colSds(rbind(mean311,mean312,mean313,mean314,mean315,
                       mean321,mean322,mean323,mean324,mean325,
                       mean331,mean332,mean333,mean334,mean335,
                       mean341,mean342,mean343,mean344,mean345,
                       mean351,mean361))/sqrt(22)

sderr4 <- colSds(rbind(mean411,mean412,mean413,mean414,mean415,
                       mean421,mean422,mean423,mean424,mean425,
                       mean431,mean432,mean433,mean434,mean435,
                       mean441,mean442,mean443,mean444,mean445,
                       mean451,mean461))/sqrt(22)


contact_mean=mean1 + mean2 + mean3 +mean4
write.csv(contact_mean,"./amber/DNA_contacts_4500ns/contact_mean_LYS.csv", row.names = FALSE)