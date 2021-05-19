#R-script to diplay DNA hbonds
library(matrixStats)
if(0) {
df11<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_sym_h2a.dat",sep=""),header=TRUE)
df12<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_sym_h2b.dat",sep=""),header=TRUE)
df13<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_sym_h3.dat",sep=""),header=TRUE)
df14<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_sym_h4.dat",sep=""),header=TRUE)

df21<-read.table(paste("./amber/PTMs_4500ns/DNA_contacts/dna_prot_contacts_ARG_ala_h2a.dat",sep=""),header=TRUE)
df22<-read.table(paste("./amber/PTMs_4500ns/DNA_contacts/dna_prot_contacts_ARG_ala_h2b.dat",sep=""),header=TRUE)
df23<-read.table(paste("./amber/PTMs_4500ns/DNA_contacts/dna_prot_contacts_ARG_ala_h3.dat",sep=""),header=TRUE)
df24<-read.table(paste("./amber/PTMs_4500ns/DNA_contacts/dna_prot_contacts_ARG_ala_h4.dat",sep=""),header=TRUE)

df31<-read.table(paste("./amber/PTMs_4500ns/DNA_contacts/dna_prot_contacts_LYS_Ac_h2a.dat",sep=""),header=TRUE)
df32<-read.table(paste("./amber/PTMs_4500ns/DNA_contacts/dna_prot_contacts_LYS_Ac_h2b.dat",sep=""),header=TRUE)
df33<-read.table(paste("./amber/PTMs_4500ns/DNA_contacts/dna_prot_contacts_LYS_Ac_h3.dat",sep=""),header=TRUE)
df34<-read.table(paste("./amber/PTMs_4500ns/DNA_contacts/dna_prot_contacts_LYS_Ac_h4.dat",sep=""),header=TRUE)

df41<-read.table(paste("./amber/PTMs_4500ns/DNA_contacts/dna_prot_contacts_LYS_Me_h2a.dat",sep=""),header=TRUE)
df42<-read.table(paste("./amber/PTMs_4500ns/DNA_contacts/dna_prot_contacts_LYS_Me_h2b.dat",sep=""),header=TRUE)
df43<-read.table(paste("./amber/PTMs_4500ns/DNA_contacts/dna_prot_contacts_LYS_Me_h3.dat",sep=""),header=TRUE)
df44<-read.table(paste("./amber/PTMs_4500ns/DNA_contacts/dna_prot_contacts_LYS_Me_h4.dat",sep=""),header=TRUE)

df51<-read.table(paste("./amber/PTMs_4500ns/DNA_contacts/dna_prot_contacts_SER_ph_h2a.dat",sep=""),header=TRUE)
df52<-read.table(paste("./amber/PTMs_4500ns/DNA_contacts/dna_prot_contacts_SER_ph_h2b.dat",sep=""),header=TRUE)
df53<-read.table(paste("./amber/PTMs_4500ns/DNA_contacts/dna_prot_contacts_SER_ph_h3.dat",sep=""),header=TRUE)
df54<-read.table(paste("./amber/PTMs_4500ns/DNA_contacts/dna_prot_contacts_SER_ph_h4.dat",sep=""),header=TRUE)

df61<-read.table(paste("./amber/PTMs_4500ns/DNA_contacts/dna_prot_contacts_select_h2a.dat",sep=""),header=TRUE)
df62<-read.table(paste("./amber/PTMs_4500ns/DNA_contacts/dna_prot_contacts_select_h2b.dat",sep=""),header=TRUE)
df63<-read.table(paste("./amber/PTMs_4500ns/DNA_contacts/dna_prot_contacts_select_h3.dat",sep=""),header=TRUE)
df64<-read.table(paste("./amber/PTMs_4500ns/DNA_contacts/dna_prot_contacts_select_h4.dat",sep=""),header=TRUE)


df11<-as.matrix(rbind(df11[1:7000,-1],df11[19000:34000,-1]))
df12<-as.matrix(rbind(df12[1:7000,-1],df12[19000:34000,-1]))
df13<-as.matrix(rbind(df13[1:7000,-1],df13[19000:34000,-1]))
df14<-as.matrix(rbind(df14[1:7000,-1],df14[19000:34000,-1]))

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

## 1ns interval of frames
df11 <- df11[seq(1, nrow(df11), 5), ] 
df12 <- df12[seq(1, nrow(df12), 5), ]
df13 <- df13[seq(1, nrow(df13), 5), ]
df14 <- df14[seq(1, nrow(df14), 5), ]

df21 <- df21[seq(1, nrow(df21), 5), ] 
df22 <- df22[seq(1, nrow(df22), 5), ]
df23 <- df23[seq(1, nrow(df23), 5), ]
df24 <- df24[seq(1, nrow(df24), 5), ]

df31 <- df31[seq(1, nrow(df31), 5), ] 
df32 <- df32[seq(1, nrow(df32), 5), ]
df33 <- df33[seq(1, nrow(df33), 5), ]
df34 <- df34[seq(1, nrow(df34), 5), ]

df41 <- df41[seq(1, nrow(df41), 5), ] 
df42 <- df42[seq(1, nrow(df42), 5), ]
df43 <- df43[seq(1, nrow(df43), 5), ]
df44 <- df44[seq(1, nrow(df44), 5), ]

df51 <- df51[seq(1, nrow(df51), 5), ] 
df52 <- df52[seq(1, nrow(df52), 5), ]
df53 <- df53[seq(1, nrow(df53), 5), ]
df54 <- df54[seq(1, nrow(df54), 5), ]

df61 <- df61[seq(1, nrow(df61), 5), ] 
df62 <- df62[seq(1, nrow(df62), 5), ]
df63 <- df63[seq(1, nrow(df63), 5), ]
df64 <- df64[seq(1, nrow(df64), 5), ]
}
### mean for WT tail

df1_half_wt=df11[,94]
df2_half_wt=df12[,94]
df3_half_wt=df13[,94]
df4_half_wt=df14[,94]

for (i in seq(1:93))
{
  df1_half_wt <- cbind(df1_half_wt,(df11[,94-i] + df11[,94+i])/2 ) 
  df2_half_wt <- cbind(df2_half_wt,(df12[,94-i] + df12[,94+i])/2 ) 
  df3_half_wt <- cbind(df3_half_wt,(df13[,94-i] + df13[,94+i])/2 ) 
  df4_half_wt <- cbind(df4_half_wt,(df14[,94-i] + df14[,94+i])/2 ) 
}
mean1_wt<-colMeans(df1_half_wt)
mean2_wt<-colMeans(df2_half_wt)
mean3_wt<-colMeans(df3_half_wt)
mean4_wt<-colMeans(df4_half_wt)

## mean for modified tail
df1_half=df61[,94]
df2_half=df62[,94]
df3_half=df63[,94]
df4_half=df64[,94]

for (i in seq(1:93))
{
  df1_half <- cbind(df1_half,(df61[,94-i] + df61[,94+i])/2 ) 
  df2_half <- cbind(df2_half,(df62[,94-i] + df62[,94+i])/2 ) 
  df3_half <- cbind(df3_half,(df63[,94-i] + df63[,94+i])/2 ) 
  df4_half <- cbind(df4_half,(df64[,94-i] + df64[,94+i])/2 ) 
}
mean1<-colMeans(df1_half)
mean2<-colMeans(df2_half)
mean3<-colMeans(df3_half)
mean4<-colMeans(df4_half)

## T test for changes in mean
h2a_N <- (mean1_wt-mean1)[30:60]
h2a_C <- c((mean1_wt-mean1)[0:30],(mean1_wt-mean1)[60:93])
h2b <- mean2_wt-mean2
h3 <- mean3_wt-mean3
h4 <- mean4_wt-mean4

t.test(h2a_N[h2a_N!=0], mu=0, alternative="great", conf.level= 0.99)
t.test(h2a_C[h2a_C!=0], mu=0, alternative="great", conf.level= 0.99)
t.test(h2b[h2b!=0], mu=0, alternative="great", conf.level= 0.99)
t.test(h3[h3!=0], mu=0, alternative="great", conf.level= 0.99)
t.test(h4[h4!=0], mu=0, alternative="great", conf.level= 0.99)

