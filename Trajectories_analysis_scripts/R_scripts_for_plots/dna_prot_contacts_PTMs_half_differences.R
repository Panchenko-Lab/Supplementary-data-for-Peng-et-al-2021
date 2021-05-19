#R-script to diplay DNA hbonds
library(matrixStats)

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

### mean for WT tail
df11 <- df11[seq(1, nrow(df11), 5), ] ## 1ns interval of frames
df12 <- df12[seq(1, nrow(df12), 5), ]
df13 <- df13[seq(1, nrow(df13), 5), ]
df14 <- df14[seq(1, nrow(df14), 5), ]

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
df61 <- df21[seq(1, nrow(df21), 5), ] ## 1ns interval of frames
df62 <- df22[seq(1, nrow(df22), 5), ]
df63 <- df23[seq(1, nrow(df23), 5), ]
df64 <- df24[seq(1, nrow(df24), 5), ]

df1_half=df21[,94]
df2_half=df22[,94]
df3_half=df23[,94]
df4_half=df24[,94]

for (i in seq(1:93))
{
  df1_half <- cbind(df1_half,(df21[,94-i] + df21[,94+i])/2 ) 
  df2_half <- cbind(df2_half,(df22[,94-i] + df22[,94+i])/2 ) 
  df3_half <- cbind(df3_half,(df23[,94-i] + df23[,94+i])/2 ) 
  df4_half <- cbind(df4_half,(df24[,94-i] + df24[,94+i])/2 ) 
}
mean1<-colMeans(df1_half)
mean2<-colMeans(df2_half)
mean3<-colMeans(df3_half)
mean4<-colMeans(df4_half)

## Autocorrelation
#h2a
mean11 <- colMeans(df1_half[1:1400,])
mean12 <- colMeans(df1_half[1400:(1400+750),])
mean13 <- colMeans(df1_half[(1400+750):(1400+750*2),])
mean14 <- colMeans(df1_half[(1400+750*2):(1400+750*3),])
mean15 <- colMeans(df1_half[(1400+750*3):(1400+750*4),])
#h2b
mean21 <- colMeans(df2_half[1:1400,])
mean22 <- colMeans(df2_half[1400:(1400+750),])
mean23 <- colMeans(df2_half[(1400+750):(1400+750*2),])
mean24 <- colMeans(df2_half[(1400+750*2):(1400+750*3),])
mean25 <- colMeans(df2_half[(1400+750*3):(1400+750*4),])
#h3
mean31 <- colMeans(df3_half[1:1400,])
mean32 <- colMeans(df3_half[1400:(1400+750),])
mean33 <- colMeans(df3_half[(1400+750):(1400+750*2),])
mean34 <- colMeans(df3_half[(1400+750*2):(1400+750*3),])
mean35 <- colMeans(df3_half[(1400+750*3):(1400+750*4),])
#h4
mean41 <- colMeans(df4_half[1:1400,])
mean42 <- colMeans(df4_half[1400:(1400+750),])
mean43 <- colMeans(df4_half[(1400+750):(1400+750*2),])
mean44 <- colMeans(df4_half[(1400+750*2):(1400+750*3),])
mean45 <- colMeans(df4_half[(1400+750*3):(1400+750*4),])


sderr1 <- colSds(rbind(mean11,mean12,mean13,mean14,mean15))/sqrt(5)
sderr2 <- colSds(rbind(mean21,mean22,mean23,mean24,mean15))/sqrt(5)
sderr3 <- colSds(rbind(mean31,mean32,mean33,mean34,mean15))/sqrt(5)
sderr4 <- colSds(rbind(mean41,mean42,mean43,mean44,mean15))/sqrt(5)

png(file="./amber/PTMs_4500ns/DNA_contacts/dna_prot_contacts_half_ARG_ala_difference.png",width=8,height=6,units="in",res=300)

plot(seq(0,93,1), mean1-mean1_wt , type="l" , bty="l" , xlab="" , ylab="", xaxt = "n",
      col="gold" , lwd=3 , pch=15 ,xaxt = "none", ylim=c(-5,7.5), #yaxt = "none",
      cex.lab=2,cex.axis=2,font.axis=2,font.lab=2)
#ytick<-seq(0, 15, by=5)
#axis(side=2, at=ytick, labels = FALSE)
arrows(seq(0,93,1), mean1 - mean1_wt-sderr1/2, seq(0,93,1), mean1 - mean1_wt+sderr1/2,
       length=0.02, angle=90, code=3, col="lightgoldenrod",font=3,lwd=4,shade = TRUE) #lightgoldenrod

lines(seq(0,93,1),mean2 -mean2_wt , col="red" , lwd=3 , pch=16 , type="l" )
arrows(seq(0,93,1), mean2 -mean2_wt-sderr2/2, seq(0,93,1), mean2 -mean2_wt+sderr2/2,
       length=0.02, angle=90, code=3, col="pink",font=3,lwd=4) #pink

lines(seq(0,93,1), mean3 -mean3_wt , col="blue" , lwd=3 , pch=17 , type="l" )
arrows(seq(0,93,1), mean3 -mean3_wt-sderr3/2, seq(0,93,1), mean3 -mean3_wt+sderr3/2, 
       length=0.02, angle=90, code=3, col="lightblue",font=3,lwd=4) #lightblue

lines(seq(0,93,1), mean4 -mean4_wt , col="green" , lwd=3 , pch=18 , type="l" )
arrows(seq(0,93,1), mean4 -mean4_wt-sderr4/2, seq(0,93,1), mean4 -mean4_wt+sderr4/2, 
       length=0.02, angle=90, code=3, col="lightgreen",font=3,lwd=4) #lightgreen

abline(h = 0,lty=1, col="black",lwd=2)
if(1) {
legend(-5,8.5, 
       legend = c("H2A","H2B","H3","H4"), 
       col = c("gold","red","blue","green"), 
       pch = c(16,16,16,16), 
       bty = "n", 
       pt.cex = 3, 
       cex = 2, text.font = 2,
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1,0.1, 0.1))
}
axis(1, at = seq(0, 73.5, by = 10.5),labels = c("0","±1","±2","±3","±4","±5","±6", "±7"), las=1, 
     font.axis=2,cex.axis=2,cex.lab=2)
#axis(1, at = seq(0, 73.5, by = 10.5),labels = FALSE, las=1, 
#     font.axis=2,cex.axis=1.1,cex.lab=1.6)
box(which = "plot",lwd =2)
dev.off()
