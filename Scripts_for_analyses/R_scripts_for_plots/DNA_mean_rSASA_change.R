library(matrixStats)
png(file="DNA_mean_rSASA_change.png",width=8,height=6,units="in",res=300)
df1<-read.table("rsasa_all_modelA.dat",header=TRUE)
df2<-read.table("rsasa_all_modelB.dat",header=TRUE)
df3<-read.table("rsasa_all_modelC.dat",header=TRUE)
df4<-read.table("rsasa_all_modelD.dat",header=TRUE)
df5<-read.table("rsasa_all_modelD_gromacs_run1.dat",header=TRUE)
df6<-read.table("rsasa_all_modelD_gromacs_run1.dat",header=TRUE)


df1<-as.matrix(df1[,-1])
df2<-as.matrix(df2[,-1])
df3<-as.matrix(df3[,-1])
df4<-as.matrix(df4[,-1])
df5<-as.matrix(df5[,-1])
df6<-as.matrix(df6[,-1])

##SASA change for +/- SHLs are combined due to the 2-fold pseudo-symmetry 

df1_half=df1[,94]
df2_half=df2[,94]
df3_half=df3[,94]
df4_half=df4[,94]
df5_half=df5[,94]
df6_half=df6[,94]

for (i in seq(1:93))
{
  df1_half <- cbind(df1_half,(df1[,94-i] + df1[,94+i])/2 ) 
  df2_half <- cbind(df2_half,(df2[,94-i] + df2[,94+i])/2 ) 
  df3_half <- cbind(df3_half,(df3[,94-i] + df3[,94+i])/2 ) 
  df4_half <- cbind(df4_half,(df4[,94-i] + df4[,94+i])/2 ) 
  df5_half <- cbind(df5_half,(df5[,94-i] + df5[,94+i])/2 ) 
  df6_half <- cbind(df6_half,(df6[,94-i] + df6[,94+i])/2 ) 
  
}

df_half<-rbind(df1_half,df2_half,df3_half,df4_half,df5_half,df6_half)
mean<-colMeans(df_half)

## Means for each run
mean11 <- colMeans(df_half[0:860,])
mean12 <- colMeans(df_half[860:(860+150),])
mean13 <- colMeans(df_half[(860+150):(860+150*2),])
mean14 <- colMeans(df_half[(860+150*2):(860+150*3),])
mean15 <- colMeans(df_half[(860+150*3):(860+150*4),])

mean21 <- colMeans(df_half[1460:(1460+960),])
mean22 <- colMeans(df_half[2420:(2420+150),])
mean23 <- colMeans(df_half[(2420+150):(2420+150*2),])
mean24 <- colMeans(df_half[(2420+150*2):(2420+150*3),])
mean25 <- colMeans(df_half[(2420+150*3):(2420+150*4),])

mean31 <- colMeans(df_half[3020:(3020+860),])
mean32 <- colMeans(df_half[3880:(3880+150),])
mean33 <- colMeans(df_half[(3880+150):(3880+150*2),])
mean34 <- colMeans(df_half[(3880+150*2):(3880+150*3),])
mean35 <- colMeans(df_half[(3880+150*3):(3880+150*4),])

mean41 <- colMeans(df_half[4480:(4480+760),])
mean42 <- colMeans(df_half[5240:(5240+150),])
mean43 <- colMeans(df_half[(5240+150):(5240+150*2),])
mean44 <- colMeans(df_half[(5240+150*2):(5240+150*3),])
mean45 <- colMeans(df_half[(5240+150*3):(5240+150*4),])

mean51 <- colMeans(df_half[5840:(5840+960),])
mean61 <- colMeans(df_half[(5840+960):(5840+960*2),])

## calculate standard error
sderr <- colSds(rbind(mean11,mean12,mean13,mean14,mean15,
                      mean21,mean22,mean23,mean24,mean25,
                      mean31,mean32,mean33,mean34,mean35,
                      mean41,mean42,mean43,mean44,mean45,
                      mean51,mean61))/sqrt(22)

plot(seq(0,93,1), -mean , type="l" , bty="l" , xlab="Superhelical Location (SHL)" , 
     ylab="Percentage of SASA Decrease", xaxt = "none",
      col="black" , lwd=3 , pch=15 , ylim=c(-0.3,0) ,xaxt = "none", 
      cex.lab=1.5,cex.axis=1.3,font.axis=2,font.lab=2)
arrows(seq(0,93,1), -mean-sderr/2, seq(0,93,1), -mean+sderr/2, 
       length=0.02, angle=90, code=3, col="darkgrey",font=2,lwd=2)


axis(1, at = seq(0, 73.5, by = 10.5),labels = c("0","±1","±2","±3","±4","±5","±6", "±7"), las=1, 
     font.axis=2,cex.axis=1.1,cex.lab=1.6)
box(which = "plot",lwd =2)
dev.off()

