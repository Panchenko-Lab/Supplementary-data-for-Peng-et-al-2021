library(matrixStats)
file = "all"
png(file="DNA_rSASA_change_0.2.png",width=8,height=6,units="in",res=300)
df1<-read.table("rsasa_all_modelA.dat",header=TRUE)
df2<-read.table("rsasa_all_modelB.dat",header=TRUE)
df3<-read.table("rsasa_all_modelC.dat",header=TRUE)
df4<-read.table("rsasa_all_modelD.dat",header=TRUE)
df5<-read.table("rsasa_all_modelD_gromacs_run1.dat",header=TRUE)
df6<-read.table("rsasa_all_modelD_gromacs_run1.dat",header=TRUE)

df11<-as.matrix(df11[,-1])
df12<-as.matrix(df12[,-1])
df13<-as.matrix(df13[,-1])
df14<-as.matrix(df14[,-1])
df15<-as.matrix(df15[,-1])
df16<-as.matrix(df16[,-1])

df1<-rbind(df11,df12,df13,df14,df15,df16)

##SASA change for +/- SHLs are combined due to the 2-fold pseudo-symmetry 

df1_half=df1[,94]

for (i in seq(1:93))
{
  df1_half <- cbind(df1_half,(df1[,94-i] + df1[,94+i])/2 ) 
}

## counts the number of frames with percentage of SASA change > 20% per DNA basepair
counts <- vector()
for(i in 1:ncol(df1)){
  counts[i] <- length(df1[df1[,i] >=0.2 ,i]) 
}
counts_half=counts[94]

## +/- SHLs are combined due to the 2-fold pseudo-symmetry 
for (i in seq(1:93))
{
  counts_half <- cbind(counts_half,(counts[94-i] + counts[94+i])/2) 
}

counts_half = counts_half/length(df1[,1])


## make plota
par(mar =c(5,5,2,1))
plot(seq(0,93,1), counts_half , type="b" , bty="l" , xlab="Superhelical Location (SHL)" ,
     ylab="Fraction of Frames", xaxt = "none",
     col="black" , lwd=3 , pch=15 , ylim=c(0,0.75) ,xaxt = "none", 
     cex.lab=2,cex.axis=2,font.axis=2,font.lab=2)

axis(1, at = seq(0, 73.5, by = 10.5),labels = c("0","±1","±2","±3","±4","±5","±6", "±7"), las=1, 
     font.axis=2,cex.axis=2,cex.lab=2)
box(which = "plot",lwd =2)
dev.off()

