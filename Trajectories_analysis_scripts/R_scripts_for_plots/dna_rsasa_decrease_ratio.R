#R-script to diplay DNA hbonds
library(matrixStats)
file = "all"
png(file=paste("./amber/sasa_4500ns/rsasa_all_half_conformtions_",file,".png",sep=""),width=8,height=6,units="in",res=300)
df11<-read.table(paste("./amber/sasa_4500ns/rsasa_",file,"_1aoi.dat",sep=""),header=TRUE)
df12<-read.table(paste("./amber/sasa_4500ns/rsasa_",file,"_1eqz.dat",sep=""),header=TRUE)
df13<-read.table(paste("./amber/sasa_4500ns/rsasa_",file,"_ext2.dat",sep=""),header=TRUE)
df14<-read.table(paste("./amber/sasa_4500ns/rsasa_",file,"_sym.dat",sep=""),header=TRUE)
df15<-read.table(paste("./amber/sasa_4500ns/rsasa_",file,"_shuxiang_ext_sym_amber_run1.dat",sep=""),header=TRUE)
df16<-read.table(paste("./amber/sasa_4500ns/rsasa_",file,"_shuxiang_ext_sym_amber_run2.dat",sep=""),header=TRUE)



df11<-as.matrix(df11[,-1])
df12<-as.matrix(df12[,-1])
df13<-as.matrix(df13[,-1])
df14<-as.matrix(df14[,-1])
df15<-as.matrix(df15[,-1])
df16<-as.matrix(df16[,-1])




df1<-rbind(df11,df12,df13,df14,df15,df16)


df1_half=df1[,94]

for (i in seq(1:93))
{
  df1_half <- cbind(df1_half,(df1[,94-i] + df1[,94+i])/2 ) 
  
}


counts <- vector()

for(i in 1:ncol(df1)){
  counts[i] <- length(df1[df1[,i] >=0.2 ,i]) 
}
counts_half=counts[94]

for (i in seq(1:93))
{
  counts_half <- cbind(counts_half,(counts[94-i] + counts[94+i])/2) 
  
}
counts_half = counts_half/length(df1[,1])

mean1<-colMeans(df1_half)
dev1<-colSds(df1_half)
sderr1 <-dev1/sqrt(length(df1_half[,1]))

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

#barx <- barplot(-mean1, beside=FALSE, 
#                names.arg=0:93, col=c("#008080"), axis.lty=1, font.axis=2,font.lab=2,
#                cex.lab=1.5,cex.axis=1.3,
#                xlab="Nucleosomal DNA position (SHL)" , ylab="Percentage of SASA Decrease ", xaxt = "none",ylim = c(-0.25,0))
par(mar =c(5,5,2,1))
plot(seq(0,93,1), counts_half , type="b" , bty="l" , xlab="Superhelical Location (SHL)" , ylab="Fraction of Frames", xaxt = "none",
     col="black" , lwd=3 , pch=15 , ylim=c(0,0.75) ,xaxt = "none", 
     cex.lab=2,cex.axis=2,font.axis=2,font.lab=2)
#arrows(seq(0,93,1), -mean1-dev1/2, seq(0,93,1), -mean1+dev1/2, 
#       length=0.02, angle=90, code=3, col="red",font=2,lwd=2)
#abline(h=0.3,lty="dashed",lwd=2)

axis(1, at = seq(0, 73.5, by = 10.5),labels = c("0","±1","±2","±3","±4","±5","±6", "±7"), las=1, 
     font.axis=2,cex.axis=2,cex.lab=2)
box(which = "plot",lwd =2)
dev.off()

