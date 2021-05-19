#R-script to diplay DNA hbonds
library(matrixStats)
file = "all"
png(file="./amber/DNA_contacts_4500ns/dna_prot_contacts_half_combined.png",width=8,height=6,units="in",res=300)

df1<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_1aoi_all.dat",sep=""),header=TRUE)
df2<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_1eqz_all.dat",sep=""),header=TRUE)
df3<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_ext2_all.dat",sep=""),header=TRUE)
df4<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_sym_all.dat",sep=""),header=TRUE)
df5<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_shuxiang_ext_sym_amber_run1_all.dat",sep=""),header=TRUE)
df6<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_shuxiang_ext_sym_amber_run2_all.dat",sep=""),header=TRUE)

df1<-as.matrix(df1[,-1])
df2<-as.matrix(df2[,-1])
df3<-as.matrix(df3[,-1])
df4<-as.matrix(df4[,-1])
df5<-as.matrix(df5[,-1])
df6<-as.matrix(df6[,-1])

df1 <- df1[seq(1, nrow(df1), 5), ] ## 1ns interval of frames
df2 <- df2[seq(1, nrow(df2), 5), ]
df3 <- df3[seq(1, nrow(df3), 5), ]
df4 <- df4[seq(1, nrow(df4), 5), ]
df5 <- df5[seq(1, nrow(df5), 5), ]
df6 <- df6[seq(1, nrow(df6), 5), ]

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
mean11 <- colMeans(df_half[0:4300,])
mean12 <- colMeans(df_half[4300:(4300+750),])
mean13 <- colMeans(df_half[(4300+750):(4300+750*2),])
mean14 <- colMeans(df_half[(4300+750*2):(4300+750*3),])
mean15 <- colMeans(df_half[(4300+750*3):(4300+750*4),])

mean21 <- colMeans(df_half[7300:(7300+4800),])
mean22 <- colMeans(df_half[12100:(12100+750),])
mean23 <- colMeans(df_half[(12100+750):(12100+750*2),])
mean24 <- colMeans(df_half[(12100+750*2):(12100+750*3),])
mean25 <- colMeans(df_half[(12100+750*3):(12100+750*4),])

mean31 <- colMeans(df_half[15100:(15100+4300),])
mean32 <- colMeans(df_half[19400:(19400+750),])
mean33 <- colMeans(df_half[(19400+750):(19400+750*2),])
mean34 <- colMeans(df_half[(19400+750*2):(19400+750*3),])
mean35 <- colMeans(df_half[(19400+750*3):(19400+750*4),])

mean41 <- colMeans(df_half[22400:(22400+3800),])
mean42 <- colMeans(df_half[26200:(26200+750),])
mean43 <- colMeans(df_half[(26200+750):(26200+750*2),])
mean44 <- colMeans(df_half[(26200+750*2):(26200+750*3),])
mean45 <- colMeans(df_half[(26200+750*3):(26200+750*4),])

mean51 <- colMeans(df_half[29200:(29200+4800),])
mean61 <- colMeans(df_half[(29200+4800):(29200+4800*2),])

sderr <- colSds(rbind(mean11,mean12,mean13,mean14,mean15,
                      mean21,mean22,mean23,mean24,mean25,
                      mean31,mean32,mean33,mean34,mean35,
                      mean41,mean42,mean43,mean44,mean45,
                      mean51,mean61))/sqrt(22)


error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

#barx <- barplot(data, beside=FALSE, 
#                names.arg=0:93, col=c("yellow","red","blue","green"), axis.lty=1, font.axis=2,font.lab=2,
#                cex.lab=1.6,cex.axis=1.3,
#                xlab="Nucleosomal DNA position (SHL)" , ylab="Mean number of contacts", xaxt = "none")
plot( seq(0,93,1), mean , type="l" , bty="l" , xlab="Superhelical Location (SHL)" , ylab="Mean number of contacts", xaxt = "none",
      col="black" , lwd=3 , pch=15 , ylim=c(0,22) ,xaxt = "none", 
      cex.lab=1.8,cex.axis=1.4,font.axis=2,font.lab=2)
arrows(seq(0,93,1), mean-sderr/2, seq(0,93,1), mean+sderr/2, length=0.02, angle=90, code=3, col="darkgrey",font=2,lwd=2)


axis(1, at = seq(0, 73.5, by = 10.5),labels = c("0","±1","±2","±3","±4","±5","±6", "±7"), las=1, 
     font.axis=2,cex.axis=1.4,cex.lab=1.8)
box(which = "plot",lwd =2)
dev.off()


write.csv(mean,"./amber/DNA_contacts_4500ns/contact_mean_half.csv", row.names = FALSE)
