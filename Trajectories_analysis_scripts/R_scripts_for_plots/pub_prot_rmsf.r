#R-script to analyze DNA-rmsd evolution
library(ggplot2)
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(reshape2)
library(gridExtra)
library(plyr)

df1<-read.table("./amber/rmsf_4000ns/1aoi_rmsf_chains.dat",skip=4,header=TRUE,check.name=FALSE)
df2<-read.table("./amber/rmsf_4000ns/1eqz_rmsf_chains.dat",skip=4,header=TRUE,check.name=FALSE)
df3<-read.table("./amber/rmsf_4000ns/ext2_rmsf_chains.dat",skip=4,header=TRUE,check.name=FALSE)
df4<-read.table("./amber/rmsf_4000ns/sym_rmsf_chains.dat",skip=4,header=TRUE,check.name=FALSE)

df1<-df1[,1:3]
df2<-df2[,2:3]
df3<-df3[,2:3]
df4<-df4[,2:3]

df_all <- cbind(df1,df2,df3,df4)
df <- df_all$Resid

H3_mean_rmsf <- rowMeans(df_all[,2:9])
H3_se_rmsf <- rowSds(as.matrix(df_all[,2:9]))/sqrt(8)
H3_mean_rmsf <- H3_mean_rmsf[95:144]
H3_se_rmsf <- H3_se_rmsf[95:144]
png(file="./amber/rmsf_4000ns/H3_RMSF.png",width=8,height=4,units="in",res=300)
## plot for different models
plot( seq(1,50,1),H3_mean_rmsf , type="l" , bty="l" , xlab="H3 Tail" , ylab="RMSF", #xaxt = "none",
      col="gold" , lwd=3 , pch=15 , ylim=c(0,25) ,#xaxt = "none", 
      cex.lab=1.5,cex.axis=1.3,font.axis=2,font.lab=2)
arrows(seq(1,50,1), H3_mean_rmsf-H3_se_rmsf/2, seq(1,50,1), H3_mean_rmsf+H3_se_rmsf/2, length=0.02, angle=90, code=3, col="lightgoldenrod",font=2,lwd=2)

legend(76,26, 
       legend = c("H3"), 
       col = c("blue"), 
       pch = c(16), 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.5, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1))
#axis(1, at = seq(0, 73.5, by = 10.5),labels = c("0","±1","±2","±3","±4","±5","±6", "±7"), las=1, 
#     font.axis=2,cex.axis=1.1,cex.lab=1.6)
box(which = "plot",lwd =2)
dev.off()