library(matrixStats)
tail = "h2a_C"
tail1 = "h2a_C"
tail2 = "h2a_G"
cut_off = 0.1
png(file=paste("./amber/residence_time_4500ns/time_binary/tail_contacts_time_",tail,"_",cut_off,".png",sep=""),width=12,height=14,units="in",res=300)
df11<-read.table(paste("./amber/tail_contacts_4500ns/tail_contacts_all_1aoi_",tail1,".dat",sep=""),header=TRUE)
df12<-read.table(paste("./amber/tail_contacts_4500ns/tail_contacts_all_1aoi_",tail2,".dat",sep=""),header=TRUE)

df21<-read.table(paste("./amber/tail_contacts_4500ns/tail_contacts_all_1eqz_",tail1,".dat",sep=""),header=TRUE)
df22<-read.table(paste("./amber/tail_contacts_4500ns/tail_contacts_all_1eqz_",tail2,".dat",sep=""),header=TRUE)

df31<-read.table(paste("./amber/tail_contacts_4500ns/tail_contacts_all_ext2_",tail1,".dat",sep=""),header=TRUE)
df32<-read.table(paste("./amber/tail_contacts_4500ns/tail_contacts_all_ext2_",tail2,".dat",sep=""),header=TRUE)

df41<-read.table(paste("./amber/tail_contacts_4500ns/tail_contacts_all_sym_",tail1,".dat",sep=""),header=TRUE)
df42<-read.table(paste("./amber/tail_contacts_4500ns/tail_contacts_all_sym_",tail2,".dat",sep=""),header=TRUE)


df11<-as.matrix(df11[,-1])
df12<-as.matrix(df12[,-1])
df21<-as.matrix(df21[,-1])
df22<-as.matrix(df22[,-1])
df31<-as.matrix(df31[,-1])
df32<-as.matrix(df32[,-1])
df41<-as.matrix(df41[,-1])
df42<-as.matrix(df42[,-1])

##H2A_N
#df11<-as.matrix(df11[0:21501,1:8])
#df12<-as.matrix(df12[0:21501,1:8])
#df21<-as.matrix(df21[0:24001,1:8])
#df22<-as.matrix(df22[0:24001,1:8])
#df31<-as.matrix(df31[0:21501,1:8])
#df32<-as.matrix(df32[0:21501,1:8])
#df41<-as.matrix(df41[0:19001,1:8])
#df42<-as.matrix(df42[0:19001,1:8])

##H2A_C
df11<-as.matrix(df11[0:21501,15:23])
df12<-as.matrix(df12[0:21501,15:23])
df21<-as.matrix(df21[0:24001,15:23])
df22<-as.matrix(df22[0:24001,15:23])
df31<-as.matrix(df31[0:21501,15:23])
df32<-as.matrix(df32[0:21501,15:23])
df41<-as.matrix(df41[0:19001,15:23])
df42<-as.matrix(df42[0:19001,15:23])

##H2B
#df11<-as.matrix(df11[0:21501,1:18])
#df12<-as.matrix(df12[0:21501,1:18])
#df21<-as.matrix(df21[0:24001,1:18])
#df22<-as.matrix(df22[0:24001,1:18])
#df31<-as.matrix(df31[0:21501,1:18])
#df32<-as.matrix(df32[0:21501,1:18])
#df41<-as.matrix(df41[0:19001,1:18])
#df42<-as.matrix(df42[0:19001,1:18])

##H3
#df11<-as.matrix(df11[0:21501,1:33])
#df12<-as.matrix(df12[0:21501,1:33])
#df21<-as.matrix(df21[0:24001,1:33])
#df22<-as.matrix(df22[0:24001,1:33])
#df31<-as.matrix(df31[0:21501,1:33])
#df32<-as.matrix(df32[0:21501,1:33])
#df41<-as.matrix(df41[0:19001,1:33])
#df42<-as.matrix(df42[0:19001,1:33])

##H4
#df11<-as.matrix(df11[1:21501,1:15])
#df12<-as.matrix(df12[1:21501,1:15])
#df21<-as.matrix(df21[1:24001,1:15])
#df22<-as.matrix(df22[1:24001,1:15])
#df31<-as.matrix(df31[1:21501,1:15])
#df32<-as.matrix(df32[1:21501,1:15])
#df41<-as.matrix(df41[1:19001,1:15])
#df42<-as.matrix(df42[1:19001,1:15])



df_all <- rbind(df11,df12,df21,df22,df31,df32,df41,df42)

ratio11 = rowSums(df11!= 0)/length(df11[1,])
ratio12 = rowSums(df12!= 0)/length(df12[1,])
ratio21 = rowSums(df21!= 0)/length(df21[1,])
ratio22 = rowSums(df22!= 0)/length(df22[1,])
ratio31 = rowSums(df31!= 0)/length(df31[1,])
ratio32 = rowSums(df32!= 0)/length(df32[1,])
ratio41 = rowSums(df41!= 0)/length(df41[1,])
ratio42 = rowSums(df42!= 0)/length(df42[1,])

all_ratio = rowSums(df_all!= 0)/length(df_all[1,])

counts_bind = length(which( all_ratio > cut_off, arr.ind=TRUE))
counts_unbind = length(which( all_ratio <= cut_off, arr.ind=TRUE))
print(counts_bind)
print(counts_unbind)


if(1){
ratio11[which(ratio11 >= (1-cut_off))] = 1.05
ratio11[which((ratio11 >cut_off) & (ratio11 <(1-cut_off)) )] = 0.5
ratio11[which(ratio11 <= cut_off)] = 0
ratio12[which(ratio12 >= (1-cut_off))] = 1
ratio12[which((ratio12 >cut_off) & (ratio12 <(1-cut_off)) )] = 0.45
ratio12[which(ratio12 <= cut_off)] = -0.05

ratio21[which(ratio21 >= (1-cut_off))] = 1.05
ratio21[which((ratio21 >cut_off) & (ratio21 <(1-cut_off)) )] = 0.5
ratio21[which(ratio21 <= cut_off)] = 0
ratio22[which(ratio22 >= (1-cut_off))] = 1
ratio22[which((ratio22 >cut_off) & (ratio22 <(1-cut_off)) )] = 0.45
ratio22[which(ratio22 <= cut_off)] = -0.05

ratio31[which(ratio31 >= (1-cut_off))] = 1.05
ratio31[which((ratio31 >cut_off) & (ratio31 <(1-cut_off)) )] = 0.5
ratio31[which(ratio31 <= cut_off)] = 0
ratio32[which(ratio32 >= (1-cut_off))] = 1
ratio32[which((ratio32 >cut_off) & (ratio32 <(1-cut_off)) )] = 0.45
ratio32[which(ratio32 <= cut_off)] = -0.05

ratio41[which(ratio41 >= (1-cut_off))] = 1.05
ratio41[which((ratio41 >cut_off) & (ratio41 <(1-cut_off)) )] = 0.5
ratio41[which(ratio41 <= cut_off)] = 0
ratio42[which(ratio42 >= (1-cut_off))] = 1
ratio42[which((ratio42 >cut_off) & (ratio42 <(1-cut_off)) )] = 0.45
ratio42[which(ratio42 <= cut_off)] = -0.05
}

par(mfrow=c(4,1),mar=c(5,2,5,1) + 0.1)


####
plot(seq(200,4500,0.2),  ratio11 , type="l" , bty="l" , xlab="" , ylab="" , yaxt="n",
     col="red" , lwd=1 , pch=15 , ylim=c(0,1) ,cex.lab=1.1, xaxt="none")
lines(seq(200,4500,0.2), ratio12 , col="blue" , lwd=1 , pch=16 , type="l" )
axis(1,c(200,1000,2000,3000,4000,4500),cex.axis=2,font=2)

#mtext(side=1, line=2, "Time(ns)", col="black", font=2,cex=1.2)
#mtext(side=2, line=3, "binding Ratio ", col="black", font=2, cex=1.2)
title(main="Model A",font=2.4,cex.main=2.4)
box(which = "plot", lty = "solid",lwd=2)
####
plot(seq(200,5000,0.2),  ratio21 , type="l" , bty="l" , xlab="" , ylab="" , yaxt="n",
     col="red" , lwd=1 , pch=15 , ylim=c(0,1) , cex.lab=1.1,xaxt="none")
lines(seq(200,5000,0.2), ratio22 , col="blue" , lwd=1 , pch=16 , type="l" )
axis(1,c(200,1000,2000,3000,4000,5000),cex.axis=2,font=2)


#mtext(side=1, line=2, "Time(ns)", col="black", font=2,cex=1.2)
#mtext(side=2, line=3, "binding Ratio ", col="black", font=2, cex=1.2)
title(main="Model B",font=2.4,cex.main=2.4)
box(which = "plot", lty = "solid",lwd=2)
####
plot(seq(200,4500,0.2),  ratio31 , type="l" , bty="l" , xlab="" , ylab="" , yaxt="n",
     col="red" , lwd=1 , pch=15 , ylim=c(0,1) , cex.lab=1.1,xaxt="none")
lines(seq(200,4500,0.2), ratio32 , col="blue" , lwd=1 , pch=16 , type="l" )
axis(1,c(200,1000,2000,3000,4000,4500),cex.axis=2,font=2)
#mtext(side=1, line=2, "Time(ns)", col="black", font=2,cex=1.2)
#mtext(side=2, line=3, "binding Ratio ", col="black", font=2, cex=1.2)
title(main="Model C",font=2.4,cex.main=2.4)
box(which = "plot", lty = "solid",lwd=2)
####

plot(seq(200,4000,0.2),  ratio41 , type="l" , bty="l" , xlab="" , ylab="" , yaxt="n",
     col="red" , lwd=1 , pch=15 , ylim=c(0,1) , cex.lab=1.1,xaxt="none")
lines(seq(200,4000,0.2), ratio42 , col="blue" , lwd=1 , pch=16 , type="l" )
axis(1,c(200,1000,2000,3000,4000),cex.axis=2,font=2)

mtext(side=1, line=2, "Time(ns)", col="black", font=2,cex=2.4,padj = 1)
#mtext(side=2, line=3, "binding Ratio ", col="black", font=2, cex=1.2)
title(main="Model D",font=2.4,cex.main=2.4)
box(which = "plot", lty = "solid",lwd=2)

dev.off()




