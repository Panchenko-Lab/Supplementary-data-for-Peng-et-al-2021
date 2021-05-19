#R-script to diplay DNA hbonds
library(matrixStats)
file = "all"
png(file="./amber/DNA_contacts_4500ns/dna_prot_contacts_all_binary.png",width=10,height=3,units="in",res=300)

df11<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_1aoi_h2a.dat",sep=""),header=TRUE)
df12<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_1aoi_h2b.dat",sep=""),header=TRUE)
df13<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_1aoi_h3.dat",sep=""),header=TRUE)
df14<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_1aoi_h4.dat",sep=""),header=TRUE)
df15<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_1aoi_all.dat",sep=""),header=TRUE)


df21<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_1eqz_h2a.dat",sep=""),header=TRUE)
df22<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_1eqz_h2b.dat",sep=""),header=TRUE)
df23<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_1eqz_h3.dat",sep=""),header=TRUE)
df24<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_1eqz_h4.dat",sep=""),header=TRUE)
df25<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_1eqz_all.dat",sep=""),header=TRUE)


df31<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_ext2_h2a.dat",sep=""),header=TRUE)
df32<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_ext2_h2b.dat",sep=""),header=TRUE)
df33<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_ext2_h3.dat",sep=""),header=TRUE)
df34<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_ext2_h4.dat",sep=""),header=TRUE)
df35<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_ext2_all.dat",sep=""),header=TRUE)

df41<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_sym_h2a.dat",sep=""),header=TRUE)
df42<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_sym_h2b.dat",sep=""),header=TRUE)
df43<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_sym_h3.dat",sep=""),header=TRUE)
df44<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_sym_h4.dat",sep=""),header=TRUE)
df45<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_sym_all.dat",sep=""),header=TRUE)

df51<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_shuxiang_ext_sym_amber_run1_h2a.dat",sep=""),header=TRUE)
df52<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_shuxiang_ext_sym_amber_run1_h2b.dat",sep=""),header=TRUE)
df53<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_shuxiang_ext_sym_amber_run1_h3.dat",sep=""),header=TRUE)
df54<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_shuxiang_ext_sym_amber_run1_h4.dat",sep=""),header=TRUE)
df55<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_shuxiang_ext_sym_amber_run1_all.dat",sep=""),header=TRUE)

df61<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_shuxiang_ext_sym_amber_run2_h2a.dat",sep=""),header=TRUE)
df62<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_shuxiang_ext_sym_amber_run2_h2b.dat",sep=""),header=TRUE)
df63<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_shuxiang_ext_sym_amber_run2_h3.dat",sep=""),header=TRUE)
df64<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_shuxiang_ext_sym_amber_run2_h4.dat",sep=""),header=TRUE)
df65<-read.table(paste("./amber/DNA_contacts_4500ns/dna_prot_contacts_shuxiang_ext_sym_amber_run2_all.dat",sep=""),header=TRUE)


df11<-as.matrix(df11[,-1])
df12<-as.matrix(df12[,-1])
df13<-as.matrix(df13[,-1])
df14<-as.matrix(df14[,-1])
df15<-as.matrix(df15[,-1])

df21<-as.matrix(df21[,-1])
df22<-as.matrix(df22[,-1])
df23<-as.matrix(df23[,-1])
df24<-as.matrix(df24[,-1])
df25<-as.matrix(df25[,-1])

df31<-as.matrix(df31[,-1])
df32<-as.matrix(df32[,-1])
df33<-as.matrix(df33[,-1])
df34<-as.matrix(df34[,-1])
df35<-as.matrix(df35[,-1])

df41<-as.matrix(df41[,-1])
df42<-as.matrix(df42[,-1])
df43<-as.matrix(df43[,-1])
df44<-as.matrix(df44[,-1])
df45<-as.matrix(df45[,-1])

df51<-as.matrix(df51[,-1])
df52<-as.matrix(df52[,-1])
df53<-as.matrix(df53[,-1])
df54<-as.matrix(df54[,-1])
df55<-as.matrix(df55[,-1])

df61<-as.matrix(df61[,-1])
df62<-as.matrix(df62[,-1])
df63<-as.matrix(df63[,-1])
df64<-as.matrix(df64[,-1])
df65<-as.matrix(df65[,-1])

df1<-rbind(df11,df21,df31,df41,df51,df61)
df2<-rbind(df12,df22,df32,df42,df52,df62)
df3<-rbind(df13,df23,df33,df43,df53,df63)
df4<-rbind(df14,df24,df34,df44,df54,df64)
df5<-rbind(df15,df25,df35,df45,df55,df65)

df1_half=df1[,94]
df2_half=df2[,94]
df3_half=df3[,94]
df4_half=df4[,94]
df5_half=df5[,94]

for (i in seq(1:93))
{
  df1_half <- cbind(df1_half,(df1[,94-i] + df1[,94+i])/2 ) 
  df2_half <- cbind(df2_half,(df2[,94-i] + df2[,94+i])/2 ) 
  df3_half <- cbind(df3_half,(df3[,94-i] + df3[,94+i])/2 ) 
  df4_half <- cbind(df4_half,(df4[,94-i] + df4[,94+i])/2 )
  df5_half <- cbind(df5_half,(df5[,94-i] + df5[,94+i])/2 ) 
  
}

mean1<-colMeans(df1_half)
mean2<-colMeans(df2_half)
mean3<-colMeans(df3_half)
mean4<-colMeans(df4_half)
mean5<-colMeans(df5_half)

dev1<-colSds(df1_half)
sderr1 <-dev1/sqrt(length(df1_half[,1]))
dev2<-colSds(df2_half)
sderr2 <-dev2/sqrt(length(df2_half[,1]))
dev3<-colSds(df3_half)
sderr3 <-dev3/sqrt(length(df3_half[,1]))
dev4<-colSds(df4_half)
sderr4 <-dev4/sqrt(length(df4_half[,1]))
dev5<-colSds(df5_half)
sderr5 <-dev5/sqrt(length(df5_half[,1]))


mean5 <- append(rev(mean5), mean5[-1]) 
dev5 <- append(rev(dev5), dev5[-1]) 
mean5[which(mean5 <5)] = 0
mean5[which(mean5 >=5)] = 5

#mean1<-colMeans(df1)
#mean2<-colMeans(df2)
#mean3<-colMeans(df3)
#mean4<-colMeans(df4)
#mean5<-colMeans(df5)

#dev1<-colSds(df1)
#sderr1 <-dev1/sqrt(length(df1[,1]))
#dev2<-colSds(df2)
#sderr2 <-dev2/sqrt(length(df2[,1]))
#dev3<-colSds(df3)
#sderr3 <-dev3/sqrt(length(df3[,1]))
#dev4<-colSds(df4)
#sderr4 <-dev4/sqrt(length(df4[,1]))

#dev5<-colSds(df5)
#sderr5 <-dev5/sqrt(length(df5[,1]))

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

#barx <- barplot(data, beside=FALSE, 
#                names.arg=0:93, col=c("yellow","red","blue","green"), axis.lty=1, font.axis=2,font.lab=2,
#                cex.lab=1.6,cex.axis=1.3,
#                xlab="Nucleosomal DNA position (SHL)" , ylab="Mean number of contacts", xaxt = "none")
plot( seq(-93,93,1),mean5 , type="l" , bty="l" , xlab="" , ylab="",
      col="blue" , lwd=2 , pch=15 , ylim=c(0,5) ,xaxt = "none", yaxt = "none", 
      cex.lab=2,cex.axis=2,font.axis=2,font.lab=2)
#arrows(seq(-93,93,1), mean5-dev5/2, seq(-93,93,1), mean5+dev5/2, length=0.02, angle=90, code=3, col="red",font=2,lwd=1.5)

#lines(mean2 , col="red" , lwd=3 , pch=16 , type="l" )
#arrows(seq(-93,93,1), mean2-dev2/2, seq(-93,93,1), mean2+dev2/2, length=0.02, angle=90, code=3, col="red",font=2,lwd=2)

#lines(mean3 , col="blue" , lwd=3 , pch=17 , type="l" )
#arrows(seq(-93,93,1), mean3-dev3/2, seq(-93,93,1), mean3+dev3/2, length=0.02, angle=90, code=3, col="blue",font=2,lwd=2)

#lines(mean4 , col="green" , lwd=3 , pch=18 , type="l" )
#arrows(seq(-93,93,1), mean4-dev4/2, seq(-93,93,1), mean4+dev4/2, length=0.02, angle=90, code=3, col="green",font=2,lwd=2)

#legend("topright", 
#       legend = c("H2A","H2B","H3","H4"), 
#       col = c("yellow","red","blue","green"), 
#       pch = c(16,16,16,16), 
#       bty = "n", 
#       pt.cex = 3, 
#       cex = 1.5, 
#       text.col = "black", 
#       horiz = F , 
#       inset = c(0.1, 0.1,0.1, 0.1))
axis(1, at = seq(-73.5, 73.5, by = 10.5),labels = c("-7","-6","-5","-4","-3","-2", "-1","0","1","2","3","4","5","6", "7"), las=1, 
     font.axis=2,cex.axis=1.4,cex.lab=2)
axis(2, at = c(0,5),labels = c("0","5"), las=3, 
     font.axis=2,cex.axis=1.6,cex.lab=2)

box(which = "plot",lwd =2)
dev.off()

#contact_mean=mean1 + mean2 + mean3 +mean4
#write.csv(contact_mean,"./amber/DNA_contacts/contact_mean.csv", row.names = FALSE)
