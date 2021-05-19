library(matrixStats)
tail = "h4"
tail1 = "h4_B"
tail2 = "h4_F"
type = "all"
png(file=paste("./amber/PTMs_4500ns/tail_contacts/tail_contacts_compare_",type,"_",tail,".png",sep=""),width=12,height=5,units="in",res=300)
df11<-read.table(paste("./amber/tail_contacts_4500ns/tail_contacts_",type,"_sym_",tail1,".dat",sep=""),header=TRUE)
df12<-read.table(paste("./amber/tail_contacts_4500ns/tail_contacts_",type,"_sym_",tail2,".dat",sep=""),header=TRUE)

df21<-read.table(paste("./amber/PTMs_4500ns/tail_contacts/tail_contacts_",type,"_ARG_ala_",tail1,".dat",sep=""),header=TRUE)
df22<-read.table(paste("./amber/PTMs_4500ns/tail_contacts/tail_contacts_",type,"_ARG_ala_",tail2,".dat",sep=""),header=TRUE)

df31<-read.table(paste("./amber/PTMs_4500ns/tail_contacts/tail_contacts_",type,"_LYS_Ac_",tail1,".dat",sep=""),header=TRUE)
df32<-read.table(paste("./amber/PTMs_4500ns/tail_contacts/tail_contacts_",type,"_LYS_Ac_",tail2,".dat",sep=""),header=TRUE)

df41<-read.table(paste("./amber/PTMs_4500ns/tail_contacts/tail_contacts_",type,"_LYS_Me_",tail1,".dat",sep=""),header=TRUE)
df42<-read.table(paste("./amber/PTMs_4500ns/tail_contacts/tail_contacts_",type,"_LYS_Me_",tail2,".dat",sep=""),header=TRUE)

df51<-read.table(paste("./amber/PTMs_4500ns/tail_contacts/tail_contacts_",type,"_SER_ph_",tail1,".dat",sep=""),header=TRUE)
df52<-read.table(paste("./amber/PTMs_4500ns/tail_contacts/tail_contacts_",type,"_SER_ph_",tail2,".dat",sep=""),header=TRUE)

df11<-as.matrix(rbind(df11[1:7000,-1],df11[19000:34000,-1]))
df12<-as.matrix(rbind(df12[1:7000,-1],df12[19000:34000,-1]))

#df11<-as.matrix(df11[,-1])
#df12<-as.matrix(df12[,-1])
df21<-as.matrix(df21[,-1])
df22<-as.matrix(df22[,-1])
df31<-as.matrix(df31[,-1])
df32<-as.matrix(df32[,-1])
df41<-as.matrix(df41[,-1])
df42<-as.matrix(df42[,-1])
df51<-as.matrix(df51[,-1])
df52<-as.matrix(df52[,-1])

df1 <-rbind(df11,df12)
df2 <-rbind(df21,df22)
df3 <-rbind(df31,df32)
df4 <-rbind(df41,df42)
df5 <-rbind(df51,df52)

mean1<-colMeans(df1)
mean2<-colMeans(df2)
mean3<-colMeans(df3)
mean4<-colMeans(df4)
mean5<-colMeans(df5)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

par(mar=c(5,5,1,1))
plot( mean1 , type="b" , bty="l" , xlab="" , ylab="" , 
    col="magenta" , lwd=3 , pch=15 , ylim=c(0,20) ,xaxt = "none", cex.lab=2,cex.axis=2,font=2, cex=2)
lines(mean2 , col="aqua marine" , lwd=3 , pch=16 , type="b", cex=2 )
lines(mean3 , col="orange" , lwd=3 , pch=17 , type="b",cex=2  )
lines(mean4 , col="cyan" , lwd=3 , pch=18 , type="b", cex=2  )
lines(mean5 , col="Gray" , lwd=3 , pch=20 , type="b", cex=2 )

#h2a
if(0){
plot( seq(1,13,1),mean1[1:13] , type="b" , bty="l" , xlab="" , ylab="" , 
      col="magenta" , lwd=3 , pch=15 ,xlim=c(0,23), ylim=c(0,20) ,xaxt = "none", cex.lab=2, cex.axis=2,font=2, cex=2)
lines(seq(1,13,1),mean2[1:13] , col="aqua marine" , lwd=3 , pch=16 , type="b",cex=2 )
lines(seq(1,13,1),mean3[1:13] , col="orange" , lwd=3 , pch=17 , type="b",cex=2 )
lines(seq(1,13,1),mean4[1:13] , col="cyan" , lwd=3 , pch=18 , type="b" ,cex=2)
lines(seq(1,13,1),mean5[1:13] , col="Gray" , lwd=3 , pch=20 , type="b" ,cex=2)


lines(seq(14,23,1),mean1[14:23] , col="magenta" , lwd=3 , pch=15 , type="b" ,cex=2)
lines(seq(14,23,1),mean2[14:23] , col="aqua marine" , lwd=3 , pch=16 , type="b",cex=2 )
lines(seq(14,23,1),mean3[14:23] , col="orange" , lwd=3 , pch=17 , type="b",cex=2 )
lines(seq(14,23,1),mean4[14:23] , col="cyan" , lwd=3 , pch=18 , type="b" ,cex=2)
lines(seq(14,23,1),mean5[14:23] , col="Gray" , lwd=3 , pch=20 , type="b",cex=2 )
abline(v = 13.5,lty=2, col="dark grey",lwd=3)

#h2a
axis(1, at = seq(1, 23, by = 1),labels = c("S1","G2","R3","G4","K5","Q6","G7","G8","K9","T10","R11","A12","K13","K119","T120","E121","S122",
                                          "S123","K124","S125","K126","S127","K128" ), las=2, cex.axis=2,padj = 0.5,font=2)
}
#h2b
#axis(1, at = seq(1, 23, by = 1),labels = c("A1","K2","S3","A4","P5","A6","P7","K8","K9","G10","S11","K12","K13","A14","V15",
#                                         "T16","K17","T18","Q19","K20","K21","D22","G23"), las=2, cex.axis=2,padj = 0.5,font=2)
#h3
#axis(1, at = seq(1, 36, by = 1),labels = c("A1","R2","T3","K4","Q5","T6","A7","R8","K9","S10","T11","G12","G13",
#                                           "K14","A15","P16","R17","K18","Q19","L20","A21","T22","K23","A24","A25","R26",
#                                          "K27","S28","A29","P30","A31","T32","G33","G34","V35","K36"), las=2, cex.axis=1.8,padj = 0.5,font=2)
#h4
axis(1, at = seq(1, 20, by = 1),labels = c("S1","G2","R3","G4","K5","G6","G7","K8","G9","L10",
                                           "G11","K12","G13","G14","A15","K16","R17","H18","R19","K20"),
                                          las=3, cex.axis=2,padj = 0.5,font=2)

mtext(side=1, line=2, "", col="blue", font=2,cex=2)
mtext(side=2, line=3, "Mean Contacts Number", col="black", font=2, cex=2)

#title(main=paste("Mean contacts between ", tail, " tail and DNA",seq =""))
if(1){
legend(15,22, box.lwd=2, box.col="black",
        col = c("magenta","aqua marine", "orange" ,"cyan","Gray"),
       legend = c("Unmodified Tail ", "Arg -> Ala", "Acetylation", "Methylation","Phosphorylation"), 
       pch = c(15,16,17,18,20), 
       # col = c("magenta", "orange" ,"cyan","Gray"),
       #legend = c("Unmodified Tail" , "Acetylation", "Methylation","Phosphorylation"), 
       #pch = c(15,17,18,20), 
       #bty = "n", 
       pt.cex = 3, 
       cex = 1.8, 
       text.col = "black", text.font = 2,
       horiz = F  )
}
box(which = "plot", lty = "solid",lwd=2)

dev.off()

