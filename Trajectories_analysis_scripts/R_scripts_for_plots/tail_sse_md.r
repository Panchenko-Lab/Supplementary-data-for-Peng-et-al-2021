#R-script to diplay DNA hbonds
library(matrixStats)
library(plyr)
tail1 = "h2a_C"
png(file=paste("./amber/SSE_4500ns/tail_SSE_all_h2a_N.png",sep=""),width=7,height=4,units="in",res=300)
df11<-read.table(paste("./amber/SSE_4500ns/tail_SSE_","1aoi_",tail1,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"), header=TRUE)
df21<-read.table(paste("./amber/SSE_4500ns/tail_SSE_","1eqz_",tail1,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)
df31<-read.table(paste("./amber/SSE_4500ns/tail_SSE_","ext2_",tail1,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)
df41<-read.table(paste("./amber/SSE_4500ns/tail_SSE_","sym_",tail1,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)
df51<-read.table(paste("./amber/SSE_4500ns/tail_SSE_","shuxiang_ext_sym_amber_run1_",tail1,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)
df61<-read.table(paste("./amber/SSE_4500ns/tail_SSE_","shuxiang_ext_sym_amber_run2_",tail1,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)

tail2 = "h2a_G"
df12<-read.table(paste("./amber/SSE_4500ns/tail_SSE_","1aoi_",tail2,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)
df22<-read.table(paste("./amber/SSE_4500ns/tail_SSE_","1eqz_",tail2,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)
df32<-read.table(paste("./amber/SSE_4500ns/tail_SSE_","ext2_",tail2,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)
df42<-read.table(paste("./amber/SSE_4500ns/tail_SSE_","sym_",tail2,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)
df52<-read.table(paste("./amber/SSE_4500ns/tail_SSE_","shuxiang_ext_sym_amber_run1_",tail2,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)
df62<-read.table(paste("./amber/SSE_4500ns/tail_SSE_","shuxiang_ext_sym_amber_run2_",tail2,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)

df = rbind(df11,df21,df31,df41,df51,df61,df12,df22,df32,df42,df52,df62)
#df = rbind(df1,df2,df3,df4,df5,df6)
#x = c(seq(1, 11, by = 1),seq(119, 128, by = 1))
df_count <- data.frame("T" = c(rep(0, ncol(df[,-1]))), "E" = c(rep(0, ncol(df[,-1]))), "B" = c(rep(0, ncol(df[,-1]))),
                       "H" = c(rep(0, ncol(df[,-1]))), "G" = c(rep(0, ncol(df[,-1]))), "I" = c(rep(0, ncol(df[,-1]))),
                       "C" = c(rep(0, ncol(df[,-1]))))
i =0
for (name in colnames(df[,-1])){
i=i+1
counts = count(df[name]) 
for (j in 1:nrow(counts)){
df_count[i,counts[j,1]] = counts[j,2]
}
}
df_count = as.matrix(t(df_count/nrow(df)))

par(mar=c(7,5,3.5,4))
barx <- barplot(df_count[,1:13], beside=FALSE, legend.text=c("T","E","B","H","G","I","C"), 
                  axis.lty=1,  
                ylab="Secondary Structure", cex.names  = 2, font.axis = 2, font.lab = 2,cex.lab=2, las=3, cex.axis = 2,
#h2a_N                
                names.arg=c("S1","G2","R3","G4","K5","Q6","G7","G8","K9","T10","R11","A12", "K13"),
#h2a_C                
#                names.arg=c("K119","T120","E121","S122","S123","K124","S125","K126","S127","K128"),

#h2b
#                names.arg=c("A1","K2","S3","A4","P5","6A","P7","K8","K9","G10","S11","K12","K13","A14","V15","T16","K17","T18","Q19","K20", "K21", "D22", "G23"),

#h3
#                names.arg=c("A1","R2","T3","K4","Q5","T6","A7","R8","K9","S10","T11","G12","G13","K14","A15","P16","R17","K18","Q19","L20","A21","T22","K23","A24","A25","R26","K27","S28","A29","P30","A31","T32","G33","G34","V35","K36"),
#h4
#                names.arg=c("S1","G2","R3","G4","K5","G6","G7","K8","G9","L10","G11","K12","G13","G14","A15","K16","R17","H18","R19","K20"),
                col=c("green","blue","red","orange","yellow","purple","cyan"),
#-8.5 7.5, 6.6, 12,
                args.legend=list(
                        x=ncol(df_count) +7.5,
                        y=max(colSums(df_count)),
                        bty = "n", cex =2))
mtext(side = 1, text = "H2A-N tail", line = 5,font = 2,cex = 2, padj = 0.5)
#title(main=paste("SSE of the ", substr(tail1,1,2), " tail "))
dev.off()

