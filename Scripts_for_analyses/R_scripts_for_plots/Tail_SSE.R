#R-script to diplay DNA hbonds
library(matrixStats)
library(plyr)
tailcopy1 = "h2a_C"
png(file="Tail_SSE.png",width=7,height=4,units="in",res=300)
df11<-read.table(paste("SSE_modelA_",tailcopy1,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"), header=TRUE)
df21<-read.table(paste("SSE_modelB_",tailcopy1,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)
df31<-read.table(paste("SSE_modelC_",tailcopy1,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)
df41<-read.table(paste("SSE_modelD_",tailcopy1,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)
df51<-read.table(paste("SSE_modelD_gromacs_run1_",tailcopy1,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)
df61<-read.table(paste("SSE_modelD_gromacs_run2_",tailcopy1,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)

tailcopy2 = "h2a_G"
df12<-read.table(paste("SSE_modelA_",tailcopy2,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)
df22<-read.table(paste("SSE_modelB_",tailcopy2,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)
df32<-read.table(paste("SSE_modelC_",tailcopy2,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)
df42<-read.table(paste("SSE_modelD_",tailcopy2,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)
df52<-read.table(paste("SSE_modelD_gromacs_run1_",tailcopy2,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)
df62<-read.table(paste("SSE_modelD_gromacs_run2_",tailcopy2,".dat",sep=""),stringsAsFactors=FALSE, 
                colClasses = c("character"),header=TRUE)

df = rbind(df11,df21,df31,df41,df51,df61,df12,df22,df32,df42,df52,df62)

df_count <- data.frame("T" = c(rep(0, ncol(df[,-1]))), "E" = c(rep(0, ncol(df[,-1]))), "B" = c(rep(0, ncol(df[,-1]))),
                       "H" = c(rep(0, ncol(df[,-1]))), "G" = c(rep(0, ncol(df[,-1]))), "I" = c(rep(0, ncol(df[,-1]))),
                       "C" = c(rep(0, ncol(df[,-1]))))

i =0
## count SSE portion  of frames for each tail residues
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
                args.legend=list(
                        x=ncol(df_count) +7.5,
                        y=max(colSums(df_count)),
                        bty = "n", cex =2))
mtext(side = 1, text = "H2A-N tail", line = 5,font = 2,cex = 2, padj = 0.5)
dev.off()

