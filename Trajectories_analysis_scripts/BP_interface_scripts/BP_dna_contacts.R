#R-script to diplay DNA hbonds
library(matrixStats)
png(file="./BP_DNA_contact.png",width=10,height=5,units="in",res=300)
path <- "/Users/yunhuipeng/Desktop/Interactome_paper/tail_paper/figure2/DNA_contacts_4500ns/"
##tail-DNA contacts from MD simulations
df11<-read.table(paste(path,"dna_prot_contacts_1aoi_all.dat",sep=""),header=TRUE)
df12<-read.table(paste(path,"dna_prot_contacts_1eqz_all.dat",sep=""),header=TRUE)
df13<-read.table(paste(path,"dna_prot_contacts_ext2_all.dat",sep=""),header=TRUE)
df14<-read.table(paste(path,"dna_prot_contacts_sym_all.dat",sep=""),header=TRUE)
df15<-read.table(paste(path,"dna_prot_contacts_shuxiang_ext_sym_amber_run1_all.dat",sep=""),header=TRUE)
df16<-read.table(paste(path,"dna_prot_contacts_shuxiang_ext_sym_amber_run2_all.dat",sep=""),header=TRUE)


df11<-as.matrix(df11[,-1])
df12<-as.matrix(df12[,-1])
df13<-as.matrix(df13[,-1])
df14<-as.matrix(df14[,-1])
df15<-as.matrix(df15[,-1])
df16<-as.matrix(df16[,-1])

df11 <- df11[seq(1, nrow(df11), 5), ] ## 1ns interval of frames
df12 <- df12[seq(1, nrow(df12), 5), ]
df13 <- df13[seq(1, nrow(df13), 5), ]
df14 <- df14[seq(1, nrow(df14), 5), ]
df15 <- df15[seq(1, nrow(df15), 5), ]
df16 <- df16[seq(1, nrow(df16), 5), ]

df1<-rbind(df11,df12,df13,df14,df15,df16)

df1_half=df1[,94]


for (i in seq(1:93))
{
  df1_half <- cbind(df1_half,(df1[,94-i] + df1[,94+i])/2 ) 
}

mean1<-colMeans(df1_half)
mean1 <- append(rev(mean1), mean1[-1]) 
mean1[which(mean1 <5)] = 0
mean1[which(mean1 >=5)] = 8

## protein-DNA mean contacts
df<-read.table("./BP_DNA_contact.dat",header=TRUE)
contacts <-abs(as.matrix(df[,-1]))
min<-colMins(contacts)
max<-colMaxs(contacts)
mean<-colMeans(contacts)
dev<-colSds(contacts)
sderr <-dev/sqrt(length(contacts[,1]))
ratio <- colSums(contacts != 0) / length(df[,1])


#barx <- barplot(mean, beside=FALSE, 
#                names.arg=-93:93, col=c("blue"), axis.lty=1,
#                xlab="DNA Position (SHL)",
#                ylab="Number of mean contacts",xaxt="none")
#                ylab="raio of structure with contacts",xaxt="none")
#                 ylab="change of SASA by binding partner",xaxt="none")
par(mar = c(5, 5, 3, 3))
plot(seq(-93,93,1), mean , type="l" , bty="l" ,pch=18, xlab="Superhelical Location (SHL)" , ylab="Mean number of contacts" , 
     col="Blue" , lwd=4 ,ylim = c(0,8),xaxt="n",cex.lab=1.8, cex.axis=1.6,font.lab=2,font.axis=2)
arrows(seq(-93,93,1), mean-sderr/2, seq(-93,93,1), mean+sderr/2, length=0.02, angle=90, code=3, col="red",font=2,lwd=2)
#lines(seq(-93,93,1), mean1 , col="green" , lwd=2 , pch=16 , lty="solid" )
axis(1, at = seq(-73.5, 73.5, by = 10.5),labels = c("-7","-6","-5","-4","-3","-2","-1",
                                                   "0","1","2","3","4","5","6", "7"), las=1,cex.axis=1.4,font.axis=2)
box(which = "plot",lwd=2)
dev.off()

