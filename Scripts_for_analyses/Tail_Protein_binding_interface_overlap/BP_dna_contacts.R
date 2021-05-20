#R-script to diplay DNA hbonds
library(matrixStats)
png(file="./BP_DNA_contact.png",width=10,height=5,units="in",res=300)
## protein-DNA mean contacts from PDB structures 
df<-read.table("./BP_DNA_contact.dat",header=TRUE)
contacts <-abs(as.matrix(df[,-1]))

mean<-colMeans(contacts)
dev<-colSds(contacts)
sderr <-dev/sqrt(length(contacts[,1]))

par(mar = c(5, 5, 3, 3))
plot(seq(-93,93,1), mean , type="l" , bty="l" ,pch=18, xlab="Superhelical Location (SHL)" , ylab="Mean number of contacts" , 
     col="Blue" , lwd=4 ,ylim = c(0,8),xaxt="n",cex.lab=1.8, cex.axis=1.6,font.lab=2,font.axis=2)
arrows(seq(-93,93,1), mean-sderr/2, seq(-93,93,1), mean+sderr/2, length=0.02, angle=90, code=3, col="red",font=2,lwd=2)
axis(1, at = seq(-73.5, 73.5, by = 10.5),labels = c("-7","-6","-5","-4","-3","-2","-1",
                                                   "0","1","2","3","4","5","6", "7"), las=1,cex.axis=1.4,font.axis=2)
box(which = "plot",lwd=2)
dev.off()

