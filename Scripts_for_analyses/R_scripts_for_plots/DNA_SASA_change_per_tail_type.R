#R-script to diplay DNA hbonds
library(matrixStats)
file = ""
png(file="DNA_SASA_change_per_tail_type.png",width=8,height=6,units="in",res=300)
df11<-read.table("dsasa_H2A_modelA.dat",header=TRUE)
df12<-read.table("dsasa_H2A_modelB.dat",header=TRUE)
df13<-read.table("dsasa_H2A_modelC.dat",header=TRUE)
df14<-read.table("dsasa_H2A_modelD.dat",header=TRUE)
df15<-read.table("dsasa_H2A_modelD_gromacs_run1.dat",header=TRUE)
df16<-read.table("dsasa_H2A_modelD_gromacs_run2.dat",header=TRUE)

df21<-read.table("dsasa_H2B_modelA.dat",header=TRUE)
df22<-read.table("dsasa_H2B_modelB.dat",header=TRUE)
df23<-read.table("dsasa_H2B_modelC.dat",header=TRUE)
df24<-read.table("dsasa_H2B_modelD.dat",header=TRUE)
df25<-read.table("dsasa_H2B_modelD_gromacs_run1.dat",header=TRUE)
df26<-read.table("dsasa_H2B_modelD_gromacs_run2.dat",header=TRUE)

df31<-read.table("dsasa_H3_modelA.dat",header=TRUE)
df32<-read.table("dsasa_H3_modelB.dat",header=TRUE)
df33<-read.table("dsasa_H3_modelC.dat",header=TRUE)
df34<-read.table("dsasa_H3_modelD.dat",header=TRUE)
df35<-read.table("dsasa_H3_modelD_gromacs_run1.dat",header=TRUE)
df36<-read.table("dsasa_H3_modelD_gromacs_run2.dat",header=TRUE)

df41<-read.table("dsasa_H4_modelA.dat",header=TRUE)
df42<-read.table("dsasa_H4_modelB.dat",header=TRUE)
df43<-read.table("dsasa_H4_modelC.dat",header=TRUE)
df44<-read.table("dsasa_H4_modelD.dat",header=TRUE)
df45<-read.table("dsasa_H4_modelD_gromacs_run1.dat",header=TRUE)
df46<-read.table("dsasa_H4_modelD_gromacs_run2.dat",header=TRUE)

df11<-as.matrix(df11[,-1])
df12<-as.matrix(df12[,-1])
df13<-as.matrix(df13[,-1])
df14<-as.matrix(df14[,-1])
df15<-as.matrix(df15[,-1])
df16<-as.matrix(df16[,-1])

df21<-as.matrix(df21[,-1])
df22<-as.matrix(df22[,-1])
df23<-as.matrix(df23[,-1])
df24<-as.matrix(df24[,-1])
df25<-as.matrix(df25[,-1])
df26<-as.matrix(df26[,-1])

df31<-as.matrix(df31[,-1])
df32<-as.matrix(df32[,-1])
df33<-as.matrix(df33[,-1])
df34<-as.matrix(df34[,-1])
df35<-as.matrix(df35[,-1])
df36<-as.matrix(df36[,-1])

df41<-as.matrix(df41[,-1])
df42<-as.matrix(df42[,-1])
df43<-as.matrix(df43[,-1])
df44<-as.matrix(df44[,-1])
df45<-as.matrix(df45[,-1])
df46<-as.matrix(df46[,-1])


df1<-rbind(df11,df12,df13,df14,df15,df16)
df2<-rbind(df21,df22,df23,df24,df25,df26)
df3<-rbind(df31,df32,df33,df34,df35,df36)
df4<-rbind(df41,df42,df43,df44,df45,df46)
dfall<-rbind(dfall1,dfall1,dfall3,dfall4,dfall5,dfall6)

## +/- SHLs are combined due to the 2-fold pseudo-symmetry 
df1_half=df1[,94]
df2_half=df2[,94]
df3_half=df3[,94]
df4_half=df4[,94]
dfall_half=dfall[,94]

for (i in seq(1:93))
{
  df1_half <- cbind(df1_half,(df1[,94-i] + df1[,94+i])/2 ) 
  df2_half <- cbind(df2_half,(df2[,94-i] + df2[,94+i])/2 ) 
  df3_half <- cbind(df3_half,(df3[,94-i] + df3[,94+i])/2 ) 
  df4_half <- cbind(df4_half,(df4[,94-i] + df4[,94+i])/2 ) 
  dfall_half <- cbind(dfall_half,(dfall[,94-i] + dfall[,94+i])/2 ) 
}

mean1<-colMeans(df1_half)
mean2<-colMeans(df2_half)
mean3<-colMeans(df3_half)
mean4<-colMeans(df4_half)

## Means for each run
#H2A
mean111 <- colMeans(df1_half[0:860,])
mean112 <- colMeans(df1_half[860:(860+150),])
mean113 <- colMeans(df1_half[(860+150):(860+150*2),])
mean114 <- colMeans(df1_half[(860+150*2):(860+150*3),])
mean115 <- colMeans(df1_half[(860+150*3):(860+150*4),])

mean121 <- colMeans(df1_half[1460:(1460+960),])
mean122 <- colMeans(df1_half[2420:(2420+150),])
mean123 <- colMeans(df1_half[(2420+150):(2420+150*2),])
mean124 <- colMeans(df1_half[(2420+150*2):(2420+150*3),])
mean125 <- colMeans(df1_half[(2420+150*3):(2420+150*4),])

mean131 <- colMeans(df1_half[3020:(3020+860),])
mean132 <- colMeans(df1_half[3880:(3880+150),])
mean133 <- colMeans(df1_half[(3880+150):(3880+150*2),])
mean134 <- colMeans(df1_half[(3880+150*2):(3880+150*3),])
mean135 <- colMeans(df1_half[(3880+150*3):(3880+150*4),])

mean141 <- colMeans(df1_half[4480:(4480+760),])
mean142 <- colMeans(df1_half[5240:(5240+150),])
mean143 <- colMeans(df1_half[(5240+150):(5240+150*2),])
mean144 <- colMeans(df1_half[(5240+150*2):(5240+150*3),])
mean145 <- colMeans(df1_half[(5240+150*3):(5240+150*4),])

mean151 <- colMeans(df1_half[5840:(5840+960),])
mean161 <- colMeans(df1_half[(5840+960):(5840+960*2),])

#H2B
mean211 <- colMeans(df2_half[0:860,])
mean212 <- colMeans(df2_half[860:(860+150),])
mean213 <- colMeans(df2_half[(860+150):(860+150*2),])
mean214 <- colMeans(df2_half[(860+150*2):(860+150*3),])
mean215 <- colMeans(df2_half[(860+150*3):(860+150*4),])

mean221 <- colMeans(df2_half[1460:(1460+960),])
mean222 <- colMeans(df2_half[2420:(2420+150),])
mean223 <- colMeans(df2_half[(2420+150):(2420+150*2),])
mean224 <- colMeans(df2_half[(2420+150*2):(2420+150*3),])
mean225 <- colMeans(df2_half[(2420+150*3):(2420+150*4),])

mean231 <- colMeans(df2_half[3020:(3020+860),])
mean232 <- colMeans(df2_half[3880:(3880+150),])
mean233 <- colMeans(df2_half[(3880+150):(3880+150*2),])
mean234 <- colMeans(df2_half[(3880+150*2):(3880+150*3),])
mean235 <- colMeans(df2_half[(3880+150*3):(3880+150*4),])

mean241 <- colMeans(df2_half[4480:(4480+760),])
mean242 <- colMeans(df2_half[5240:(5240+150),])
mean243 <- colMeans(df2_half[(5240+150):(5240+150*2),])
mean244 <- colMeans(df2_half[(5240+150*2):(5240+150*3),])
mean245 <- colMeans(df2_half[(5240+150*3):(5240+150*4),])

mean251 <- colMeans(df2_half[5840:(5840+960),])
mean261 <- colMeans(df2_half[(5840+960):(5840+960*2),])


#H3
mean311 <- colMeans(df3_half[0:860,])
mean312 <- colMeans(df3_half[860:(860+150),])
mean313 <- colMeans(df3_half[(860+150):(860+150*2),])
mean314 <- colMeans(df3_half[(860+150*2):(860+150*3),])
mean315 <- colMeans(df3_half[(860+150*3):(860+150*4),])

mean321 <- colMeans(df3_half[1460:(1460+960),])
mean322 <- colMeans(df3_half[2420:(2420+150),])
mean323 <- colMeans(df3_half[(2420+150):(2420+150*2),])
mean324 <- colMeans(df3_half[(2420+150*2):(2420+150*3),])
mean325 <- colMeans(df3_half[(2420+150*3):(2420+150*4),])

mean331 <- colMeans(df3_half[3020:(3020+860),])
mean332 <- colMeans(df3_half[3880:(3880+150),])
mean333 <- colMeans(df3_half[(3880+150):(3880+150*2),])
mean334 <- colMeans(df3_half[(3880+150*2):(3880+150*3),])
mean335 <- colMeans(df3_half[(3880+150*3):(3880+150*4),])

mean341 <- colMeans(df3_half[4480:(4480+760),])
mean342 <- colMeans(df3_half[5240:(5240+150),])
mean343 <- colMeans(df3_half[(5240+150):(5240+150*2),])
mean344 <- colMeans(df3_half[(5240+150*2):(5240+150*3),])
mean345 <- colMeans(df3_half[(5240+150*3):(5240+150*4),])

mean351 <- colMeans(df3_half[5840:(5840+960),])
mean361 <- colMeans(df3_half[(5840+960):(5840+960*2),])


#H4
mean411 <- colMeans(df4_half[0:860,])
mean412 <- colMeans(df4_half[860:(860+150),])
mean413 <- colMeans(df4_half[(860+150):(860+150*2),])
mean414 <- colMeans(df4_half[(860+150*2):(860+150*3),])
mean415 <- colMeans(df4_half[(860+150*3):(860+150*4),])

mean421 <- colMeans(df4_half[1460:(1460+960),])
mean422 <- colMeans(df4_half[2420:(2420+150),])
mean423 <- colMeans(df4_half[(2420+150):(2420+150*2),])
mean424 <- colMeans(df4_half[(2420+150*2):(2420+150*3),])
mean425 <- colMeans(df4_half[(2420+150*3):(2420+150*4),])

mean431 <- colMeans(df4_half[3020:(3020+860),])
mean432 <- colMeans(df4_half[3880:(3880+150),])
mean433 <- colMeans(df4_half[(3880+150):(3880+150*2),])
mean434 <- colMeans(df4_half[(3880+150*2):(3880+150*3),])
mean435 <- colMeans(df4_half[(3880+150*3):(3880+150*4),])

mean441 <- colMeans(df4_half[4480:(4480+760),])
mean442 <- colMeans(df4_half[5240:(5240+150),])
mean443 <- colMeans(df4_half[(5240+150):(5240+150*2),])
mean444 <- colMeans(df4_half[(5240+150*2):(5240+150*3),])
mean445 <- colMeans(df4_half[(5240+150*3):(5240+150*4),])

mean451 <- colMeans(df4_half[5840:(5840+960),])
mean461 <- colMeans(df4_half[(5840+960):(5840+960*2),])

## calculate standard error
sderr1 <- colSds(rbind(mean111,mean112,mean113,mean114,mean115,
                       mean121,mean122,mean123,mean124,mean125,
                       mean131,mean132,mean133,mean134,mean135,
                       mean141,mean142,mean143,mean144,mean145,
                       mean151,mean161))/sqrt(22)

sderr2 <- colSds(rbind(mean211,mean212,mean213,mean214,mean215,
                       mean221,mean222,mean223,mean224,mean225,
                       mean231,mean232,mean233,mean234,mean235,
                       mean241,mean242,mean243,mean244,mean245,
                       mean251,mean261))/sqrt(22)

sderr3 <- colSds(rbind(mean311,mean312,mean313,mean314,mean315,
                       mean321,mean322,mean323,mean324,mean325,
                       mean331,mean332,mean333,mean334,mean335,
                       mean341,mean342,mean343,mean344,mean345,
                       mean351,mean361))/sqrt(22)

sderr4 <- colSds(rbind(mean411,mean412,mean413,mean414,mean415,
                       mean421,mean422,mean423,mean424,mean425,
                       mean431,mean432,mean433,mean434,mean435,
                       mean441,mean442,mean443,mean444,mean445,
                       mean451,mean461))/sqrt(22)

## make plots
par(mar =c(5,5,2,1))
plot( seq(0,93,1),-mean1 , type="l" , bty="l" , xlab="Superhelical Location (SHL)" , ylab="Mean Decrease of DNA accessibility", xaxt = "none",
      col="gold" , lwd=3 , pch=15 , ylim=c(-100,0) ,xaxt = "none", yaxt = "none",
      cex.lab=1.8,cex.axis=2,font.axis=2,font.lab=2)
arrows(seq(0,93,1), -mean1-sderr1/2, seq(0,93,1), -mean1+sderr1/2, length=0.02, angle=90, code=3, col="lightgoldenrod",font=2,lwd=2)

lines(seq(0,93,1),-mean2 , col="red" , lwd=3 , pch=16 , type="l" )
arrows(seq(0,93,1), -mean2-sderr2/2, seq(0,93,1), -mean2+sderr2/2, length=0.02, angle=90, code=3, col="pink",font=2,lwd=2)

lines(seq(0,93,1),-mean3 , col="blue" , lwd=3 , pch=17 , type="l" )
arrows(seq(0,93,1), -mean3-sderr3/2, seq(0,93,1), -mean3+sderr3/2, length=0.02, angle=90, code=3, col="lightblue",font=2,lwd=2)

lines(seq(0,93,1),-mean4 , col="green" , lwd=3 , pch=18 , type="l" )
arrows(seq(0,93,1), -mean4-sderr4/2, seq(0,93,1), -mean4+sderr4/2, length=0.02, angle=90, code=3, col="lightgreen",font=2,lwd=2)

abline(h = 0,lty=1, col="black",lwd=2)
legend(70,-40, 
       legend = c("H2A","H2B","H3","H4"), 
       col = c("gold","red","blue","green"), 
       pch = c(16,16,16,16), 
       bty = "n", 
       pt.cex = 3, 
       cex = 2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1,0.1, 0.1))

axis(1, at = seq(0, 73.5, by = 10.5),labels = c("0","±1","±2","±3","±4","±5","±6", "±7"), las=1, 
     font.axis=2,cex.axis=2,cex.lab=2)
axis(2, at = seq(-80, 0, by = 20),labels = c("-80","-60","-40","-20","0"), las=3,
     font.axis=2,cex.axis=2,cex.lab=2)

box(which = "plot",lwd =2)
dev.off()
