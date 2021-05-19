library(VennDiagram)
png(file=paste("./structure_classfy.png",sep=""),width=8,height=6,units="in",res=300)
draw.triple.venn(area1 = 116, area2 = 86, area3 = 43, n12 = 71, n23 = 43, n13 = 43,n123 = 43, 
                   #                   category = c("Interact with DNA", "Interact with histone"),
                   category = c("", "", ""),
                   fill = c("skyblue", "light coral", "orange"), alpha = rep(0.5, 3),  
                   # Circles
                   lwd = 4,
                   #                   lty = 'blank',
                   # Numbers
                   cex = 2.5,
                   fontface = "bold",
                   fontfamily = "sans",
                   cat.cex = 2,
                   cat.fontface = "bold",
                   cat.default.pos = "outer",
   #                cat.pos = c(-27, 153),
   #                cat.dist = c(0.04, 0),
                   cat.fontfamily = "sans",
                   #                   rotation = 1
                   
)
dev.off()