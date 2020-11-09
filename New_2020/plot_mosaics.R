# Plot mosaics
# Pregnancy data
require(dplyr)
require(magrittr)
require(mgcv)
require(caret)

outcomes3 <- read.csv("Outcomes.csv")

gL <- table(outcomes3[,c("Location","OutcomeGpNumeric")])
colnames(gL) <- c("Birth","Ongoing","Miscarriage","Ectopic")

pdf("mosaic_Loc_out.pdf",height = 4.75,width = 5.25)
mosaicplot(t(gL),main="",col=c('#fee5d9','#fcae91','#fb6a4a','#cb181d'),xlab="Pregnancy Outcome",ylab="Location",
           cex = 0.85,border="black",
           off=5,las=1)
dev.off()


gF <- table(outcomes3[,c("Fibrosis","OutcomeGpNumeric")])
colnames(gF) <- c("Birth","Ongoing","Miscarriage","Ectopic")

pdf("mosaic_fib_out.pdf",height = 4.75,width = 5.25)
mosaicplot(t(gF),main="",col=rev(c('#bf9b30','#f1b502','#ffc004','#ffe28a')),xlab="Pregnancy Outcome",ylab="Fibrosis",
           cex = 0.85,border="black",
           off=5,las=1)
dev.off()

gD <- table(outcomes3[,c("Diameter","OutcomeGpNumeric")])
colnames(gD) <- c("Birth","Ongoing","Miscarriage","Ectopic")
rownames(gD) <- c("Similar",'~Dissimilar',"Dissimilar")

pdf("mosaic_diam_out.pdf",height = 4.75,width = 5.25)
mosaicplot(t(gD),main="",col=c('#eff3ff','#bdd7e7','#6baed6','#2171b5'),xlab="Pregnancy Outcome",ylab="Diameter",
           cex = 0.85,border="black",
           off=5,las=1)
dev.off()

