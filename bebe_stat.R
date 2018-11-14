# Pregnancy data
require(dplyr)
require(magrittr)
require(mgcv)
require(visreg)
library(FactoMineR)
library("factoextra")

preg <- read.csv("BTA-Patients-MAW.csv") %>% select(c("BECOME_PREGNANT","TUBELENGTH_R_DISTAL","PREGNANT_NUMERIC",
                                                                 "TUBELENGTH_L_DISTAL","LIGATION_GROUP",
                                                                "AGE","RIGHT_TUBE_LENGTH","LEFT_TUBE_LENGTH",
                                                                "TUBELENGTH_R_PROX","TUBELENGTH_L_PROX","AV_TUBELENGTH_GP")) %>%
  filter(AGE != "Yes") %>%  mutate(AGE = as.numeric(as.character(AGE)) )  %>%
  filter(AGE > 10) %>%
  na.omit() %>% mutate(LEFT_TUBE_LENGTH = as.numeric(as.character(LEFT_TUBE_LENGTH)) ) %>%
  mutate(PCA_TUBE_LENGTH = -prcomp(data.frame(RIGHT_TUBE_LENGTH,LEFT_TUBE_LENGTH))$x[,1])




pca_dat <- data.frame(preg$LEFT_TUBE_LENGTH,preg$RIGHT_TUBE_LENGTH)
res.pca <- PCA(pca_dat, scale.unit = TRUE, ncp = 2, graph = TRUE)


fviz_pca_biplot(res.pca, 
                # Fill individuals by groups
                geom.ind = "point",
                pointshape = 21,
                pointsize = 2.5,
                fill.ind = preg$LIGATION_GROUP,
                col.ind = "black",
                repel = TRUE,        # Avoid label overplotting
                legend.title = list(fill = "Ligation Group")
) +
  ggpubr::fill_palette("jco")+      # Indiviual fill color
  ggpubr::color_palette("npg")      # Variable colors




fit <- gam(BECOME_PREGNANT~s(AGE,bs="cr",k=15)  + AV_TUBELENGTH_GP  + LIGATION_GROUP,data=preg,family= binomial(link="logit"),method="REML")


fit2 <- gam(BECOME_PREGNANT~s(AGE,bs="cr",k=15)  + s(PCA_TUBE_LENGTH,bs="cr",k=15)  + LIGATION_GROUP,data=preg,family= binomial(link="logit"),method="REML")




visreg(fit,"AV_TUBELENGTH_GP",by="LIGATION_GROUP",
       ylab = "BECOME PREGNANT", xlab="AV_TUBELENGTH_GP",scale="response")


visreg(fit2,"PCA_TUBE_LENGTH",by="LIGATION_GROUP",
       ylab = "BECOME PREGNANT", xlab="PC1",scale="response")

visreg(fit2,"AGE",by="LIGATION_GROUP",
       ylab = "BECOME PREGNANT", xlab="AGE",scale="response")


visreg2d(fit,"AGE","AV_TUBELENGTH_GP",
       ylab = "Average Tube length", xlab="AGE",zlab="Probability Pregnancy",scale="response",plot.type = "persp",
       phi=35,theta=47.5)

visreg2d(fit2,"AGE","PCA_TUBE_LENGTH",
         ylab = "PC1", xlab="AGE",zlab="Probability Pregnancy",scale="response",plot.type = "persp",
         phi=35,theta=47.5)
