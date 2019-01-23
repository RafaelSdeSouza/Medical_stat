# Pregnancy data
require(dplyr)
require(magrittr)
require(mgcv)
require(visreg)
library(FactoMineR)
require(caret)
library("factoextra")
library(pROC)
require(PRROC)
library(lime) 
preg <- read.csv("BTA-Patients-MAW.csv") %>% select(c("BECOME_PREGNANT","TUBELENGTH_R_DISTAL","PREGNANT_NUMERIC",
                                                                 "TUBELENGTH_L_DISTAL","LIGATION_GROUP",
                                                                "AGE","RIGHT_TUBE_LENGTH","LEFT_TUBE_LENGTH",
                                                                "TUBELENGTH_R_PROX","TUBELENGTH_L_PROX","AV_TUBELENGTH_GP")) %>%
  filter(AGE != "Yes") %>%  mutate(AGE = as.numeric(as.character(AGE)) )  %>%
  filter(AGE > 10) %>%
  na.omit() %>% mutate(LEFT_TUBE_LENGTH = as.numeric(as.character(LEFT_TUBE_LENGTH)) ) %>%
  mutate(PCA_TUBE_LENGTH = -prcomp(data.frame(RIGHT_TUBE_LENGTH,LEFT_TUBE_LENGTH))$x[,1]) %>%
#  mutate(assymetry = (RIGHT_TUBE_LENGTH + LEFT_TUBE_LENGTH)/apply(data.frame(RIGHT_TUBE_LENGTH,LEFT_TUBE_LENGTH), 1, max)) %>%
  mutate(assymetry = (RIGHT_TUBE_LENGTH + LEFT_TUBE_LENGTH)) %>%
  filter(BECOME_PREGNANT %in% c("Yes","No")) %>% droplevels()


trainIndex <- createDataPartition(preg$BECOME_PREGNANT, p = .5, 
                                  list = FALSE, 
                                  times = 1)

Train <- preg[trainIndex,]
Test  <- preg[-trainIndex,]


# Case 0
fit0 <- gam(BECOME_PREGNANT~s(AGE,bs="cr",k=10)   + LIGATION_GROUP,data=preg,family = binomial(link="logit"))

pdf("case0.pdf",height = 5,width = 6.5)
visreg(fit0,"AGE",by="LIGATION_GROUP",
       ylab = "Pregnancy probability", xlab="Age",scale="response")
dev.off()

# Case 1

fit <- gam(BECOME_PREGNANT~s(AGE,bs="cr",k=10)  + assymetry  + LIGATION_GROUP,data=preg,family= binomial(link="logit"))

pdf("case1_0.pdf",height = 5,width = 6.5)
visreg(fit,"assymetry",by="LIGATION_GROUP",cond = list(AGE = 25),
       ylab = "Pregnancy probability", xlab="(L + R)/max(L,R)",scale="response")

pdf("case1_0.pdf",height = 10,width = 12)
par(mfrow=c(2,2))
visreg2d(fit,"AGE","assymetry",cond = list(LIGATION_GROUP = c("Clip")),plot.type = "persp",scale="response",
         zlab = "Pregnancy probability", xlab="AGE",ylab="L/R asymmetry (cm)",main="Clip",
         theta=30,phi=13.5,color="#e41a1c")
visreg2d(fit,"AGE","assymetry",cond = list(LIGATION_GROUP = c("Coagulation")),plot.type = "persp",scale="response",
         zlab = "Pregnancy probability", xlab="AGE",ylab="Average tube lenght (cm)",main="Coagulation",
         theta=30,phi=13.5,color="#e41a1c")
visreg2d(fit,"AGE","assymetry",cond = list(LIGATION_GROUP = c("Ligation/Resection")),plot.type = "persp",scale="response",
         zlab = "Pregnancy probability", xlab="AGE",ylab="Average tube lenght (cm)",main="Ligation/Resection",
         theta=30,phi=13.5,color="#e41a1c")
visreg2d(fit,"AGE","assymetry",cond = list(LIGATION_GROUP = c("Ring")),plot.type = "persp",scale="response",
         zlab = "Pregnancy probability", xlab="AGE",ylab="Average tube lenght (cm)",main="Ring",
         theta=30,phi=13.5,color="#e41a1c")
dev.off()




pdf("case1_1.pdf",height = 5,width = 6.5)
visreg(fit,"AV_TUBELENGTH_GP",by="LIGATION_GROUP",
       ylab = "Pregnancy probability", xlab="Average tube lenght (cm)",scale="response")
dev.off()

fit2 <- gam(BECOME_PREGNANT~s(AGE,bs="cr",k=15)  + s(PCA_TUBE_LENGTH,bs="cr",k=15)  + LIGATION_GROUP,data=Train,family= binomial(link="logit"))


x1=anova(fit2)$chi.sq-anova(fit2)$edf
x2=anova(Beta_GAM)$chi.sq-anova(Beta_GAM)$edf
x3=x1+x2


pred0 <- predict(fit0,newdata = Test[,-1],type="response")
pred1 <- predict(fit,newdata = Test[,-1],type="response")

pred2 <- predict(fit2,newdata = Test[,-1],type="response")




plot(roc(Test$BECOME_PREGNANT, pred0, direction="<"),print.auc=TRUE,
     lwd=3)


plot(roc(Test$BECOME_PREGNANT, pred1, direction="<"),print.auc=TRUE,
      lwd=3)



plot.roc(Test$BECOME_PREGNANT,pred0, 
         percent=TRUE,  ci=TRUE, print.auc=TRUE,
         thresholds="best",
         print.thres="best")


plot.roc(Test$BECOME_PREGNANT,pred2, 
         percent=TRUE,  ci=TRUE, print.auc=TRUE,
           thresholds="best",
          print.thres="best")
         


pr <- pr.curve(scores.class0 = Test$BECOME_PREGNANT, weights.class0 = pred2,curve = TRUE)


fitglm0 <- glm(BECOME_PREGNANT~ AGE   + LIGATION_GROUP,data=Train,family = binomial(link="logit"))
predglm0 <- predict(fitglm0,newdata = Test[,-1],type="response")

PRAUC(y_pred = as.numeric(predglm0), y_true = Test$BECOME_PREGNANT)

ggplot(data.frame(pr$curve),aes(x=X1,y=X2,color=X3)) + geom_line() 





pdf("case1.pdf",height = 6,width = 6)
visreg(fit,"AV_TUBELENGTH_GP",by="LIGATION_GROUP",
       ylab = "BECOME PREGNANT", xlab="AV_TUBELENGTH_GP",scale="response")
dev.off()

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



rocdf <- function(pred, obs, data=NULL, type=NULL) {
  # plot_type is "roc" or "pr"
  if (!is.null(data)) {
    pred <- eval(substitute(pred), envir=data)
    obs  <- eval(substitute(obs), envir=data)
  }
  
  rocr_xy <- switch(type, roc=c("tpr", "fpr"), pr=c("prec", "rec"))
  rocr_df <- prediction(pred, obs)
  rocr_pr <- performance(rocr_df, rocr_xy[1], rocr_xy[2])
  xy <- data.frame(rocr_pr@x.values[[1]], rocr_pr@y.values[[1]])
  
  # If PR, designate first (undefined) point as recall = 0, precision = x
  if (type=="pr") {
    xy[1, 2] <- 0
    #xy <- xy[!(rowSums(xy)==0), ]
  }
  
  colnames(xy) <- switch(type, roc=c("tpr", "fpr"), pr=c("rec", "prec"))
  return(xy)
}

