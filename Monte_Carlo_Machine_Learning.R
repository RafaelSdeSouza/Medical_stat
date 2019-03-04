# Pregnancy data
require(dplyr)
require(magrittr)
require(mgcv)
require(visreg)
require(caret)
library(pROC)
library(lime) 
require(reshape)
require(corrplot)
require(randomForest)
library(lime)       # ML local interpretation
library(vip)        # ML global interpretation
library(pdp)        # ML global interpretation
library(ggplot2)    # visualization pkg leveraged by above packages
library(caret)      # ML model building
library(h2o)        # ML model building
require(xgboost)
require(DALEX)
# Auxiliar function to randomly select a given column 



crand <- function(x,seed){
  if (!missing(seed)) 
    set.seed(seed) 
  xr <- x[1+rbinom(1,1,0.5)]
  return(xr)
}

#"PREGNANT_NUMERIC", 
preg <- read.csv("BTA-Patients-MAW.csv") %>% select(c("BECOME_PREGNANT",  "AGE",    
                                                        "LIGATION_GROUP", 
                                                      "TUBELENGTH_L_DISTAL",  "TUBELENGTH_R_DISTAL", 
                                                      "LEFT_TUBE_LENGTH",     "RIGHT_TUBE_LENGTH",
                                                      "TUBELENGTH_L_PROX",    "TUBELENGTH_R_PROX",        
                                                      "L_DIAMETER_NUMERIC",   "R_DIAMETER_NUMERIC", 
                                                      "L_FIBROSIS_NUMERIC",   "R_FIBROSIS_NUMERIC",
                                                      "ANASTOMOSIS2_NUMERIC","ANASTOMOSIS1_NUMERIC"
                                                      )) %>%
   filter(AGE != "Yes") %>%  mutate(AGE = as.numeric(as.character(AGE)) )  %>%
   filter(AGE > 10) %>%
   na.omit() %>% mutate(LEFT_TUBE_LENGTH = as.numeric(as.character(LEFT_TUBE_LENGTH)) ) %>%
   mutate(TUBELENGTH_L_DISTAL = as.numeric(as.character(TUBELENGTH_L_DISTAL)) ) %>%
   filter(BECOME_PREGNANT %in% c("Yes","No")) %>% droplevels()  


rlist <- rbinom(nrow(preg),1,0.5) + 1

temp1 <- preg[,c("LEFT_TUBE_LENGTH","RIGHT_TUBE_LENGTH")]
temp2 <- preg[,c( "TUBELENGTH_L_DISTAL",  "TUBELENGTH_R_DISTAL")]
temp3 <- preg[,c("TUBELENGTH_L_PROX", "TUBELENGTH_R_PROX")]
temp4 <- preg[,c("ANASTOMOSIS2_NUMERIC","ANASTOMOSIS1_NUMERIC")]
temp5 <- preg[,c("L_FIBROSIS_NUMERIC","R_FIBROSIS_NUMERIC")]



TL_rand <- c()
TLD_rand <- c()
TLP_rand <- c()
ANAS_rand <- c()
Fibr_rand <- c()
for (i in 1:nrow(preg)) {
TL_rand <-  append(TL_rand,temp1[i,rlist[i]])
TLD_rand <- append(TLD_rand,temp2[i,rlist[i]])
TLP_rand <- append(TLP_rand,temp3[i,rlist[i]])
ANAS_rand <- append(ANAS_rand,temp4[i,rlist[i]])
Fibr_rand <- append(Fibr_rand,temp5[i,rlist[i]])
}

preg2 <- preg %>%
  mutate(TL_rand = TL_rand) %>%
  mutate(TLD_rand = TLD_rand) %>%
  mutate(TLP_rand = TLP_rand) %>%
  mutate(ANAS_rand = ANAS_rand) %>%
  mutate(Fibr_rand = Fibr_rand) 



  
  







M <- preg2[,c("TL_rand","TLD_rand","TLP_rand","ANAS_rand","Fibr_rand")]

corrplot(cor(M),method="number",type="lower")


trainIndex <- createDataPartition(preg$BECOME_PREGNANT, p = .5, 
                                  list = FALSE, 
                                  times = 1)
Train <- preg2[trainIndex,]
Test  <- preg2[-trainIndex,]



fit.caret <- train(
  BECOME_PREGNANT~AGE+TL_rand + TLD_rand + TLP_rand  + ANAS_rand + Fibr_rand, 
  data = Train, 
  method = 'ranger',
  trControl = trainControl(method = "cv", number = 5, classProbs = TRUE),
  tuneLength = 1,
  importance = 'impurity'
)

vip(fit.caret) + ggtitle("ranger: RF")



global_obs <- Train %>% select(c( "AGE",
                                  "TL_rand","TLD_rand","TLP_rand","ANAS_rand","Fibr_rand"
)) 
explainer_caret <- lime(global_obs, fit.caret, n_bins = 5)

class(explainer_caret)

local_obs <- Train %>% select(c( "AGE",
                                 "TL_rand","TLD_rand","TLP_rand","ANAS_rand","Fibr_rand"
))  %>%  filter(AGE >= 47)

explanation_caret <- explain(
  x = local_obs, 
  explainer = explainer_caret, 
  n_permutations = 5000,
  dist_fun = "gower",
  kernel_width = .75,
  n_features = 10, 
  feature_select = "highest_weights",
  labels = "Yes"
)

plot_features(explanation_caret)

# Case 1 Age + ligation method

fit0 <- gam(BECOME_PREGNANT~s(AGE,bs="cr",k=10)   + LIGATION_GROUP,data=Train,family = binomial(link="logit"))


# Case 2 Age + Fallopian tube anatomy 


fit <- gam(BECOME_PREGNANT~s(AGE,bs="cr",k=10)  + TL_rand + TLD_rand + TLP_rand  + ANAS_rand + Fibr_rand,data=Train,family= binomial(link="logit"))


pred0 <- predict(fit0,newdata = Test[,-1],type="response")
pred1 <- predict(fit,newdata = Test[,-1],type="response")


plot(roc(Test$BECOME_PREGNANT, pred0, direction="<"),print.auc=TRUE,
     lwd=3)


plot(roc(Test$BECOME_PREGNANT, pred1, direction="<"),print.auc=TRUE,
     lwd=3)



tree_model <- randomForest(BECOME_PREGNANT~AGE+TL_rand + TLD_rand + TLP_rand  + ANAS_rand + Fibr_rand,data=Train,ntree=5000)
predRF <- predict(tree_model,newdata = Test[,-1],type="prob")
plot(roc(Test$BECOME_PREGNANT, predRF[,2], direction="<"),print.auc=TRUE,
     lwd=3)


pred2 <- predict(fit2,newdata = Test[,-1],type="response")








tubesum <- preg %>%  select(c("LIGATION_GROUP","TUBELENGTH_R_DISTAL",
                         "TUBELENGTH_L_DISTAL",
                         "RIGHT_TUBE_LENGTH","LEFT_TUBE_LENGTH",
                         "TUBELENGTH_R_PROX","TUBELENGTH_L_PROX",
                         "R_DIAMETER_NUMERIC",
                         "L_DIAMETER_NUMERIC",
                         "R_FIBROSIS_NUMERIC",
                         "L_FIBROSIS_NUMERIC"
                           )) %>% 
                           melt(.,id.vars="LIGATION_GROUP") %>%
                           mutate(value=as.numeric(value))



ggplot(tubesum,aes(x=variable,y=value,group=variable,
                   fill=variable)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(.~LIGATION_GROUP) +
  ylab("Lenght") + xlab("") +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        text=element_text(family="serif"),
        strip.text = element_text(size=10),
        strip.background = element_rect(fill = "white"),
        axis.title = element_text(size=10),
        axis.text  = element_text(size=7),
        axis.ticks = element_line(size = 0.45),
        axis.text.x = element_text(angle = 60,hjust = 1),
        legend.position = "none",
        legend.background = element_rect(colour = "white", fill = "white"),
        legend.text = element_text(size=25),
        legend.text.align = 0,
        legend.key = element_rect(colour = "white", fill = "white")) 


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

fit <- gam(BECOME_PREGNANT~s(AGE,bs="cr",k=10)  + TL + DL + FL  + LIGATION_GROUP,data=preg,family= binomial(link="logit"))

pdf("case1_0.pdf",height = 5,width = 6.5)
visreg(fit,"assymetry",by="LIGATION_GROUP",cond = list(AGE = 25),
       ylab = "Pregnancy probability", xlab="L + R",scale="response")

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

