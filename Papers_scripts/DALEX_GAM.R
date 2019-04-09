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
source("my_style.R")
require(ggthemes)
require(kernlab)
require(forcats)
require(kernlab)
require(directlabels)
require(MLmetrics)
# Auxiliar function to randomly select a given column 


# Data-processing
preg <- read.csv("BTA-Patients-MAW.csv") %>% select(c("PREGNANT_NUMERIC",  "AGE",    
                                                        "LIGATION_GROUP", 
                                                      "TUBELENGTH_L_DISTAL",  "TUBELENGTH_R_DISTAL", 
                                                      "LEFT_TUBE_LENGTH",     "RIGHT_TUBE_LENGTH",
                                                      "TUBELENGTH_L_PROX",    "TUBELENGTH_R_PROX",  
                                                       "L_DIAMETER_NUMERIC","R_DIAMETER_NUMERIC",
                                                      "L_DIAMETER_NUMERIC",   "R_DIAMETER_NUMERIC", 
                                                      "L_FIBROSIS_NUMERIC",   "R_FIBROSIS_NUMERIC",
                                                      "ANASTOMOSIS2_NUMERIC","ANASTOMOSIS1_NUMERIC"
                                                      )) %>%
   filter(AGE != "Yes") %>%  mutate(AGE = as.numeric(as.character(AGE)) )  %>%
   filter(AGE > 10) %>%
   na.omit() %>% mutate(LEFT_TUBE_LENGTH = as.numeric(as.character(LEFT_TUBE_LENGTH)) ) %>%
   mutate(TUBELENGTH_L_DISTAL = as.numeric(as.character(TUBELENGTH_L_DISTAL)) ) %>%
   filter(PREGNANT_NUMERIC %in% c(0,1)) %>% 
   mutate(PREGNANT_NUMERIC = as.numeric(as.character(PREGNANT_NUMERIC))) %>% 
   droplevels()  


# Sort left or right for each woman via bernoulli process
set.seed(42)
rlist <- rbinom(nrow(preg),1,0.5) + 1

temp1 <- preg[,c("LEFT_TUBE_LENGTH","RIGHT_TUBE_LENGTH")]
temp2 <- preg[,c( "TUBELENGTH_L_DISTAL",  "TUBELENGTH_R_DISTAL")]
temp3 <- preg[,c("TUBELENGTH_L_PROX", "TUBELENGTH_R_PROX")]
temp4 <- preg[,c("ANASTOMOSIS2_NUMERIC","ANASTOMOSIS1_NUMERIC")]
temp5 <- preg[,c("L_FIBROSIS_NUMERIC","R_FIBROSIS_NUMERIC")]
temp6 <- preg[,c("L_DIAMETER_NUMERIC",   "R_DIAMETER_NUMERIC")]

TL_rand <- c()
TLD_rand <- c()
TLP_rand <- c()
ANAS_rand <- c()
Fibr_rand <- c()
Diam_rand <- c()
for (i in 1:nrow(preg)) {
TL_rand <-  append(TL_rand,temp1[i,rlist[i]])
TLD_rand <- append(TLD_rand,temp2[i,rlist[i]])
TLP_rand <- append(TLP_rand,temp3[i,rlist[i]])
ANAS_rand <- append(ANAS_rand,temp4[i,rlist[i]])
Fibr_rand <- append(Fibr_rand,temp5[i,rlist[i]])
Diam_rand <- append(Diam_rand,temp6[i,rlist[i]])
}

# Create new dataset with choosen features
preg2 <- preg %>%
  mutate(TL_rand = TL_rand) %>%
  mutate(TLD_rand = TLD_rand) %>%
  mutate(TLP_rand = TLP_rand) %>%
  mutate(ANAS_rand = ANAS_rand) %>%
  mutate(Fibr_rand = Fibr_rand) %>%
  mutate(PREGNANT_NUMERIC = as.factor(PREGNANT_NUMERIC)) %>%
  mutate(Diam_rand  = Diam_rand) %>%
  select(c("PREGNANT_NUMERIC",
#           "LIGATION_GROUP",
           "AGE", "TL_rand", "ANAS_rand","Fibr_rand",
           "Diam_rand")) %>%
  mutate(Fibr_rand = recode(Fibr_rand, "0" = "None",
                      "1" = "Mild",
                      "2" = "Moderate",
                      "3" = "Severe")) %>%
  mutate(Fibr_rand = as.factor(Fibr_rand)) %>%
  mutate(Fibr_rand = factor(Fibr_rand,levels=c("None","Mild","Moderate","Severe"))) %>%
  mutate(ANAS_rand = recode(ANAS_rand, "0" = "Identical",
                      "1" = "1-SPD",
                      "2" = "2-SPD",
                      "3" = "3-SPD")) %>% 
  mutate(ANAS_rand = factor(ANAS_rand,levels=c("Identical","1-SPD","2-SPD","3-SPD")))  %>% 
  mutate(Diam_rand = recode(Diam_rand, "1" = "Similar","2" = "Somewhat dissimilar","3" = "Dissimilar")) %>% 
  mutate(Diam_rand = factor(Diam_rand,levels=c("Similar","Somewhat dissimilar","Dissimilar"))) 
 
 colnames(preg2) <- c("PREGNANT_NUMERIC", "Age", "Length","Location","Fibrosis",       
                   "Diameter")

 


ga <- table(preg2[,c("Location","PREGNANT_NUMERIC")])
colnames(ga) <- c("No","Yes")

gD <- table(preg2[,c("Diameter","PREGNANT_NUMERIC")])
colnames(gD) <- c("No","Yes")

gf <- table(preg2[,c("Fibrosis","PREGNANT_NUMERIC")])
colnames(gf) <- c("No","Yes")




pdf("mosaic_fib.pdf",height = 6,width = 7)
mosaicplot(t(gf),main="",col=c('#bf9b30','#f1b502','#ffc004','#ffe28a'),xlab="Pregnancy",ylab="Fibrosis",
           cex = 1,border="white",
           off=1,las=1)
dev.off()

pdf("mosaic_anas.pdf",height = 6,width = 7) 
mosaicplot(t(ga),main="",col=c('#fee5d9','#fcae91','#fb6a4a','#cb181d'),xlab="Pregnancy",
           ylab="Location",
           cex = 1,border="white",
           off=1,las=1)
dev.off()

pdf("mosaic_diam.pdf",height = 6,width = 7) 
mosaicplot(t(gD),main="",col=c('#eff3ff','#bdd7e7','#6baed6','#2171b5'),xlab="Pregnancy",
           ylab="Diameter",
           cex = 1,border="white",
           off=1,las=1)
dev.off()



Fg <- ggplot(preg2,aes(x=PREGNANT_NUMERIC,y = 100*(..count..)/sum(..count..),alpha=PREGNANT_NUMERIC,fill=Fibr_rand)) +
  geom_bar(position='dodge') + my_style() +
  scale_alpha_manual(guide="none",values=c(0.6,1)) +
  ylab("") + xlab("Fibrosis") + scale_x_discrete(labels=c("",""))+
  scale_fill_wsj(name="") + theme(legend.spacing.x = unit(0.15, 'cm')) 


Dg <- ggplot(preg2,aes(x=PREGNANT_NUMERIC,y = 100*(..count..)/sum(..count..),alpha=PREGNANT_NUMERIC,fill=Diam_rand)) +
  geom_bar(position='dodge') + my_style() +
  scale_alpha_manual(name="",guide="none",values=c(0.6,1)) +
  ylab("") + xlab("Diameter") + scale_x_discrete(labels=c("",""))+
  scale_fill_wsj(name="") + theme(legend.spacing.x = unit(0.15, 'cm'))


Ag <- ggplot(preg2,aes(x=PREGNANT_NUMERIC,alpha=PREGNANT_NUMERIC,y = 100*(..count..)/sum(..count..),fill=ANAS_rand)) +
  geom_bar(position='dodge') + my_style() +
  ylab("") + xlab("Location") +  scale_x_discrete(labels=c("","")) +
  scale_alpha_manual(guide="none",values=c(0.6,1)) +
  scale_fill_wsj(name="") + theme(legend.spacing.x = unit(0.15, 'cm'))

Tg <- ggplot(preg2,aes(x = TL_rand,  alpha=PREGNANT_NUMERIC, fill=PREGNANT_NUMERIC,group=PREGNANT_NUMERIC)) +
  geom_histogram(binwidth = 0.5,aes(group=PREGNANT_NUMERIC,y = 100*(..count..)/sum(..count..)))  +  my_style() +
  scale_fill_wsj(name = "") + scale_x_discrete(labels=c("",""))+
  scale_alpha_manual(guide="none",values=c(0.6,1)) +
  xlab("Length (cm)") + theme(legend.position = "none") + ylab("") 


pdf("preg_ana.pdf",height = 9.5,width = 11)
grid.arrange(Ag,Dg,Fg, Tg, ncol = 2,nrow=2,
             left =  text_grob("Percent in each group", size=18,rot=90) )
dev.off()


# Split train vs test sample
trainIndex <- createDataPartition(preg2$PREGNANT_NUMERIC, p = .95, 
                                  list = FALSE, 
                                  times = 1)
Train <- preg2[trainIndex,] 
#%>% mutate(LIGATION_GROUP = as.numeric(LIGATION_GROUP))


Test  <- preg2[-trainIndex,]  
#%>% mutate(LIGATION_GROUP = as.numeric(LIGATION_GROUP))




# Train the models: GLM, GAM, RF

classif_gam <- gam(PREGNANT_NUMERIC~
                     #LIGATION_GROUP
                   s(Age,bs="cr",k=12)+s(Length,bs="cr",k=12) +  Location + Fibrosis + Diameter, data = Train, 
                   family=binomial(link = "logit"))

####################
pdf("AgebyFibr.pdf",height = 5.5,width = 6.5)
  visreg(classif_gam,"Age", by = "Fibrosis",scale="response",
         aes(fill=Fibrosis),ylab = "Pregnancy probability",
         gg=TRUE,xlab="Age (yrs)",layout=c(2,2),
         line=list(col=c("black")))  + 
     facet_wrap(.~Fibrosis, ncol=2) +
     my_style() +  scale_fill_wsj() + scale_color_wsj() +
    theme(legend.position = "none") 
dev.off()

pdf("AgebyDiam.pdf",height = 5.5,width = 6.5)
visreg(classif_gam,"AGE", "Diam_rand",scale="response",
       ylab = "Pregnancy probability",
       gg=TRUE,xlab="Age (yrs)",layout=c(2,2),
       line=list(col=c("black")))  + 
  facet_wrap(.~Diam_rand, ncol=2) +
  my_style() +  scale_fill_wsj() + scale_color_wsj() +
  theme(legend.position = "none") 
dev.off()

pdf("AgebyAna.pdf",height = 5.5,width = 6.5)
visreg(classif_gam,"AGE", "ANAS_rand",scale="response",
       ylab = "Pregnancy probability",
       gg=TRUE,xlab="Age (yrs)",layout=c(2,2),
       line=list(col=c("black")))  + 
  facet_wrap(.~ANAS_rand, ncol=2) +
  my_style() +  scale_fill_wsj() + scale_color_wsj() +
  theme(legend.position = "none") 
dev.off()
##############################


##############################
pdf("TLbyFibr.pdf",height = 5.5,width = 6.5)
visreg(classif_gam,"TL_rand", by = "Fibr_rand",scale="response",
       aes(fill=Fibr_rand),ylab = "Pregnancy probability",
       gg=TRUE,xlab="Length (cm)",layout=c(2,2),
       line=list(col=c("black")))  + 
  facet_wrap(.~Fibr_rand, ncol=2) +
  my_style() +  scale_fill_wsj() + scale_color_wsj() +
  theme(legend.position = "none") + scale_x_continuous(breaks = c(0,3,6,9,12))
dev.off()

pdf("TLbyDiam.pdf",height = 5.5,width = 6.5)
visreg(classif_gam,"TL_rand", by = "Diam_rand",scale="response",
       aes(fill=Fibr_rand),ylab = "Pregnancy probability",
       gg=TRUE,xlab="Length (cm)",layout=c(2,2),
       line=list(col=c("black")))  + 
  facet_wrap(.~Diam_rand, ncol=2) +
  my_style() +  scale_fill_wsj() + scale_color_wsj() +
  theme(legend.position = "none") + scale_x_continuous(breaks = c(0,3,6,9,12))
dev.off()


pdf("TLbyAna.pdf",height = 5.5,width = 6.5)
visreg(classif_gam,"TL_rand", by = "ANAS_rand",scale="response",
       ylab = "Pregnancy probability",
       gg=TRUE,xlab="Length (cm)",layout=c(2,2),
       line=list(col=c("black")))  + 
  facet_wrap(.~ANAS_rand, ncol=2) +
  my_style() +  scale_fill_wsj() + scale_color_wsj() +
  theme(legend.position = "none") + scale_x_continuous(breaks = c(0,3,6,9,12))
dev.off()
##############################



p_gam  <- function(object, newdata){predict(object, newdata=newdata, type="response")}
yTest <- as.numeric(as.character(Test$PREGNANT_NUMERIC))




explain_gam  <- explain(classif_gam, label = "GAM", 
                                        data = Test[,-1], y = yTest,
                                 predict_function = p_gam)


# Model Performance


# residuals
mp_gam <- model_performance(explain_gam) %>% as.data.frame()
#mp_rf <- model_performance(explain_rf) %>% as.data.frame()
#mp_all <- rbind(mp_gam,mp_rf)


pdf("performance.pdf",height = 4,width = 5)
ggplot(mp_gam,aes(x=label,y=abs(diff),group=label,
                  fill=label)) +
  geom_boxplot() +
 my_style() +
  coord_flip()+
  scale_color_fivethirtyeight() + ylab("ECDF") +
  scale_fill_fivethirtyeight() +
  xlab("|Residuals|")
dev.off()


# Model comparison 
classif_gam_A <- gam(PREGNANT_NUMERIC~
                       LIGATION_GROUP +
                       s(AGE,bs="cr",k=12), data = Train, 
                     family=binomial(link = "logit"))

classif_gam_B <- gam(PREGNANT_NUMERIC~
                     #LIGATION_GROUP
                     s(AGE,bs="cr",k=12)+s(TL_rand,bs="cr",k=12) +  ANAS_rand + Fibr_rand + Diam_rand, data = Train, 
                   family=binomial(link = "logit"))

classif_gam_C <- gam(PREGNANT_NUMERIC~
                     LIGATION_GROUP +
                     s(AGE,bs="cr",k=12)+s(TL_rand,bs="cr",k=12) +  ANAS_rand + Fibr_rand + Diam_rand, data = Train, 
                   family=binomial(link = "logit"))




PRAUC(y_pred = classif_gam_A$fitted.values, y_true = Train[,1])
PRAUC(y_pred = classif_gam_B$fitted.values, y_true = Train[,1])
PRAUC(y_pred = classif_gam_C$fitted.values, y_true = Train[,1])

loss_entropy <- function(y, p){
  y <- as.numeric(y)
  -(1/length(y))*sum((y*log(p) + (1-y)*log(1-p))
  )
}

y_true = as.numeric(Train[,1])-1

LEA <- loss_entropy(y_true,yhat(classif_gam_A))
LEB <- loss_entropy(y_true,yhat(classif_gam_B))
LEC <- loss_entropy(y_true,yhat(classif_gam_C))



# AUC


pred_glm <- p_glm(classif_glm,newdata = Test[,-1])
#pred_rf <- p_rf(classif_rf ,newdata = Test[,-1]) 


plot.roc(yTest,pred_glm, 
         percent=TRUE,  ci=TRUE, print.auc=TRUE,
         thresholds="best",
         print.thres="best")


plot.roc(yTest,pred_rf, 
         percent=TRUE,  ci=TRUE, print.auc=TRUE,
         thresholds="best",
         print.thres="best")

# Variable Importance




vi_gam <- variable_importance(explain_gam,  n_sample = -1,loss_function = loss_entropy,type = "difference")

pdf("VIR_PregnancyLik.pdf",height = 4.5,width = 5.5)
plot(vi_gam,bar_width = 4)
dev.off()


`%not_in%` <- purrr::negate(`%in%`)


vi_all <- vi_gam %>%
filter(variable %not_in% c("_baseline_","_full_model_"))


pdf("vi.pdf",height = 6,width = 9)
ggplot(vi_all,aes(x=variable,y=dropout_loss,fill = label)) +
  coord_flip() +
  geom_bar(stat="identity") +
  scale_fill_fivethirtyeight() +
  my_style() +
  facet_wrap(.~label,scale="free_x")
 dev.off() 











