# Pregnancy data
require(dplyr)
require(magrittr)
require(mgcv)
require(caret)
library(pROC)
require(reshape)
require(randomForest)
library(vip)        # ML global interpretation
library(pdp)        # ML global interpretation
library(ggplot2)    # visualization pkg leveraged by above packages
library(h2o)        # ML model building
require(xgboost)
require(DALEX)
source("my_style.R")
require(ggthemes)
require(kernlab)
require(forcats)
# Auxiliar function to randomly select a given column 



outcomes <- read.csv("BTA-Pregnancies-anonymized.csv") %>% select(c("OutcomeGpNumeric",
                                                                    "LIGATION_GROUP",
                                                                    "AGE",    
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
  droplevels()   


# Sort left or right for each woman via bernoulli process
set.seed(42)
rlist <- rbinom(nrow(outcomes),1,0.5) + 1

temp1 <- outcomes[,c("LEFT_TUBE_LENGTH","RIGHT_TUBE_LENGTH")]
temp2 <- outcomes[,c( "TUBELENGTH_L_DISTAL",  "TUBELENGTH_R_DISTAL")]
temp3 <- outcomes[,c("TUBELENGTH_L_PROX", "TUBELENGTH_R_PROX")]
temp4 <- outcomes[,c("ANASTOMOSIS2_NUMERIC","ANASTOMOSIS1_NUMERIC")]
temp5 <- outcomes[,c("L_FIBROSIS_NUMERIC","R_FIBROSIS_NUMERIC")]
temp6 <- outcomes[,c("L_DIAMETER_NUMERIC",   "R_DIAMETER_NUMERIC")]

TL_rand <- c()
TLD_rand <- c()
TLP_rand <- c()
ANAS_rand <- c()
Fibr_rand <- c()
Diam_rand <- c()
for (i in 1:nrow(outcomes)) {
  TL_rand <-  append(TL_rand,temp1[i,rlist[i]])
  TLD_rand <- append(TLD_rand,temp2[i,rlist[i]])
  TLP_rand <- append(TLP_rand,temp3[i,rlist[i]])
  ANAS_rand <- append(ANAS_rand,temp4[i,rlist[i]])
  Fibr_rand <- append(Fibr_rand,temp5[i,rlist[i]])
  Diam_rand <- append(Diam_rand,temp6[i,rlist[i]])
}

# Create new dataset with choosen features
outcomes2 <- outcomes %>%
  mutate(TL_rand = TL_rand) %>%
  mutate(TLD_rand = TLD_rand) %>%
  mutate(TLP_rand = TLP_rand) %>%
  mutate(ANAS_rand = ANAS_rand) %>%
  mutate(Fibr_rand = Fibr_rand) %>%
  mutate(OutcomeGpNumeric = as.factor(OutcomeGpNumeric)) %>%
  mutate(Diam_rand  = Diam_rand) %>%
  select(c("OutcomeGpNumeric",
#"LIGATION_GROUP",
           "AGE", "TL_rand","ANAS_rand","Fibr_rand",
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

#colnames(outcomes2) <- c("OutcomeGpNumeric","Sterilization_Method", "Age", "Length","Location","Fibrosis",       
#                     "Diameter")

colnames(outcomes2) <- c("OutcomeGpNumeric", "Age", "Length","Location","Fibrosis",       
                     "Diameter")



gL <- table(outcomes2[,c("Location","OutcomeGpNumeric")])
colnames(gL) <- c("Birth","Ongoing","Miscarriage","Ectopic")

pdf("mosaic_Loc_out.pdf",height = 6,width = 7) 
mosaicplot(t(gL),main="",col=c('#fee5d9','#fcae91','#fb6a4a','#cb181d'),xlab="Pregnancy Outcome",ylab="Location",
           cex = 0.85,border="black",
           off=5,las=1)
dev.off()


gF <- table(outcomes2[,c("Fibrosis","OutcomeGpNumeric")])
colnames(gF) <- c("Birth","Ongoing","Miscarriage","Ectopic")

pdf("mosaic_fib_out.pdf",height = 6,width = 7) 
mosaicplot(t(gF),main="",col=rev(c('#bf9b30','#f1b502','#ffc004','#ffe28a')),xlab="Pregnancy Outcome",ylab="Fibrosis",
           cex = 0.85,border="black",
           off=5,las=1)
dev.off()

gD <- table(outcomes2[,c("Diameter","OutcomeGpNumeric")])
colnames(gD) <- c("Birth","Ongoing","Miscarriage","Ectopic")
rownames(gD) <- c("Similar",'~Dissimilar',"Dissimilar")

pdf("mosaic_diam_out.pdf",height = 6,width = 7) 
mosaicplot(t(gD),main="",col=c('#eff3ff','#bdd7e7','#6baed6','#2171b5'),xlab="Pregnancy Outcome",ylab="Diameter",
           cex = 0.85,border="black",
           off=5,las=1)
dev.off()


# Split train vs test sample
trainIndex <- createDataPartition(outcomes2$OutcomeGpNumeric, p = .95, 
                                  list = FALSE, 
                                  times = 1)
Train <- outcomes2[trainIndex,] 
#%>% mutate(LIGATION_GROUP = as.numeric(LIGATION_GROUP))


Test  <- outcomes2[-trainIndex,]  
#%>% mutate(LIGATION_GROUP = as.numeric(LIGATION_GROUP))



#### Model comparison 

classif_gam_A <- vgam(OutcomeGpNumeric~
                        s(Age,bs="cr",k=10)+
                        Sterilization_Method,multinomial, data = Train)

classif_gam_B <- vgam(OutcomeGpNumeric~
                        s(Age,bs="cr",k=10)+
                        s(Length,bs="cr",k=10) + 
                        Location + Fibrosis + Diameter,multinomial, data = Train)

classif_gam_C <- vgam(OutcomeGpNumeric~
                        s(Age,bs="cr",k=10)+
                        Sterilization_Method +
                        s(Length,bs="cr",k=10) + 
                        Location + Fibrosis + Diameter,multinomial, data = Train)



y_true = as.numeric(Train[,1])

LEA <- MultiLogLoss(y_pred = predict(classif_gam_A,type="response"), y_true= y_true)
LEB <- MultiLogLoss(y_pred = predict(classif_gam_B,type="response"), y_true= y_true)
LEC <- MultiLogLoss(y_pred = predict(classif_gam_C,type="response"), y_true= y_true)
####


classif_gam <- vgam(OutcomeGpNumeric~
       s(Age,bs="cr",k=10)+
       s(Length,bs="cr",k=10) + 
         Location + Fibrosis + Diameter,multinomial, data = Train)




#1=birth
#2=ongoing
#3=miscarriage
#4=ectopic

predgam1 <- function(m, x)   predict(m, x, type = "response")[,1]
predgam2 <- function(m, x)   predict(m, x, type = "response")[,2]
predgam3 <- function(m, x)   predict(m, x, type = "response")[,3]
predgam4 <- function(m, x)   predict(m, x, type = "response")[,4]


explain_gam_1 <- explain(classif_gam, data = Train[,-1], 
                        y = Train$OutcomeGpNumeric == "1", 
                        predict_function = predgam1, label = "Birth")
explain_gam_2 <- explain(classif_gam, data = Train[,-1], 
                        y = Train$OutcomeGpNumeric == "2",
                        predict_function = predgam2, label = "Ongoing")
explain_gam_3 <- explain(classif_gam, data = Train[,-1], 
                        y = Train$OutcomeGpNumeric == "3", 
                        predict_function = predgam3, label = "Miscarriage")

explain_gam_4 <- explain(classif_gam, data = Train[,-1], 
                        y = Train$OutcomeGpNumeric == "4",
                        predict_function = predgam4 , label = "Ectopic")



mp_gam_1 <- model_performance(explain_gam_1) %>% as.data.frame()
mp_gam_2 <- model_performance(explain_gam_2) %>% as.data.frame()
mp_gam_3 <- model_performance(explain_gam_3) %>% as.data.frame()
mp_gam_4 <- model_performance(explain_gam_4) %>% as.data.frame()

mp_gam_all <- rbind(mp_gam_1,mp_gam_2,mp_gam_3,mp_gam_4) %>%
  mutate(method="GAM")

mp_multi_all <- mp_gam_all %>% mutate(diff = abs(diff))
  



ggplot(mp_multi_all,aes(x=diff,group=label,
                  color=label)) +
  stat_ecdf(geom = "step",pad=FALSE) +
  my_style() +
  coord_flip()+
  scale_color_wsj() + ylab("ECDF") +
  scale_fill_wsj() +
  xlab("|Residuals|") +
  facet_wrap(.~method)




# AGE


ale_gam1_age   <- variable_response(explain_gam_1, variable =  "Age", type = "pdp")
ale_gam2_age   <- variable_response(explain_gam_2, variable =  "Age", type = "pdp")
ale_gam3_age   <- variable_response(explain_gam_3, variable =  "Age", type = "pdp")
ale_gam4_age   <- variable_response(explain_gam_4, variable =  "Age", type = "pdp")
ale_gam_all_age <- rbind(ale_gam1_age,ale_gam2_age,ale_gam3_age,ale_gam4_age ) %>%
  mutate(method = "GAM")


pdf("ale_age_outcome.pdf",height = 5.5,width = 6.5)
g <- ggplot(ale_gam_all_age,aes(x=x,y=y,group=label,color=label,fill=label)) + 
  geom_smooth(method = 'loess',span=0.3) +
  #+ geom_point() +
  my_style() +
  scale_linetype_stata(name="") +
  scale_color_wsj(name="") + ylab("Outcome likelihood") +
  scale_fill_wsj(name = "") +
  coord_cartesian(xlim=c(20,50),ylim=c(0.1,0.7)) + 
  xlab("Age (yr)") 
direct.label(g,"angled.endpoints")
dev.off()






## Tube lenght properties

ale_gam1_Length   <- variable_response(explain_gam_1, variable =  "Length", type = "pdp")
ale_gam2_Length   <- variable_response(explain_gam_2, variable =  "Length", type = "pdp")
ale_gam3_Length   <- variable_response(explain_gam_3, variable =  "Length", type = "pdp")
ale_gam4_Length   <- variable_response(explain_gam_4, variable =  "Length", type = "pdp")
ale_gam_all_Length <- rbind(ale_gam1_Length,ale_gam2_Length,ale_gam3_Length,ale_gam4_Length ) 






pdf("ale_TL_outcome.pdf",height = 5.5,width = 6.5)
g2 <- ggplot(ale_gam_all_Length,aes(x=x,y=y,group=label,color=label,fill=label)) +
  geom_smooth(method = 'loess',span=0.3) + 
 # geom_point() +  
  my_style() +
  scale_linetype_stata(name="") +
  scale_color_wsj(name="") + ylab("Outcome likelihood") +
  scale_fill_wsj(name = "") + 
  ylab("Outcome likelihood") + xlab("Length (cm)") +
  coord_cartesian(xlim=c(0.5,10))  
direct.label(g2,"lasso.labels")
dev.off()


# AGE


ale_gam1_fib   <- variable_response(explain_gam_1, variable =  "Fibrosis", type = "factor")
ale_gam2_fib   <- variable_response(explain_gam_2, variable =  "Fibrosis", type = "factor")
ale_gam3_fib   <- variable_response(explain_gam_3, variable =  "Fibrosis", type = "factor")
ale_gam4_fib   <- variable_response(explain_gam_4, variable =  "Fibrosis", type = "factor")
ale_gam_all_fib <- rbind(ale_gam1_fib,ale_gam2_fib,ale_gam3_fib,ale_gam4_fib ) %>%
  mutate(method = "GAM")






loss_entropy <- function(y, p){
  y <- as.numeric(y)
  -(1/length(y))*sum((y*log(p) + (1-y)*log(1-p))
  )
}

y_true = 
classif_gam$fitted.values



LS1 <- function(y,p){
  p <-  predict(classif_gam,type = "response")[,1]
  y <- as.numeric(Train[,1])
  y_pred <- p
  y_true <- y
  MultiLogLoss(y_pred, y_true)
}

vi_gam <- variable_importance(explain_gam_1,
                              n_sample = -1,loss_function = loss_entropy,
                              type = "difference")

vi_gam1 <- variable_importance(explain_gam_1,
                              n_sample = -1,loss_function = loss_entropy,
                              type = "difference")
vi_gam2 <- variable_importance(explain_gam_2,
                               n_sample = -1,loss_function = loss_entropy,
                               type = "difference")
vi_gam3 <- variable_importance(explain_gam_3,
                               n_sample = -1,loss_function = loss_entropy,
                               type = "difference")
vi_gam4 <- variable_importance(explain_gam_4,loss_function = loss_entropy,
                               n_sample = -1,
                               type = "difference")
pdf("vir_outcome.pdf",height = 10,width = 7.5)
plot(vi_gam1,vi_gam2,vi_gam3,vi_gam4)
dev.off()

ale_rf_age   <- variable_response(explain_rf, variable =  "AGE", type = "pdp")



ggplot(ale_rf_age,aes( x = x, y = y, group=label,color=label,linetype=label)) +  
  geom_smooth(method = 'loess',n = 500,show.legend=F) + geom_point() +  my_style() +
  scale_color_fivethirtyeight() + ylab("Pregnancy likelihood") +
  scale_fill_fivethirtyeight() +
  xlab("AGE (yr)") + theme(legend.position = "none") +
  coord_cartesian(xlim=c(21,50))



