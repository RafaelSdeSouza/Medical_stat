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
# Auxiliar function to randomly select a given column 


# Data-processing
preg <- read.csv("BTA-Patients-MAW.csv") %>% select(c("PREGNANT_NUMERIC",  "AGE",    
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
   filter(PREGNANT_NUMERIC %in% c(0,1)) %>% 
   mutate(PREGNANT_NUMERIC = as.numeric(as.character(PREGNANT_NUMERIC))) %>% 
   droplevels()  

# Sort left or right for each woman via bernoulli process
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

# Create new dataset with choosen features
preg2 <- preg %>%
  mutate(TL_rand = TL_rand) %>%
  mutate(TLD_rand = TLD_rand) %>%
  mutate(TLP_rand = TLP_rand) %>%
  mutate(ANAS_rand = ANAS_rand) %>%
  mutate(Fibr_rand = Fibr_rand) %>%
  mutate(PREGNANT_NUMERIC = as.factor(PREGNANT_NUMERIC)) %>%
  select(c("PREGNANT_NUMERIC","LIGATION_GROUP", "AGE", "TL_rand", "TLD_rand","TLP_rand","ANAS_rand","Fibr_rand"))


# Split train vs test sample
trainIndex <- createDataPartition(preg2$PREGNANT_NUMERIC, p = .9, 
                                  list = FALSE, 
                                  times = 1)
Train <- preg2[trainIndex,]
Test  <- preg2[-trainIndex,]


classif_gbm <- train(PREGNANT_NUMERIC~AGE+LIGATION_GROUP+TL_rand + TLD_rand + TLP_rand  + ANAS_rand + Fibr_rand, data = Train , method = "gbm",  tuneLength = 1)

classif_glm <- train(PREGNANT_NUMERIC~AGE+LIGATION_GROUP+TL_rand + TLD_rand + TLP_rand  + ANAS_rand + Fibr_rand, data = Train , method="glm", family="binomial")

classif_xgb<- train(PREGNANT_NUMERIC~AGE+LIGATION_GROUP+TL_rand + TLD_rand + TLP_rand  + ANAS_rand + Fibr_rand, data = Train , method="xgbTree")


p_fun <- function(object, newdata){predict(object, newdata=newdata, type="prob")[,2]}
yTest <- as.numeric(as.character(Test$PREGNANT_NUMERIC))

explainer_classif_gbm <- explain(classif_gbm, label = "gbm",
                                       data = Test[,-1], y = yTest,
                                 predict_function = p_fun)

explainer_classif_glm <- explain(classif_glm, label = "glm", 
                                        data = Test[,-1], y = yTest,
                                 predict_function = p_fun)


explainer_classif_xgb <- explain(classif_xgb, label = "xgb", 
                                        data = Test, y = yTest,
                                        predict_function = p_fun)

explain(HR_xgb_model,model_matrix_train,
        y=Train$PREGNANT_NUMERIC,
        label="modelxgb")   


mp_classif_gbm <- model_performance(explainer_classif_gbm)
mp_classif_xgb <- model_performance(explainer_classif_xgb)
mp_classif_glm <- model_performance(explainer_classif_glm)

plot(mp_classif_gbm, mp_classif_glm,mp_classif_xgb) +  my_style() 





vi_classif_gbm <- variable_importance(explainer_classif_gbm, loss_function = loss_cross_entropy,type = "difference")
vi_classif_glm <- variable_importance(explainer_classif_glm, loss_function = loss_root_mean_square,type = "difference")
vi_classif_xgb <- variable_importance(explainer_classif_xgb, loss_function = loss_root_mean_square,type = "difference")

plot(vi_classif_gbm, vi_classif_glm,vi_classif_xgb)





pdp_classif_xgb  <- variable_response(explainer_classif_xgb, variable =  "AGE", type = "ale")
pdp_classif_gbm  <- variable_response(explainer_classif_gbm, variable =  "AGE", type = "ale")
pdp_classif_glm  <- variable_response(explainer_classif_glm, variable =  "AGE", type = "ale")
plot(pdp_classif_xgb, pdp_classif_gbm, pdp_classif_glm ) +  my_style() 














model_matrix_train <- model.matrix(PREGNANT_NUMERIC~AGE+TL_rand + 
                                     TLD_rand + TLP_rand  + 
                                     ANAS_rand + Fibr_rand-1,Train)

model_matrix_test <- model.matrix(PREGNANT_NUMERIC~AGE+TL_rand + 
                                     TLD_rand + TLP_rand  + 
                                     ANAS_rand + Fibr_rand-1,Test)

data_train <- xgb.DMatrix(model_matrix_train,
                          label = Train$PREGNANT_NUMERIC)

data_test <- xgb.DMatrix(model_matrix_test,
                         label = Test$PREGNANT_NUMERIC)


param <- list(max_depth=2,objective="binary:logistic")

HR_xgb_model <- xgb.train(param,data_train,nrounds=50)

HR_glm_model <- glm(PREGNANT_NUMERIC~AGE+TL_rand + 
                      TLD_rand + TLP_rand  + 
                      ANAS_rand + Fibr_rand,Train,family=binomial(link = "logit"))





explainer_glm <- explain(HR_glm_model,Train)

expl_glm1 <- variable_response(explainer_glm,"AGE",
                              "pdp")
expl_glm2 <- variable_response(explainer_glm,"TL_rand",
                               "pdp")
expl_glm3 <- variable_response(explainer_glm,"TLD_rand",
                               "pdp")
expl_glm4 <- variable_response(explainer_glm,"TLP_rand",
                               "pdp")
expl_glm5 <- variable_response(explainer_glm,"ANAS_rand",
                               "pdp")
expl_glm6 <- variable_response(explainer_glm,"Fibr_rand",
                               "pdp")

expl_glm <- rbind(expl_glm1,expl_glm2,expl_glm3,expl_glm4,expl_glm5,
                  expl_glm6)


explainer_xgb <- explain(HR_xgb_model,model_matrix_train,
                         y=Train$PREGNANT_NUMERIC,
                         label="modelxgb")                        

expl_xgb1 <- variable_response(explainer_xgb,"AGE",
                              "pdp")
expl_xgb2 <- variable_response(explainer_xgb,"TL_rand",
                              "pdp")
expl_xgb3 <- variable_response(explainer_xgb,"TLD_rand",
                              "pdp")
expl_xgb4 <- variable_response(explainer_xgb,"TLP_rand",
                              "pdp")
expl_xgb5 <- variable_response(explainer_xgb,"ANAS_rand",
                              "pdp")
expl_xgb6 <- variable_response(explainer_xgb,"Fibr_rand",
                               "pdp")

expl_xgb <- rbind(expl_xgb1,expl_xgb2,expl_xgb3,expl_xgb4,expl_xgb5,
                  expl_xgb6)

expl <- rbind(expl_glm,expl_xgb) %>%
mutate(var = recode(var, TL_rand = "Tube length",
                    TLD_rand = "Tube length distal",
                    TLP_rand = "Tube length prox",
                    ANAS_rand = "ANASTOMOSIS",
                    Fibr_rand = "FIBROSIS")) %>%
  mutate(label = recode(label, lm = "glm",
                        xgb.Booster = "xgb"))


pdf("partial_dep.pdf",height = 7,width = 10)
ggplot(expl,aes(x=x,y=y,group=label,color=label)) +
  geom_line() + geom_point() + scale_colour_tableau()+
  facet_wrap(.~var,scales="free") + my_style() +
  ylab("Pregnancy probability") + xlab("")
dev.off()



vd_xgb <- variable_importance(explainer_xgb,type="difference")

plot(vd_xgb)

nobs <- model_matrix_train[1, , drop = FALSE]
pred_xgb <- prediction_breakdown(explainer_xgb,observation = nobs)
plot(pred_xgb )








