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



#cor(Train[,3:8], use="complete.obs", method="kendall") 


# Train the models: GLM, GAM, RF

classif_glm <- glm(PREGNANT_NUMERIC~AGE+LIGATION_GROUP+TL_rand + TLD_rand + 
                   TLP_rand  + ANAS_rand + Fibr_rand, data = Train, 
                   family=binomial(link = "logit"))


classif_svm <- ksvm(PREGNANT_NUMERIC~AGE+LIGATION_GROUP+TL_rand + TLD_rand + TLP_rand  + 
                    ANAS_rand + Fibr_rand,data = Train,
                    kernel="rbfdot",type="C-bsvc",prob.model = TRUE)


classif_rf <-  randomForest(PREGNANT_NUMERIC~AGE+LIGATION_GROUP+TL_rand + 
                            TLD_rand + TLP_rand  + ANAS_rand + Fibr_rand, 
                            data = Train, ntree=2500)



p_rf <- function(object, newdata){predict(object, newdata=newdata, type="prob")[,2]}
p_glm  <- function(object, newdata){predict(object, newdata=newdata, type="response")}
yTest <- as.numeric(as.character(Test$PREGNANT_NUMERIC))




explain_glm  <- explain(classif_glm, label = "GLM", 
                                        data = Test[,-1], y = yTest,
                                 predict_function = p_glm)

explain_svm  <- explain(classif_svm, label = "SVM", 
                                        data = Test[,-1], y = yTest,
                                        predict_function = p_rf)


explain_rf <- explain(classif_rf, label = "RF",
                      data = Test[,-1], y = yTest,
                      predict_function = p_rf)




# Model Performance



mp_glm <- model_performance(explain_glm) %>% as.data.frame()
mp_svm  <- model_performance(explain_svm) %>% as.data.frame()
mp_rf <- model_performance(explain_rf) %>% as.data.frame()
mp_all <- rbind(mp_glm,mp_svm,mp_rf)



pdf("performance.pdf",height = 4,width = 5)
ggplot(mp_all,aes(x=label,y=abs(diff),group=label,
                  fill=label)) +
  geom_boxplot() +
 my_style() +
  coord_flip()+
  scale_color_fivethirtyeight() + ylab("ECDF") +
  scale_fill_fivethirtyeight() +
  xlab("|Residuals|")
dev.off()

#  Acumulated Local Effects plot


# Age

ale_glm_age  <- variable_response(explain_glm , variable =  "AGE", type = "pdp")
ale_svm_age  <- variable_response(explain_svm, variable =  "AGE", type = "pdp")
ale_rf_age   <- variable_response(explain_rf, variable =  "AGE", type = "pdp")

ale_age_all <- rbind(ale_glm_age ,ale_svm_age,ale_rf_age  )


pdf("ale_age.pdf",height = 4,width = 5)
ggplot(ale_age_all,aes( x = x, y = y, group=label,color=label,linetype=label)) +  
  geom_line() + geom_point() +  my_style() +
 scale_color_fivethirtyeight() + ylab("Pregnancy likelihood") +
  scale_fill_fivethirtyeight() +
  xlab("AGE (yr)")
dev.off()


## Tube lenght properties

ale_glm_TL1  <- variable_response(explain_glm , variable =  "TL_rand", type = "pdp")
ale_glm_TL2   <- variable_response(explain_glm , variable =  "TLD_rand", type = "pdp")
ale_glm_TL3  <- variable_response(explain_glm , variable =  "TLP_rand", type = "pdp")

ale_svm_TL1   <- variable_response(explain_svm, variable =  "TL_rand", type = "pdp")
ale_svm_TL2   <- variable_response(explain_svm, variable =  "TLD_rand", type = "pdp")
ale_svm_TL3  <- variable_response(explain_svm, variable =  "TLP_rand", type = "pdp")


ale_rf_TL1   <- variable_response(explain_rf, variable =  "TL_rand", type = "pdp")
ale_rf_TL2   <- variable_response(explain_rf, variable =  "TLD_rand", type = "pdp")
ale_rf_TL3   <- variable_response(explain_rf, variable =  "TLP_rand", type = "pdp")

ale_TL_all <- rbind(ale_glm_TL1,ale_glm_TL2,ale_glm_TL3,
                    ale_svm_TL1,ale_svm_TL2,ale_svm_TL3,
                    ale_rf_TL1,ale_rf_TL2,ale_rf_TL3 ) %>%
  mutate(var = recode(var, TL_rand = "Tube length",
                      TLD_rand = "Tube length distal",
                      TLP_rand = "Tube length prox"))

pdf("ale_TL.pdf",height = 4.5,width = 13)
ggplot(ale_TL_all,aes(x=x,y=y,group=label,color=label)) +
  geom_line() + geom_point() +  my_style() +
  scale_color_fivethirtyeight() +
  facet_wrap(.~var) + my_style() +
  ylab("Pregnancy likelihood") + xlab("Length (cm)") +
  coord_cartesian(xlim=c(0,11.9))

dev.off()


## Ligation Group

ale_glm_lg  <- variable_response(explain_glm , variable =  "LIGATION_GROUP", type = "factor")
ale_svm_lg  <- variable_response(explain_svm, variable =  "LIGATION_GROUP", type = "factor")
ale_rf_lg   <- variable_response(explain_rf, variable =  "LIGATION_GROUP", type = "factor")


pdf("ale_LG.pdf",height = 7.5,width = 8.75)
plot(ale_glm_lg,ale_svm_lg,ale_rf_lg) + my_style()
dev.off()

## Anamastosis




## FIbrosis






# Variable Importance
lfu = function(y, p){
  sum(
    -(y*log(p) + (1-y)*log(1-p))
    )
}
lfu()

log(1-p_rf(classif_rf))
lfu(as.numeric(as.character(Train$PREGNANT_NUMERIC)),p_rf(classif_rf))

vi_glm <- variable_importance(explain_glm,  n_sample = -1,loss_function = lfu,type = "difference")

vi_svm <- variable_importance(explain_svm,  n_sample = -1,loss_function = lfu,type = "difference") 

vi_rf <- variable_importance(explain_rf,n_sample = -1,loss_function = lfu,type = "difference") 

`%not_in%` <- purrr::negate(`%in%`)


vi_all <- rbind(vi_glm,vi_svm,vi_rf) %>%
  mutate(dropout_loss = abs(dropout_loss)) %>%
  mutate(variable = recode(variable, TL_rand = "Tube length",
                      TLD_rand = "Tube length distal",
                      TLP_rand = "Tube length prox",
                      ANAS_rand = "ANASTOMOSIS",
                      Fibr_rand = "FIBROSIS"))  %>%
  filter(variable %not_in% c("_baseline_","_full_model_"))


ggplot(vi_all,aes(x=variable,y=dropout_loss,fill = label)) +
  coord_flip() +
  geom_bar(stat="identity") +
  scale_fill_fivethirtyeight() +
  my_style() +
  facet_wrap(.~label)
  










ggplot(vi_classif_gbm, aes(x=variable, y=dropout_loss)) +
  geom_segment( aes(x=variable, xend=variable, y=1, yend=dropout_loss), color="skyblue") +
  geom_point( color="blue", size=4, alpha=0.6) +
  theme_light() +
  coord_flip() +
  my_style() 


pdp_classif_xgb  <- variable_response(explainer_classif_xgb, variable =  "TLP_rand", type = "ale")
pdp_classif_gbm  <- variable_response(explainer_classif_gbm, variable =  "TLP_rand", type = "ale")
pdp_classif_glm  <- variable_response(explainer_classif_glm, variable = "TLP_rand", type = "ale")
plot(pdp_classif_xgb, pdp_classif_gbm, pdp_classif_glm ) +  my_style() 





svd_rf  <- single_variable(explainer_classif_xgb , variable = "LIGATION_GROUP", type = "factor")

plot(svd_rf) +  my_style() 






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








