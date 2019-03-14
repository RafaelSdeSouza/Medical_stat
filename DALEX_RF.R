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


# Population summary
aged <- cut(preg$AGE, breaks = c(20,30, 35, 40, 51))
agev <- data.frame(age = aged,LG=preg$LIGATION_GROUP) %>%
  mutate(LG = factor(LG,levels=c("Clip","Ring","Coagulation","Ligation/Resection"))) %>%

  na.omit()

pdf("Age.pdf",height = 5.5,width = 6.5)
ggplot(agev,aes(x=age,y = 100*(..count..)/sum(..count..),fill=LG)) +
  geom_bar() + my_style() +
  ylab("Per cent in each group") + xlab("Age group (yrs)") +
  scale_fill_wsj(name="") + theme(legend.spacing.x = unit(0.15, 'cm'))
dev.off()


# Anastomosis
PyrAN <-  preg[,c("ANASTOMOSIS2_NUMERIC","ANASTOMOSIS1_NUMERIC")]  %>%
  melt() %>%  mutate(variable = recode(variable, ANASTOMOSIS2_NUMERIC = "Left",
                                       ANASTOMOSIS1_NUMERIC = "Right")) %>%
  mutate(value = recode(value, "0" = "Identical",
                        "1" = "1-SPD",
                        "2" = "2-SPD",
                        "3" = "3-SPD")) %>% 
  mutate(class = "Anastomosis") 

pdf("Anastomosis.pdf",height = 5.5,width = 6.5)
ggplot(PyrAN, aes(x = value, y =  100*(..count..)/sum(..count..), fill = variable)) +   # Fill column
  geom_bar()   +  my_style() +
  scale_fill_wsj(name = "") + 
  xlab("Anastomosis segment position difference") +
  ylab(" Per cent in each group") + 
  theme(legend.spacing.x = unit(0.15, 'cm'))
dev.off()


# Diameter 
PyrDiam <-  preg[,c("L_DIAMETER_NUMERIC","R_DIAMETER_NUMERIC")]  %>%
  melt() %>%  mutate(variable = recode(variable, L_DIAMETER_NUMERIC = "Left",
                                       R_DIAMETER_NUMERIC = "Right")) %>%
              mutate(value = recode(value, "1" = "Similar",
                                           "2" = "Somewhat dissimilar",
                                           "3" = "Dissimilar")) %>% 
           mutate(class = "Diameter")  %>%
          mutate(value = factor(value,levels=c("Similar","Somewhat dissimilar","Dissimilar"))) 


pdf("diameter.pdf",height = 5.5,width = 6.5)
ggplot(PyrDiam, aes(x = value, y =  100*(..count..)/sum(..count..), fill = variable)) +   # Fill column
  geom_bar()   +  my_style() +
  scale_fill_wsj(name = "") + 
  xlab("Diameter") +
  ylab(" Per cent in each group") + 
  theme(legend.spacing.x = unit(0.15, 'cm'))
dev.off()


# Fibrosis
PyrFib <-  preg[,c("L_FIBROSIS_NUMERIC", "R_FIBROSIS_NUMERIC")]  %>%
  melt() %>%  mutate(variable = recode(variable, L_FIBROSIS_NUMERIC = "Left",
                                       R_FIBROSIS_NUMERIC = "Right")) %>%
  mutate(value = recode(value, "0" = "None",
                        "1" = "Mild",
                        "2" = "Moderate",
                        "3" = "Severe"))  %>% 
   mutate(class = "Fibrosis") %>% 
  mutate(value = factor(value,levels=c("None","Mild","Moderate","Severe"))) 



pdf("fibrosis.pdf",height = 5.5,width = 6.5)
ggplot(PyrFib, aes(x = value, y =  100*(..count..)/sum(..count..), fill = variable)) +   # Fill column
  geom_bar()   +  my_style() +
  scale_fill_wsj(name = "") + 
  xlab("Fibrosis") +
  ylab(" Per cent in each group") + 
  theme(legend.spacing.x = unit(0.15, 'cm'))
dev.off()





Py_all <- rbind(PyrAN,PyrDiam,PyrFib) %>% 
  mutate(value = factor(value,levels=c("Identical","1-SPD","2-SPD","3-SPD",
                                       "Similar","Somewhat dissimilar","Dissimilar","None","Mild","Moderate","Severe"))) 


pdf("anatomy.pdf",height = 5.5,width = 15)
ggplot(Py_all, aes(x = value, y =  100*(..count..)/sum(..count..), fill = class,alpha=variable)) +   # Fill column
  geom_bar(position = "dodge")   +  my_style() +
  scale_fill_wsj(name = "") +
  scale_alpha_manual(values=c(0.6,1)) +
  xlab("") +
  ylab(" Per cent in each group") + 
  theme(legend.position = "nome") +
  facet_wrap(.~class,scale="free_x")
dev.off()



outcomes <- read.csv("BTA-Pregnancies-anonymized.csv")


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
           #"LIGATION_GROUP",
           "AGE", "TL_rand", "TLD_rand","TLP_rand","ANAS_rand","Fibr_rand",
           "Diam_rand"))


# Split train vs test sample
trainIndex <- createDataPartition(preg2$PREGNANT_NUMERIC, p = .7, 
                                  list = FALSE, 
                                  times = 1)
Train <- preg2[trainIndex,] 
#%>% mutate(LIGATION_GROUP = as.numeric(LIGATION_GROUP))


Test  <- preg2[-trainIndex,]  
#%>% mutate(LIGATION_GROUP = as.numeric(LIGATION_GROUP))



#cor(Train[,3:8], use="complete.obs", method="kendall") 


# Train the models: GLM, GAM, RF

classif_glm <- glm(PREGNANT_NUMERIC~
                     #LIGATION_GROUP
                   AGE+TL_rand + TLD_rand + 
                  TLP_rand  + ANAS_rand + Fibr_rand + Diam_rand, data = Train, 
                   family=binomial(link = "logit"))


classif_rf <-  randomForest(PREGNANT_NUMERIC~
                              #LIGATION_GROUP+
                               AGE+
                               TL_rand + 
                              TLD_rand + 
  TLP_rand  + 
                               ANAS_rand + Fibr_rand + Diam_rand, 
                              data = Train, ntree=2000,nodesize=50)

varImpPlot(classif_rf)





p_rf <- function(object, newdata){predict(object, newdata=newdata, type="prob")[,2]}
p_glm  <- function(object, newdata){predict(object, newdata=newdata, type="response")}
yTest <- as.numeric(as.character(Test$PREGNANT_NUMERIC))




explain_glm  <- explain(classif_glm, label = "GLM", 
                                        data = Test[,-1], y = yTest,
                                 predict_function = p_glm)

explain_rf <- explain(classif_rf, label = "RF",
                      data = Test[,-1], y = yTest,
                      predict_function = p_rf)





# Model Performance


# residuals
mp_glm <- model_performance(explain_glm) %>% as.data.frame()
mp_rf <- model_performance(explain_rf) %>% as.data.frame()
mp_all <- rbind(mp_glm,mp_rf)


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


# AUC


pred_glm <- p_glm(classif_glm,newdata = Test[,-1])
pred_rf <- p_rf(classif_rf ,newdata = Test[,-1]) 


plot.roc(yTest,pred_glm, 
         percent=TRUE,  ci=TRUE, print.auc=TRUE,
         thresholds="best",
         print.thres="best")


plot.roc( yTest,pred_rf, 
         percent=TRUE,  ci=TRUE, print.auc=TRUE,
         thresholds="best",
         print.thres="best")



#  Acumulated Local Effects plot


# Age

ale_glm_age  <- variable_response(explain_glm , variable =  "AGE", type = "pdp")
ale_rf_age   <- variable_response(explain_rf, variable =  "AGE", type = "pdp")

ale_age_all <- rbind(ale_glm_age,ale_rf_age)


pdf("ale_age.pdf",height = 4,width = 5)
ggplot(ale_rf_age,aes( x = x, y = y, group=label,color=label,linetype=label)) +  
  geom_smooth(method = 'loess',n = 500,show.legend=F) + geom_point() +  my_style() +
 scale_color_fivethirtyeight() + ylab("Pregnancy likelihood") +
  scale_fill_fivethirtyeight() +
  xlab("AGE (yr)") + theme(legend.position = "none") +
  coord_cartesian(xlim=c(21,50))
dev.off()


## Tube lenght properties

ale_glm_TL1  <- variable_response(explain_glm , variable =  "TL_rand", type = "pdp")
ale_glm_TL2   <- variable_response(explain_glm , variable =  "TLD_rand", type = "pdp")
ale_glm_TL3  <- variable_response(explain_glm , variable =  "TLP_rand", type = "pdp")

ale_rf_TL1   <- variable_response(explain_rf, variable =  "TL_rand", type = "pdp")
ale_rf_TL2   <- variable_response(explain_rf, variable =  "TLD_rand", type = "pdp")
ale_rf_TL3   <- variable_response(explain_rf, variable =  "TLP_rand", type = "pdp")

ale_TL_all <- rbind(ale_glm_TL1,ale_glm_TL2,ale_glm_TL3,
                    ale_ann_TL1,ale_ann_TL2,ale_ann_TL3,
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
ale_rf_lg   <- variable_response(explain_rf, variable =  "LIGATION_GROUP", type = "factor")


pdf("ale_LG.pdf",height = 7.5,width = 8.75)
plot(ale_glm_lg,ale_rf_lg) + my_style()
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

vi_rf <- variable_importance(explain_rf,n_sample = -1,loss_function = lfu,type = "difference") 

`%not_in%` <- purrr::negate(`%in%`)


vi_all <- rbind(vi_glm,vi_ann,vi_rf) %>%
  mutate(dropout_loss = abs(dropout_loss)) %>%
  mutate(variable = recode(variable, TL_rand = "Tube length",
                      TLD_rand = "Tube length distal",
                      TLP_rand = "Tube length prox",
                      ANAS_rand = "ANASTOMOSIS",
                      Fibr_rand = "FIBROSIS"))  %>%
  filter(variable %not_in% c("_baseline_","_full_model_"))


pdf("vi.pdf",height = 6,width = 9)
ggplot(vi_all,aes(x=variable,y=dropout_loss,fill = label)) +
  coord_flip() +
  geom_bar(stat="identity") +
  scale_fill_fivethirtyeight() +
  my_style() +
  facet_wrap(.~label,scale="free_x")
 dev.off() 



svd_rf  <- single_variable(explainer_classif_xgb , variable = "LIGATION_GROUP", type = "factor")

plot(svd_rf) +  my_style() 









