# Prep data
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
require(xgboost)
require(DALEX)
source("my_style.R")
require(ggthemes)
require(kernlab)
require(forcats)
require(VGAM)

Train0 <- read.csv("Outcomes.csv")
#%>%
# mutate(OutcomeGpNumeric = recode(OutcomeGpNumeric,Ectopic = "Miscarriage")) %>%
#  filter(OutcomeGpNumeric != "Ongoing") %>% mutate(OutcomeGpNumeric=
 #droplevels(OutcomeGpNumeric))

# Split train vs test sample
trainIndex <- createDataPartition(Train0$OutcomeGpNumeric, p = .75,
                                  list = FALSE,
                                  times = 1)

Train <- Train0[ trainIndex,]
Test  <- Train0[-trainIndex,]


fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 5,
  ## Estimate class probabilities
  repeats = 1)






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

#### Model comparison
pred_A <- predict(classif_gam_A ,newdata=Test[,-1], type = "response")
pclass_A <- max.col(pred_A) %>% as.factor()
pclass_A <- recode(pclass_A, "1" = "Birth","2" = "Ectopic","3"="Miscarriage","4" = "Ongoing")
table(pclass_A,Test[,1])


pred_B <- predict(classif_gam_B ,newdata=Test[,-1], type = "response")
pclass_B <- max.col(pred_B) %>% as.factor()
pclass_B <- recode(pclass_B, "1" = "Birth","2" = "Ectopic","3"="Miscarriage","4" = "Ongoing")
table(pclass_B,Test[,1])


pred_C <- predict(classif_gam_C ,newdata=Test[,-1], type = "response")
pclass_C <- max.col(pred_C) %>% as.factor()
pclass_C <- recode(pclass_C, "1" = "Birth","2" = "Ectopic","3"="Miscarriage","4" = "Ongoing")
table(pclass_C,Test[,1])



confusionMatrix(pclass_C,
                Test[,1],
                mode = "everything")






y_true = as.numeric(Train[,1])

LEA <- MLmetrics::MultiLogLoss(y_pred = predict(classif_gam_A,type="response"), y_true= y_true)
LEB <- MLmetrics::MultiLogLoss(y_pred = predict(classif_gam_B,type="response"), y_true= y_true)
LEC <- MLmetrics::MultiLogLoss(y_pred = predict(classif_gam_C,type="response"), y_true= y_true)
####


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


# AGE


ale_gam1_age   <- variable_effect(explain_gam_1, variable =  "Age", type = "partial_dependency")
ale_gam2_age   <- variable_effect(explain_gam_2, variable =  "Age", type = "partial_dependency")
ale_gam3_age   <- variable_effect(explain_gam_3, variable =  "Age", type = "partial_dependency")
ale_gam4_age   <- variable_effect(explain_gam_4, variable =  "Age", type = "partial_dependency")
ale_gam_all_age <- rbind(ale_gam1_age,ale_gam2_age,ale_gam3_age,ale_gam4_age ) %>%
  mutate(method = "GAM")




pdf("ale_age_outcome.pdf",height = 5.5,width = 6.5)
g <- ggplot(ale_gam_all_age,aes(x=`_x_`,y=`_yhat_`,group=`_label_`,
                                color=`_label_`,fill=`_label_`)) +
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

#y_true = classif_gam$fitted.values



logit <- function(x) exp(x)/(1+exp(x))

LS1 <- function(y,p){
  p <-  predict(classif_gam,type = "response")[,1]
  y <- as.numeric(Train[,1])
  y_pred <- p
  y_true <- y
  MultiLogLoss(y_pred, y_true)
}

vi_gam <- variable_importance(explain_gam_1,
                              type = "raw")

vi_gam1 <- variable_importance(explain_gam_1,type = "raw",loss_function = function(observed, predicted)
  sum((observed - logit(predicted))^2))
vi_gam2 <- variable_importance(explain_gam_2,type = "raw",loss_function = function(observed, predicted)
  sum((observed - logit(predicted))^2))
vi_gam3 <- variable_importance(explain_gam_3,type = "raw",loss_function = function(observed, predicted)
  sum((observed - logit(predicted))^2))
vi_gam4 <- variable_importance(explain_gam_4,type = "raw",
                               loss_function = function(observed, predicted)
                                 sum((observed - logit(predicted))^2))



pdf("vir_outcome1.pdf",height = 3.5,width = 4)
plot(vi_gam1)
dev.off()

pdf("vir_outcome2.pdf",height = 3.5,width = 4)
plot(vi_gam2)
dev.off()

pdf("vir_outcome3.pdf",height = 3.5,width = 4)
plot(vi_gam3)
dev.off()

pdf("vir_outcome4.pdf",height = 3.5,width = 4)
plot(vi_gam4)
dev.off()



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



