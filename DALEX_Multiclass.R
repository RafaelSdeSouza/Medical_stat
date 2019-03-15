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



outcomes <- read.csv("BTA-Pregnancies-anonymized.csv") %>% select(c("OutcomeGpNumeric",  "AGE",    
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
           "AGE", "TL_rand","ANAS_rand","Fibr_rand",
           "Diam_rand")) 


# Split train vs test sample
trainIndex <- createDataPartition(outcomes2$OutcomeGpNumeric, p = .7, 
                                  list = FALSE, 
                                  times = 1)
Train <- outcomes2[trainIndex,] 
#%>% mutate(LIGATION_GROUP = as.numeric(LIGATION_GROUP))


Test  <- outcomes2[-trainIndex,]  
#%>% mutate(LIGATION_GROUP = as.numeric(LIGATION_GROUP))


classif_rf <-  randomForest(OutcomeGpNumeric~ 1+ 
                              AGE+
                              TL_rand + 
                              ANAS_rand + Fibr_rand + Diam_rand, 
                            data = Train, ntree=2000,nodesize=50)



classif_gam <- vgam(OutcomeGpNumeric~
       s(AGE,bs="cr",k=10)+
       s(TL_rand,bs="cr",k=10) + 
       ANAS_rand + Fibr_rand + Diam_rand,multinomial, data = Train)


loss_cross_entropy(Train[,1], yhat(classif_rf))

loss_cross_entropy(Train[,1],  predict(classif_gam,type="response"))



#1=birth
#2=ongoing
#3=miscarriage
#4=ectopic

pred1 <- function(m, x)   predict(m, x, type = "prob")[,1]
pred2 <- function(m, x)   predict(m, x, type = "prob")[,2]
pred3 <- function(m, x)   predict(m, x, type = "prob")[,3]
pred4 <- function(m, x)   predict(m, x, type = "prob")[,4]


explain_rf_1 <- explain(classif_rf, data = Test[,-1], 
                              y = Test$OutcomeGpNumeric == "1", 
                              predict_function = pred1, label = "Birth")
explain_rf_2 <- explain(classif_rf, data = Test[,-1], 
                           y = Test$OutcomeGpNumeric == "2",
                           predict_function = pred2, label = "Ongoing")
explain_rf_3 <- explain(classif_rf, data = Test[,-1], 
                                 y = Test$OutcomeGpNumeric == "3", 
                                 predict_function = pred3, label = "Miscarriage")

explain_rf_4 <- explain(classif_rf, data = Test[,-1], 
                          y = Test$OutcomeGpNumeric == "4",
                          predict_function = pred4, label = "Ectopic")







predgam1 <- function(m, x)   predict(m, x, type = "response")[,1]
predgam2 <- function(m, x)   predict(m, x, type = "response")[,2]
predgam3 <- function(m, x)   predict(m, x, type = "response")[,3]
predgam4 <- function(m, x)   predict(m, x, type = "response")[,4]


explain_gam_1 <- explain(classif_gam, data = Test[,-1], 
                        y = Test$OutcomeGpNumeric == "1", 
                        predict_function = predgam1, label = "Birth")
explain_gam_2 <- explain(classif_gam, data = Test[,-1], 
                        y = Test$OutcomeGpNumeric == "2",
                        predict_function = predgam2, label = "Ongoing")
explain_gam_3 <- explain(classif_gam, data = Test[,-1], 
                        y = Test$OutcomeGpNumeric == "3", 
                        predict_function = predgam3, label = "Miscarriage")

explain_gam_4 <- explain(classif_gam, data = Test[,-1], 
                        y = Test$OutcomeGpNumeric == "4",
                        predict_function = predgam4 , label = "Ectopic")



mp_rf_1 <- model_performance(explain_rf_1) %>% as.data.frame()
mp_rf_2 <- model_performance(explain_rf_2) %>% as.data.frame()
mp_rf_3 <- model_performance(explain_rf_3) %>% as.data.frame()
mp_rf_4 <- model_performance(explain_rf_4) %>% as.data.frame()

mp_rf_all <- rbind(mp_rf_1,mp_rf_2,mp_rf_3,mp_rf_4) %>%
  mutate(method="RF")

mp_gam_1 <- model_performance(explain_gam_1) %>% as.data.frame()
mp_gam_2 <- model_performance(explain_gam_2) %>% as.data.frame()
mp_gam_3 <- model_performance(explain_gam_3) %>% as.data.frame()
mp_gam_4 <- model_performance(explain_gam_4) %>% as.data.frame()

mp_gam_all <- rbind(mp_gam_1,mp_gam_2,mp_gam_3,mp_gam_4) %>%
  mutate(method="GAM")

mp_multi_all <- rbind(mp_rf_all,mp_gam_all) %>% mutate(diff = abs(diff))
  



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
ale_rf1_age   <- variable_response(explain_rf_1, variable =  "AGE", type = "ale")
ale_rf2_age   <- variable_response(explain_rf_2, variable =  "AGE", type = "ale")
ale_rf3_age   <- variable_response(explain_rf_3, variable =  "AGE", type = "ale")
ale_rf4_age   <- variable_response(explain_rf_4, variable =  "AGE", type = "ale")
ale_all_age <- rbind(ale_rf1_age,ale_rf2_age,ale_rf3_age,ale_rf4_age ) %>%
  mutate(method = "RF")



ale_gam1_age   <- variable_response(explain_gam_1, variable =  "AGE", type = "ale")
ale_gam2_age   <- variable_response(explain_gam_2, variable =  "AGE", type = "ale")
ale_gam3_age   <- variable_response(explain_gam_3, variable =  "AGE", type = "ale")
ale_gam4_age   <- variable_response(explain_gam_4, variable =  "AGE", type = "ale")
ale_gam_all_age <- rbind(ale_gam1_age,ale_gam2_age,ale_gam3_age,ale_gam4_age ) %>%
  mutate(method = "GAM")


ale_all <- rbind(ale_all_age,ale_gam_all_age )


pdf("ale_age_outcome.pdf",height = 5,width = 12)
g <- ggplot(ale_all,aes(x=x,y=y,group=label,color=label,fill=label)) + 
  geom_smooth(method = 'loess',span=0.25) +
  #+ geom_point() +
  my_style() +
  scale_linetype_stata(name="") +
  scale_color_wsj(name="") + ylab("Outcome likelihood") +
  scale_fill_wsj(name = "") +
  xlab("AGE (yr)") +
#+ 
#  coord_cartesian(xlim=c(21,50),ylim=c(0,0.8)) +
  facet_wrap(.~method)
direct.label(g,"angled.endpoints")
dev.off()


pdf("ale_age_outcome.pdf",height = 4,width = 6)
g <- ggplot(ale_all_age,aes(x=x,y=y,group=label,color=label,fill=label)) + 
  geom_smooth(method = 'loess',span=0.25) +
#+ geom_point() +
  my_style() +
  scale_linetype_stata(name="") +
  scale_color_wsj(name="") + ylab("Outcome likelihood") +
  scale_fill_wsj(name = "") +
  xlab("AGE (yr)") + 
  coord_cartesian(xlim=c(21,50),ylim=c(0,1))
direct.label(g,"angled.endpoints")
#+
#  theme(legend.position = "none" )
dev.off()



gg <- ggplot(ale_gam_all_age,aes(x=x,y=y,group=label,color=label,fill=label)) + 
  geom_smooth(method = 'loess',span=0.25) +
  #+ geom_point() +
  my_style() +
  scale_linetype_stata(name="") +
  scale_color_wsj(name="") + ylab("Outcome likelihood") +
  scale_fill_wsj(name = "") +
  xlab("AGE (yr)") + 
  coord_cartesian(xlim=c(21,50),ylim=c(0,1))
direct.label(gg,"angled.endpoints")




## Tube lenght properties

ale_rf1_TL1   <- variable_response(explain_rf_1, variable =  "TL_rand", type = "pdp")
ale_rf2_TL1   <- variable_response(explain_rf_2, variable =  "TL_rand", type = "pdp")
ale_rf3_TL1   <- variable_response(explain_rf_3, variable =  "TL_rand", type = "pdp")
ale_rf4_TL1   <- variable_response(explain_rf_4, variable =  "TL_rand", type = "pdp")


ale_rf1_TLD1   <- variable_response(explain_rf_1, variable =  "TLD_rand", type = "pdp")
ale_rf2_TLD1   <- variable_response(explain_rf_2, variable =  "TLD_rand", type = "pdp")
ale_rf3_TLD1   <- variable_response(explain_rf_3, variable =  "TLD_rand", type = "pdp")
ale_rf4_TLD1   <- variable_response(explain_rf_4, variable =  "TLD_rand", type = "pdp")


ale_rf1_TLP1   <- variable_response(explain_rf_1, variable =  "TLP_rand", type = "pdp")
ale_rf2_TLP1   <- variable_response(explain_rf_2, variable =  "TLP_rand", type = "pdp")
ale_rf3_TLP1   <- variable_response(explain_rf_3, variable =  "TLP_rand", type = "pdp")
ale_rf4_TLP1   <- variable_response(explain_rf_4, variable =  "TLP_rand", type = "pdp")

ale_TL_all <- rbind(ale_rf1_TL1,ale_rf2_TL1,ale_rf3_TL1,ale_rf4_TL1,
                    ale_rf1_TLD1,ale_rf2_TLD1,ale_rf3_TLD1,ale_rf4_TLD1,
                    ale_rf1_TLP1,ale_rf2_TLP1,ale_rf3_TLP1,ale_rf4_TLP1) %>%
  mutate(var = recode(var, TL_rand = "Tube length",
                      TLD_rand = "Tube length distal",
                      TLP_rand = "Tube length prox"))




pdf("ale_TL_multi.pdf",height = 4.5,width = 14)
g2 <- ggplot(ale_TL_all,aes(x=x,y=y,group=label,color=label,fill=label)) +
  geom_smooth(method = 'loess',span=0.275) + 
 # geom_point() +  
  my_style() +
  scale_linetype_stata(name="") +
  scale_color_wsj(name="") + ylab("Outcome likelihood") +
  scale_fill_wsj(name = "") + 
  ylab("Outcome likelihood") + xlab("Length (cm)") +
  coord_cartesian(xlim=c(-0.675,12)) + facet_wrap(.~var) 
direct.label(g2,"lasso.labels")
dev.off()




explain_rf <- explain(classif_rf, label = "RF",
                      data = Test[,-1], y = yTest,
                      predict_function = p_rf)






ale_rf_age   <- variable_response(explain_rf, variable =  "AGE", type = "pdp")



ggplot(ale_rf_age,aes( x = x, y = y, group=label,color=label,linetype=label)) +  
  geom_smooth(method = 'loess',n = 500,show.legend=F) + geom_point() +  my_style() +
  scale_color_fivethirtyeight() + ylab("Pregnancy likelihood") +
  scale_fill_fivethirtyeight() +
  xlab("AGE (yr)") + theme(legend.position = "none") +
  coord_cartesian(xlim=c(21,50))



