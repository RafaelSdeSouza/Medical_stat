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
# Auxiliar function to randomly select a given column



outcomes <- read.csv("BTA-Pregnancies-anonymized.csv") %>% dplyr::select(c("OutcomeGpNumeric",
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
  dplyr::select(c("OutcomeGpNumeric",
                  "LIGATION_GROUP",
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

colnames(outcomes2) <- c("OutcomeGpNumeric","Sterilization_Method", "Age", "Length","Location","Fibrosis",
                     "Diameter")

#colnames(outcomes2) <- c("OutcomeGpNumeric", "Age", "Length","Location","Fibrosis",
#                         "Diameter")

outcomes3 <- outcomes2 %>%
  mutate(OutcomeGpNumeric = recode(OutcomeGpNumeric,
                                   "1"="Birth","2"="Ongoing","3" = "Miscarriage",
                                   "4" = "Ectopic"))

write.csv(outcomes3,"Outcomes.csv",row.names = F)
