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
require(ggpubr)
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


write.csv(preg,"Pregnancy.csv",row.names = F)


# Population summary
aged <- cut(preg$AGE, breaks = c(20,30, 35, 40, 51))
agev <- data.frame(age = aged,LG=preg$LIGATION_GROUP) %>%
  mutate(LG = factor(LG,levels=c("Clip","Ring","Coagulation","Ligation/Resection"))) %>%

  na.omit()

pdf("Age.pdf",height = 5.5,width = 6.5)
ggplot(agev,aes(x=age,y = 100*(..count..)/sum(..count..),fill=LG)) +
  geom_bar() + my_style() +
  ylab("Percent in each group") + xlab("Age group (yrs)") +
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
  mutate(value = factor(value,levels=c("Identical","1-SPD","2-SPD","3-SPD"))) %>% 

  mutate(class = "Anastomosis") 



#pdf("Anastomosis.pdf",height = 5.5,width = 6.5)
gana <- ggplot(PyrAN, aes(x = value, y =  100*(..count..)/sum(..count..), fill = variable,alpha=variable)) +   # Fill column
  geom_bar(position = "dodge",fill="#c72e29")   +  my_style() +
  scale_alpha_manual(values=c(0.6,1)) +
  scale_fill_wsj(name = "") + 
  xlab("Location") +
  ylab("") + 
  theme(legend.spacing.x = unit(0.15, 'cm'),legend.position = "none") 
#dev.off()


# Diameter 
PyrDiam <-  preg[,c("L_DIAMETER_NUMERIC","R_DIAMETER_NUMERIC")]  %>%
  melt() %>%  mutate(variable = recode(variable, L_DIAMETER_NUMERIC = "Left",
                                       R_DIAMETER_NUMERIC = "Right")) %>%
              mutate(value = recode(value, "1" = "Similar",
                                           "2" = "Somewhat dissimilar",
                                           "3" = "Dissimilar")) %>% 
           mutate(class = "Diameter")  %>%
          mutate(value = factor(value,levels=c("Similar","Somewhat dissimilar","Dissimilar"))) 


#pdf("diameter.pdf",height = 5.5,width = 6.5)
gdiam <- ggplot(PyrDiam, aes(x = value, y =  100*(..count..)/sum(..count..),alpha=variable)) +   # Fill column
  geom_bar(position = "dodge",fill = "#016392")   +  my_style() +
  scale_alpha_manual(values=c(0.6,1)) +
  xlab("Diameter") +
  ylab("") + 
  theme(legend.spacing.x = unit(0.15, 'cm'),legend.position = "none") 
#dev.off()


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



#pdf("fibrosis.pdf",height = 5.5,width = 6.5)
gfib <- ggplot(PyrFib, aes(x = value, y =  100*(..count..)/sum(..count..), fill = variable,alpha=variable)) +   # Fill column
  geom_bar(position = "dodge",fill="#be9c2e")   +  my_style() +
  scale_alpha_manual(values=c(0.6,1)) +
  scale_fill_wsj(name = "") + 
  xlab("Fibrosis") +
  ylab("") + 
  theme(legend.spacing.x = unit(0.15, 'cm'),legend.position = "none") 
#dev.off()


# Tube length
PyrTL <-  preg[,c("LEFT_TUBE_LENGTH","RIGHT_TUBE_LENGTH")]  %>%
  melt() %>%  mutate(variable = recode(variable, LEFT_TUBE_LENGTH = "Left",
                                       RIGHT_TUBE_LENGTH = "Right"))  %>% 
  mutate(class = "Length (cm)") %>% 
  mutate(variable = factor(variable,levels=c("Left","Right"))) 






gTL <- ggplot(data= PyrTL,aes(x=value, alpha=variable,group=variable)) +
  geom_histogram(position='dodge',fill="#098154",binwidth = 0.5,aes(group=variable,y = 100*(..count..)/sum(..count..)))  +  my_style() +
  scale_fill_wsj(name = "") +
  scale_alpha_manual(values=c(0.6,1)) +
  xlab("Length (cm)") + theme(legend.position = "none") + ylab("") 


pdf("anatomy.pdf",height = 9.5,width = 11)
grid.arrange(gana,gdiam , gfib ,gTL,  ncol = 2,nrow=2,
             left =  text_grob("Percent in each group", size=18,rot=90) )

dev.off()



Py_all <- rbind(PyrAN,PyrDiam,PyrFib) %>% 
  mutate(value = factor(value,levels=c("Identical","1-SPD","2-SPD","3-SPD",
                                       "Similar","Somewhat dissimilar","Dissimilar","None","Mild","Moderate","Severe"))) 




pdf("anatomy.pdf",height = 11.5,width = 12)
gg1 <- ggplot(Py_all, aes(x = value, y =  100*(..count..)/sum(..count..), fill = class,alpha=variable)) +   # Fill column
  geom_bar(position = "dodge")   +  my_style() +
  scale_fill_wsj(name = "") +
  scale_alpha_manual(values=c(0.6,1)) +
  xlab("") +
  ylab(" Per cent in each group") + 
  theme(legend.position = "none") +
  facet_wrap(.~class,scale="free",nrow=2)

dev.off()

