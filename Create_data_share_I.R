# Pregnancy data
require(dplyr)
require(magrittr)
require(reshape)
require(forcats)

# Auxiliar function to randomly select a given column 


# Data-processing
preg <- read.csv("BTA-Patients-MAW.csv") %>% select(c("PREGNANT_NUMERIC",  "AGE", "LIGATION_GROUP", 
                                                      "LEFT_TUBE_LENGTH",     "RIGHT_TUBE_LENGTH",
                                                      "L_DIAMETER_NUMERIC",   "R_DIAMETER_NUMERIC",
                                                      "L_FIBROSIS_NUMERIC",   "R_FIBROSIS_NUMERIC",
                                                      "ANASTOMOSIS2_NUMERIC","ANASTOMOSIS1_NUMERIC"
)) %>%
  filter(AGE != "Yes") %>%  mutate(AGE = as.numeric(as.character(AGE)) )  %>%
  filter(AGE > 10) %>%
  na.omit() %>% mutate(LEFT_TUBE_LENGTH = as.numeric(as.character(LEFT_TUBE_LENGTH)) ) %>%
  filter(PREGNANT_NUMERIC %in% c(0,1)) %>% 
  mutate(PREGNANT_NUMERIC = as.numeric(as.character(PREGNANT_NUMERIC))) %>% 
  droplevels()  



# Create new dataset with choosen features
preg2 <- preg %>%
  mutate(PREGNANT_NUMERIC = as.factor(PREGNANT_NUMERIC)) %>%
  mutate(L_FIBROSIS_NUMERIC = recode(L_FIBROSIS_NUMERIC, "0" = "None",
                            "1" = "Mild",
                            "2" = "Moderate",
                            "3" = "Severe")) %>%
  
  mutate(R_FIBROSIS_NUMERIC = recode(R_FIBROSIS_NUMERIC, "0" = "None",
                                     "1" = "Mild",
                                     "2" = "Moderate",
                                     "3" = "Severe")) %>%
  
  mutate(ANASTOMOSIS2_NUMERIC = recode(ANASTOMOSIS2_NUMERIC, "0" = "Identical",
                            "1" = "1-SPD",
                            "2" = "2-SPD",
                            "3" = "3-SPD")) %>% 
  
  mutate(ANASTOMOSIS1_NUMERIC = recode(ANASTOMOSIS1_NUMERIC, "0" = "Identical",
                                         "1" = "1-SPD",
                                         "2" = "2-SPD",
                                         "3" = "3-SPD")) %>%  
  
   mutate(L_DIAMETER_NUMERIC = recode(L_DIAMETER_NUMERIC, "1" = "Similar","2" = "Somewhat dissimilar","3" = "Dissimilar")) %>% 
   mutate(R_DIAMETER_NUMERIC = recode(R_DIAMETER_NUMERIC, "1" = "Similar","2" = "Somewhat dissimilar","3" = "Dissimilar"))  %>% 
   mutate(PREGNANT_NUMERIC = recode(PREGNANT_NUMERIC, "1" = "Yes","0" = "No")) 


colnames(preg2) <- c("Pregnancy", "Age", "Ligation_Group", "Left_Length","Right_Length","Left_Diameter",
                     "Right_Diameter","Left_Fibrosis", "Right_Fibrosis",
                     "Left_Location" ,"Right_Location"     
                     )



write.csv(preg2,"Pregnancy.csv",row.names = F)
