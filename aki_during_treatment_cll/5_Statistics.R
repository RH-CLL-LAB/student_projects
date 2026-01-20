# This script looks at statistics (cum. risk, cum. inc., cox uni and cox multi) for the CLL and SLL cohort
# Date: 25.06.2025
# Author: Emma Siig 
# Comment: Run 1cohort.R first, then 2Baseline_characteristics.R, then 3AdverseEvents.R. Or downloade table: cohort_riskfactors_AE1

source('/ngc/projects2/dalyca_r/clean_r/load_dalycare_package.R')
library(cmprsk)
#cohort_riskfactors_AE1 <- readRDS("/ngc/projects2/dalyca_r/emmsii_r/sandbox/tables/3cohort_riskfactors_AE1.rds")
cohort_riskfactors_AE2 <- readRDS("/ngc/projects2/dalyca_r/emmsii_r/sandbox/tables/4cohort_riskfactors_AE2.rds")

cohort_done_statistics <- cohort_riskfactors_AE2

#### Making age groups ####
cohort_done_statistics <- cohort_done_statistics %>%        #Making age groups using age at treatment
  mutate(age_group = case_when(
    age_at_treatment < 60 ~ "<60",
    age_at_treatment >= 60 & age_at_treatment < 70 ~ "60-69",
    age_at_treatment >= 70 & age_at_treatment < 80 ~ "70-79",
    age_at_treatment >= 80 ~ "80+"
  ))




#### Kaplan-Meier estimated cumulative risk of adverse events ####
#Anemia
anemia_data_cum_inc <- cohort_done_statistics %>%            #Filters rows with NA
  filter(!is.na(date_treatment_1st_line),
         !is.na(anemia_event_or_censor_date),
         !is.na(anemia_binary))

anemia_data_cum_inc <- anemia_data_cum_inc %>%                #Makes variable for time to anemia or event
  mutate(
    anemia_time = as.numeric(anemia_event_or_censor_date - date_treatment_1st_line)
  )

print(anemia_data_cum_inc)

surv_anemia_obj <- Surv(                          #Makes a survival object
  time = anemia_data_cum_inc$anemia_time,
  event = anemia_data_cum_inc$anemia_binary
)

fit_anemia <- survfit(surv_anemia_obj ~ 1)        #Making a survival object for the whole data set without stratification 

ggsurvplot(                                       #Makes curv for cumulative risk for anemia
  fit_anemia,
  data = anemia_data_cum_inc,
  fun = "event",                                  #Makes it cumulative
  conf.int = TRUE,                                #Includes confidence interval 
  risk.table = TRUE,                              #Includes table under plot with number of patients at risk
  xlab = "Days since treatment start",
  ylab = "Estimated cumulative risk of anemia (Kaplan-Meier)",
  title = "Kaplan-Meier estimated cumulative risk of anemia within 180 days after first treatment"
)

#Neutropenia
neutropenia_data_cum_inc <- cohort_done_statistics %>%            #Filters rows with NA
  filter(!is.na(date_treatment_1st_line),
         !is.na(neutropenia_event_or_censor_date),
         !is.na(neutropenia_binary))

neutropenia_data_cum_inc <- neutropenia_data_cum_inc %>%                #Makes variable for time to neutropenia or event
  mutate(
    neutropenia_time = as.numeric(neutropenia_event_or_censor_date - date_treatment_1st_line)
  )

surv_neutropenia_obj <- Surv(                          #Makes surv object
  time = neutropenia_data_cum_inc$neutropenia_time,
  event = neutropenia_data_cum_inc$neutropenia_binary
)

fit_neutropenia <- survfit(surv_neutropenia_obj ~ 1)        #Fits model

ggsurvplot(                                       #Makes cumulative risk curve for neutropenia
  fit_neutropenia,
  data = neutropenia_data_cum_inc,
  fun = "event",                                  #Makes it cumulative
  conf.int = TRUE,                                #Includes confidence interval 
  risk.table = TRUE,                              #Includes table under plot with number of patients at risk
  xlab = "Days since treatment start",
  ylab = "Estimated cumulative risk of neutropenia (Kaplan-Meier)",
  title = "Kaplan-Meier estimated cumulative risk of neutropenia within 180 days after first treatment"
)

#Lymphopenia
lymphopenia_data_cum_inc <- cohort_done_statistics %>%            #Filters rows with NA
  filter(!is.na(date_treatment_1st_line),
         !is.na(lymphopenia_event_or_censor_date),
         !is.na(lymphopenia_binary))

lymphopenia_data_cum_inc <- lymphopenia_data_cum_inc %>%                #Makes variable for time to lymphopenia or event
  mutate(
    lymphopenia_time = as.numeric(lymphopenia_event_or_censor_date - date_treatment_1st_line)
  )

surv_lymphopenia_obj <- Surv(                          #Makes surv object
  time = lymphopenia_data_cum_inc$lymphopenia_time,
  event = lymphopenia_data_cum_inc$lymphopenia_binary
)

fit_lymphopenia <- survfit(surv_lymphopenia_obj ~ 1)        #Fits model

ggsurvplot(                                       #Makes cumulative risk curve for lymphopenia
  fit_lymphopenia,
  data = lymphopenia_data_cum_inc,
  fun = "event",                                  #Makes it cumulative
  conf.int = TRUE,                                #Includes confidence interval 
  risk.table = TRUE,                              #Includes table under plot with number of patients at risk
  xlab = "Days since treatment start",
  ylab = "Estimated cumulative risk of lymphopenia (Kaplan-Meier)",
  title = "Kaplan-Meier estimated cumulative risk of lymphopenia within 180 days after first treatment"
)

#Thrombocytopenia (platelets)
thrombocytopenia_data_cum_inc <- cohort_done_statistics %>%            #Filters rows with NA
  filter(!is.na(date_treatment_1st_line),
         !is.na(thrombocytopenia_event_or_censor_date),
         !is.na(thrombocytopenia_binary))

thrombocytopenia_data_cum_inc <- thrombocytopenia_data_cum_inc %>%                #Makes variable for time to thrombocytopenia or event
  mutate(
    thrombocytopenia_time = as.numeric(thrombocytopenia_event_or_censor_date - date_treatment_1st_line)
  )

surv_thrombocytopenia_obj <- Surv(                          #Makes surv object
  time = thrombocytopenia_data_cum_inc$thrombocytopenia_time,
  event = thrombocytopenia_data_cum_inc$thrombocytopenia_binary
)

fit_thrombocytopenia <- survfit(surv_thrombocytopenia_obj ~ 1)        #Fits model

ggsurvplot(                                       #Makes cumulative risk curve for thrombocytopenia
  fit_thrombocytopenia,
  data = thrombocytopenia_data_cum_inc,
  fun = "event",                                  #Makes it cumulative
  conf.int = TRUE,                                #Includes confidence interval 
  risk.table = TRUE,                              #Includes table under plot with number of patients at risk
  xlab = "Days since treatment start",
  ylab = "Estimated cumulative risk of thrombocytopenia (Kaplan-Meier)",
  title = "Kaplan-Meier estimated cumulative risk of thrombocytopenia within 180 days after first treatment"
)

#Leukopenia (WBC)
leukopenia_data_cum_inc <- cohort_done_statistics %>%            #Filters rows with NA
  filter(!is.na(date_treatment_1st_line),
         !is.na(leukopenia_event_or_censor_date),
         !is.na(leukopenia_binary))

leukopenia_data_cum_inc <- leukopenia_data_cum_inc %>%                #Makes variable for time to leukopenia or event
  mutate(
    leukopenia_time = as.numeric(leukopenia_event_or_censor_date - date_treatment_1st_line)
  )

surv_leukopenia_obj <- Surv(                          #Makes surv object
  time = leukopenia_data_cum_inc$leukopenia_time,
  event = leukopenia_data_cum_inc$leukopenia_binary
)

fit_leukopenia <- survfit(surv_leukopenia_obj ~ 1)        #Fits model

ggsurvplot(                                       #Makes cumulative risk curve for leukopenia
  fit_leukopenia,
  data = leukopenia_data_cum_inc,
  fun = "event",                                  #Makes it cumulative
  conf.int = TRUE,                                #Includes confidence interval 
  risk.table = TRUE,                              #Includes table under plot with number of patients at risk
  xlab = "Days since treatment start",
  ylab = "Estimated cumulative risk of leukopenia (Kaplan-Meier)",
  title = "Kaplan-Meier estimated cumulative risk of leukopenia within 180 days after first treatment"
)


#Creatinine
creatinine_data_cum_inc <- cohort_done_statistics %>%            #Filters rows with NA
  filter(!is.na(date_treatment_1st_line),
         !is.na(creatinine_event_or_censor_date),
         !is.na(creatinine_binary))

creatinine_data_cum_inc <- creatinine_data_cum_inc %>%                #Makes variable for time to creatinine or event
  mutate(
    creatinine_time = as.numeric(creatinine_event_or_censor_date - date_treatment_1st_line)
  )

surv_creatinine_obj <- Surv(                          #Makes surv object
  time = creatinine_data_cum_inc$creatinine_time,
  event = creatinine_data_cum_inc$creatinine_binary
)

fit_creatinine <- survfit(surv_creatinine_obj ~ 1)        #Fits model

ggsurvplot(                                       #Makes cumulative risk curve for creatinine increase (AKI)
  fit_creatinine,
  data = creatinine_data_cum_inc,
  fun = "event",                                  #Makes it cumulative
  conf.int = TRUE,                                #Includes confidence interval 
  risk.table = TRUE,                              #Includes table under plot with number of patients at risk
  xlab = "Days since treatment start",
  ylab = "Estimated cumulative risk of creatinine increase (AKI) (Kaplan-Meier)",
  title = "Kaplan-Meier estimated cumulative risk of creatinine increase (AKI) within 180 days after first treatment"
)


#Plot with all cumulative risk curves in one

anemia_data_cum_inc$ae_type <- "Anemia"
neutropenia_data_cum_inc$ae_type <- "Neutropenia"
lymphopenia_data_cum_inc$ae_type <- "Lymphopenia"
thrombocytopenia_data_cum_inc$ae_type <- "Thrombocytopenia"
leukopenia_data_cum_inc$ae_type <- "Leukopenia"
creatinine_data_cum_inc$ae_type <- "Creatinine"

anemia_data_prep <- anemia_data_cum_inc %>%                          #Anemia prep
  dplyr::rename(time = anemia_time, event = anemia_binary)

neutropenia_data_prep <- neutropenia_data_cum_inc %>%                #Neutropenia prep
  dplyr::rename(time = neutropenia_time, event = neutropenia_binary)

lymphopenia_data_prep <- lymphopenia_data_cum_inc %>%                #lymphopenia prep
  dplyr::rename(time = lymphopenia_time, event = lymphopenia_binary)

thrombocytopenia_data_prep <- thrombocytopenia_data_cum_inc %>%      #thrombocytopenia prep
  dplyr::rename(time = thrombocytopenia_time, event = thrombocytopenia_binary)

leukopenia_data_prep <- leukopenia_data_cum_inc %>%                  #leukopenia prep
  dplyr::rename(time = leukopenia_time, event = leukopenia_binary)

creatinine_data_prep <- creatinine_data_cum_inc %>%                  #creatinine prep
  dplyr::rename(time = creatinine_time, event = creatinine_binary)


ae_data <- bind_rows(anemia_data_prep, neutropenia_data_prep, lymphopenia_data_prep, thrombocytopenia_data_prep, leukopenia_data_prep, creatinine_data_prep)

ae_data_surv_obj <- Surv(time = ae_data$time, event = ae_data$event) #Making survival-object stratified for ae_type

fit_ae <- survfit(ae_data_surv_obj ~ ae_type, data = ae_data) #Kaplan-Meier fit with stratification for ae_type

ggsurvplot(fit_ae,                         #Plot with all curves
           data = ae_data,
           fun = "event",                  #Cumulative risk
           conf.int = TRUE, 
           risk.table = TRUE, 
           xlab = "Days since treatment start", 
           ylab = "Estimated cumulative risk of adverse event", 
           title = "Kaplan-Meier estimated cumulative risk of adverse events")


#### Cumulative incidence curves (Aalen-Johansen) ####
#_________Anemia _____________
nrow(dead_without_anemia_in_followup)              #73 patients died within 180 days and without Anemia


anemia_aj_data <- anemia_data_cum_inc %>%          #Makes new variable: competing risk status (cr_status)
  mutate(anemia_cr_status = case_when(
    anemia_binary == 1 ~ 1,                                                                  # 1 = Event (anemia)
    anemia_binary == 0 & status == 1 & date_death_fu <= date_treatment_1st_line + 180 ~ 2,   # 2 = Compering risk (death)
    TRUE ~ 0                                                                                 # 0 = censored (no anemia and no death)
  ))

anemia_aj_data %>% count(anemia_cr_status) #Checking the numbers add up
table(cohort_riskfactors_AE1$anemia_binary)

anemia_aj_ci <- cmprsk::cuminc(                 #Cumulative incidence for anemia, with death as a competing risk        
  ftime = anemia_aj_data$anemia_time,           #Time to event or censor
  fstatus = anemia_aj_data$anemia_cr_status,    #Status: 1 = anemia, 2 = death, 0 = censored
  cencode = 0
)

names(anemia_aj_ci)                             #Renames
names(anemia_aj_ci)[names(anemia_aj_ci) == "1 1"] <- "Anemia"
names(anemia_aj_ci)[names(anemia_aj_ci) == "1 2"] <- "Death"


plot(anemia_aj_ci,                             #Plotting curves
     xlab = "Days since treatment start",
     ylab = "Estimated cumulative incidence of anemia",
     col = c("red", "grey"),
     lwd = 2,
     ylim = c(0, 0.25))                        #Making the y-axis go to 0.25


#____________Neutropenia _____________
nrow(dead_without_neutropenia_in_followup)          #75 patients died within 180 days and without neutropenia


neutropenia_aj_data <- neutropenia_data_cum_inc %>%          #Makes new variable: competing risk status (cr_status)
  mutate(neutropenia_cr_status = case_when(
    neutropenia_binary == 1 ~ 1,                                                                  # 1 = Event (neutropenia)
    neutropenia_binary == 0 & status == 1 & date_death_fu <= date_treatment_1st_line + 180 ~ 2,   # 2 = Compering risk (death)
    TRUE ~ 0                                                                                      # 0 = censored (no neutropenia and no death)
  ))

neutropenia_aj_data %>% count(neutropenia_cr_status) #Checking the numbers add up
table(cohort_riskfactors_AE1$neutropenia_binary)

neutropenia_aj_ci <- cmprsk::cuminc(                      #Cumulative incidence for neutropenia, with death as a competing risk        
  ftime = neutropenia_aj_data$neutropenia_time,           #Time to event or censor
  fstatus = neutropenia_aj_data$neutropenia_cr_status,    #Status: 1 = neutropenia, 2 = death, 0 = censored
  cencode = 0
)

names(neutropenia_aj_ci)                             #Renames
names(neutropenia_aj_ci)[names(neutropenia_aj_ci) == "1 1"] <- "Neutropenia"
names(neutropenia_aj_ci)[names(neutropenia_aj_ci) == "1 2"] <- "Death"


plot(neutropenia_aj_ci,                             #Plotting curves
     xlab = "Days since treatment start",
     ylab = "Estimated cumulative incidence of neutropenia",
     col = c("red", "grey"),
     lwd = 2,
     ylim = c(0, 0.5))                        #Making the y-axis go to 0.5


#_________Lymphopenia _____________
nrow(dead_without_lymphopenia_in_followup)          #76 patients died within 180 days and without lymphopenia


lymphopenia_aj_data <- lymphopenia_data_cum_inc %>%          #Makes new variable: competing risk status (cr_status)
  mutate(lymphopenia_cr_status = case_when(
    lymphopenia_binary == 1 ~ 1,                                                                  # 1 = Event (lymphopenia)
    lymphopenia_binary == 0 & status == 1 & date_death_fu <= date_treatment_1st_line + 180 ~ 2,   # 2 = Compering risk (death)
    TRUE ~ 0                                                                                      # 0 = censored (no anemia and no death)
  ))

lymphopenia_aj_data %>% count(lymphopenia_cr_status) #Checking the numbers add up
table(cohort_riskfactors_AE1$lymphopenia_binary)

lymphopenia_aj_ci <- cmprsk::cuminc(                      #Cumulative incidence for lymphopenia, with death as a competing risk        
  ftime = lymphopenia_aj_data$lymphopenia_time,           #Time to event or censor
  fstatus = lymphopenia_aj_data$lymphopenia_cr_status,    #Status: 1 = lymphopenia, 2 = death, 0 = censored
  cencode = 0
)

names(lymphopenia_aj_ci)                 #Renames
names(lymphopenia_aj_ci)[names(lymphopenia_aj_ci) == "1 1"] <- "Lymphopenia"
names(lymphopenia_aj_ci)[names(lymphopenia_aj_ci) == "1 2"] <- "Death"


plot(lymphopenia_aj_ci,                             #Plotting curves
     xlab = "Days since treatment start",
     ylab = "Estimated cumulative incidence of lymphopenia",
     col = c("red", "grey"),
     lwd = 2,
     ylim = c(0, 0.5))                        #Making the y-axis go to 0.5


#__________Thrombocytopenia_____________
nrow(dead_without_thrombocytopenia_in_followup)          #89 patients died within 180 days and without thrombocytopenia


thrombocytopenia_aj_data <- thrombocytopenia_data_cum_inc %>%          #Makes new variable: competing risk status (cr_status)
  mutate(thrombocytopenia_cr_status = case_when(
    thrombocytopenia_binary == 1 ~ 1,                                                                  # 1 = Event (thrombocytopenia)
    thrombocytopenia_binary == 0 & status == 1 & date_death_fu <= date_treatment_1st_line + 180 ~ 2,   # 2 = Compering risk (death)
    TRUE ~ 0                                                                                           # 0 = censored (no anemia and no death)
  ))

thrombocytopenia_aj_data %>% count(thrombocytopenia_cr_status) #Checking the numbers add up
table(cohort_riskfactors_AE1$thrombocytopenia_binary)

thrombocytopenia_aj_ci <- cmprsk::cuminc(                           #Cumulative incidence for thrombocytopenia, with death as a competing risk        
  ftime = thrombocytopenia_aj_data$thrombocytopenia_time,           #Time to event or censor
  fstatus = thrombocytopenia_aj_data$thrombocytopenia_cr_status,    #Status: 1 = thrombocytopenia, 2 = death, 0 = censored
  cencode = 0
)

names(thrombocytopenia_aj_ci)                 #Renames
names(thrombocytopenia_aj_ci)[names(thrombocytopenia_aj_ci) == "1 1"] <- "Thrombocytopenia"
names(thrombocytopenia_aj_ci)[names(thrombocytopenia_aj_ci) == "1 2"] <- "Death"


plot(thrombocytopenia_aj_ci,                             #Plotting curves
     xlab = "Days since treatment start",
     ylab = "Estimated cumulative incidence of thrombocytopenia",
     col = c("red", "grey"),
     lwd = 2,
     ylim = c(0, 0.2))                        #Making the y-axis go to 0.2


#__________Leukopenia_____________
nrow(dead_without_leukopenia_in_followup)          #99 patients died within 180 days and without leukopenia


leukopenia_aj_data <- leukopenia_data_cum_inc %>%          #Makes new variable: competing risk status (cr_status)
  mutate(leukopenia_cr_status = case_when(
    leukopenia_binary == 1 ~ 1,                                                                  # 1 = Event (leukopenia)
    leukopenia_binary == 0 & status == 1 & date_death_fu <= date_treatment_1st_line + 180 ~ 2,   # 2 = Compering risk (death)
    TRUE ~ 0                                                                                           # 0 = censored (no anemia and no death)
  ))

leukopenia_aj_data %>% count(leukopenia_cr_status) #Checking the numbers add up
table(cohort_riskfactors_AE1$leukopenia_binary)

leukopenia_aj_ci <- cmprsk::cuminc(                     #Cumulative incidence for leukopenia, with death as a competing risk        
  ftime = leukopenia_aj_data$leukopenia_time,           #Time to event or censor
  fstatus = leukopenia_aj_data$leukopenia_cr_status,    #Status: 1 = leukopenia, 2 = death, 0 = censored
  cencode = 0
)

names(leukopenia_aj_ci)                 #Renames
names(leukopenia_aj_ci)[names(leukopenia_aj_ci) == "1 1"] <- "Leukopenia"
names(leukopenia_aj_ci)[names(leukopenia_aj_ci) == "1 2"] <- "Death"


plot(leukopenia_aj_ci,                             #Plotting curves
     xlab = "Days since treatment start",
     ylab = "Estimated cumulative incidence of leukopenia",
     col = c("red", "grey"),
     lwd = 2,
     ylim = c(0, 0.35))                        #Making the y-axis go to 0.35


#__________Creatinine increase (AKI)_____________
nrow(dead_without_creatinine_in_followup)          #82 patients died within 180 days and without creatinine increase (AKI) 


creatinine_aj_data <- creatinine_data_cum_inc %>%          #Makes new variable: competing risk status (cr_status)
  mutate(creatinine_cr_status = case_when(
    creatinine_binary == 1 ~ 1,                                                                  # 1 = Event (creatinine increase (AKI))
    creatinine_binary == 0 & status == 1 & date_death_fu <= date_treatment_1st_line + 180 ~ 2,   # 2 = Compering risk (death)
    TRUE ~ 0                                                                                           # 0 = censored (no anemia and no death)
  ))

creatinine_aj_data %>% count(creatinine_cr_status) #Checking the numbers add up
table(cohort_riskfactors_AE1$creatinine_binary)

creatinine_aj_ci <- cmprsk::cuminc(                     #Cumulative incidence for creatinine, with death as a competing risk        
  ftime = creatinine_aj_data$creatinine_time,           #Time to event or censor
  fstatus = creatinine_aj_data$creatinine_cr_status,    #Status: 1 = creatinine increase (AKI), 2 = death, 0 = censored
  cencode = 0
)

names(creatinine_aj_ci)                 #Renames
names(creatinine_aj_ci)[names(creatinine_aj_ci) == "1 1"] <- "Creatinine increase (AKI)"
names(creatinine_aj_ci)[names(creatinine_aj_ci) == "1 2"] <- "Death"


plot(creatinine_aj_ci,                             #Plotting curves
     xlab = "Days since treatment start",
     ylab = "Estimated cumulative incidence of creatinine increase (AKI)",
     col = c("red", "grey"),
     lwd = 2,
     ylim = c(0, 0.2))                        #Making the y-axis go to 0.2


#Grouping all cumulative incidence plots in one

#Making function
cuminc_to_df <- function(cuminc_obj, varname, element_name) {     #Function to converting the cuminc-object to a data frame
  est <- cuminc_obj[[element_name]]                         #Status 1 for the variable
  tibble(
    time = est$time, 
    est = est$est, 
    var = varname
  )
}
names(creatinine_aj_ci) #Checking names to know which element to extract 

#Prepping variables
df_anemia <- cuminc_to_df(anemia_aj_ci, "Anemia", "Anemia")                 #Converting anemia cuminc-object to data frame
df_neutropenia <- cuminc_to_df(neutropenia_aj_ci, "Neutropenia", "Neutropenia")  #Converting neutropenia cuminc-object to data frame
df_lymphopenia <- cuminc_to_df(lymphopenia_aj_ci, "Lymphopenia", "Lymphopenia")  #Converting lymphopenia cuminc-object to data frame
df_thrombocytopenia <- cuminc_to_df(thrombocytopenia_aj_ci, "Thrombocytopenia", "Thrombocytopenia")  #Converting thrombocytopenia cuminc-object to data frame
df_leukopenia <- cuminc_to_df(leukopenia_aj_ci, "Leukopenia", "Leukopenia")  #Converting leukopenia cuminc-object to data frame
df_creatinine <- cuminc_to_df(creatinine_aj_ci, "Acute Kidney Injury", "Creatinine increase (AKI)")  #Converting creatinine cuminc-object to data frame

#Combining all variables
all_ci_data <- bind_rows(df_anemia, df_neutropenia, df_lymphopenia, df_thrombocytopenia, df_leukopenia, df_creatinine)   #Combining all the data frames of the cumulative incidences

#Making the ggplot
ggplot(all_ci_data, aes(x=time, y=est, color=var)) +
  geom_step(linewidth = 2) + 
  labs(
    title = "Cumulative incidence of adverse events",
    x = "Days since treatment start", 
    y = "Estimated cumulative incidence of adverse events",
    color = "Adverse event"
  ) +
  theme_minimal(base_size=18) +
  theme(
    plot.title = element_text(size=24, face="bold", hjust=0,5),
    axis.title = element_text(size=20, face="bold"),
    axis.text = element_text(size=18),
    legend.title = element_text(size=18, face="bold"),
    legend.text = element_text(size=16),
    legend.position = "bottom")


#### Cox proportional hazard models (univariable) ####
#Prepping for IGHV, TP53_ab and del17p_at_treatment
#IGHV
cohort_done_statistics$IGHV_grouped <- as.character(cohort_done_statistics$IGHV)
cohort_done_statistics$IGHV_grouped[is.na(cohort_done_statistics$IGHV_grouped)] <- "Unknown"
cohort_done_statistics$IGHV_grouped <- factor(cohort_done_statistics$IGHV_grouped, levels = c("Mutated", "Unmutated", "Unknown"))

#TP53_ab
cohort_done_statistics$TP53_ab_grouped <- as.character(cohort_done_statistics$TP53_ab)
cohort_done_statistics$TP53_ab_grouped[is.na(cohort_done_statistics$TP53_ab_grouped)] <- "Unknown"
cohort_done_statistics$TP53_ab_grouped <- factor(cohort_done_statistics$TP53_ab_grouped, levels = c("Yes", "No", "Unknown"))

#del17p_at_treatment
cohort_done_statistics$del17p_at_treatment_grouped <- as.character(cohort_done_statistics$del17p_at_treatment)
cohort_done_statistics$del17p_at_treatment_grouped[is.na(cohort_done_statistics$del17p_at_treatment_grouped)] <- "Unknown"
cohort_done_statistics$del17p_at_treatment_grouped <- factor(cohort_done_statistics$del17p_at_treatment_grouped, levels = c("Yes", "No", "Unknown"))


#Anemia
cox_anemia_data <- cohort_done_statistics %>%
  filter(!is.na(anemia_event_or_censor_date),
         !is.na(date_treatment_1st_line),
         !is.na(anemia_binary),
         !is.na(sex))

cox_uni_anemia_sex <- coxph(             #Sex
  Surv(as.numeric(anemia_event_or_censor_date - date_treatment_1st_line), anemia_binary) ~ sex,
  data = cox_anemia_data
)
summary(cox_uni_anemia_sex)

cox_anemia_data$age_group <- factor(cox_anemia_data$age_group)

cox_uni_anemia_age <- coxph(            #Age
  Surv(as.numeric(anemia_event_or_censor_date - date_treatment_1st_line), anemia_binary) ~ age_group, data = cox_anemia_data)

summary(cox_uni_anemia_age)


cox_anemia_data$treatment_type <- factor(cox_anemia_data$treatment_type)  #Treatment

table(cox_anemia_data$treatment_type, cox_anemia_data$anemia_binary)
cox_anemia_data <- cox_anemia_data %>%     #Filters ibrutinib-venetoclax away since it only has one AE
  filter(treatment_type != "ibrutinib-venetoclax") 

cox_anemia_data$treatment_type <- relevel(factor(cox_anemia_data$treatment_type), ref = "btk-inhibitor")

cox_uni_anemia_treatment <- coxph(      
  Surv(as.numeric(anemia_event_or_censor_date - date_treatment_1st_line), anemia_binary) ~ treatment_type, data = cox_anemia_data)

summary(cox_uni_anemia_treatment)

#Univariable analysis for IGHV, TP53_ab and del17p_at_treatment
#IGHV
cox_anemia_data$IGHV_grouped <- relevel(factor(cox_anemia_data$IGHV_grouped), ref = "Unmutated")
cox_uni_anemia_IGHV <- coxph(            
  Surv(as.numeric(anemia_event_or_censor_date - date_treatment_1st_line), anemia_binary) ~ IGHV_grouped, data = cox_anemia_data)
summary(cox_uni_anemia_IGHV)

#TP53_ab
cox_anemia_data$TP53_ab_grouped <- relevel(factor(cox_anemia_data$TP53_ab_grouped), ref = "No")
cox_uni_anemia_TP53 <- coxph(            
  Surv(as.numeric(anemia_event_or_censor_date - date_treatment_1st_line), anemia_binary) ~ TP53_ab_grouped, data = cox_anemia_data)
summary(cox_uni_anemia_TP53)

#del17p_at_treatment
cox_anemia_data$del17p_at_treatment_grouped <- relevel(factor(cox_anemia_data$del17p_at_treatment_grouped), ref = "No")
cox_uni_anemia_del17p <- coxph(            
  Surv(as.numeric(anemia_event_or_censor_date - date_treatment_1st_line), anemia_binary) ~ del17p_at_treatment_grouped, data = cox_anemia_data)
summary(cox_uni_anemia_del17p)


#Neutropenia
cox_neutropenia_data <- cohort_done_statistics %>%
  filter(!is.na(neutropenia_event_or_censor_date),
         !is.na(date_treatment_1st_line),
         !is.na(neutropenia_binary),
         !is.na(treatment_type),
         !is.na(sex))

cox_uni_neutropenia_sex <- coxph(            #Sex
  Surv(as.numeric(neutropenia_event_or_censor_date - date_treatment_1st_line), neutropenia_binary) ~ sex,
  data = cox_neutropenia_data
)
summary(cox_uni_neutropenia_sex)

cox_neutropenia_data$age_group <- factor(cox_neutropenia_data$age_group)

cox_uni_neutropenia_age <- coxph(          #Age
  Surv(as.numeric(neutropenia_event_or_censor_date - date_treatment_1st_line), neutropenia_binary) ~ age_group, data = cox_neutropenia_data)

summary(cox_uni_neutropenia_age)

cox_neutropenia_data$treatment_type <- factor(cox_neutropenia_data$treatment_type)  #Treatment

cox_neutropenia_data$treatment_type <- relevel(factor(cox_neutropenia_data$treatment_type), ref = "btk-inhibitor")

cox_uni_neutropenia_treatment <- coxph(      
  Surv(as.numeric(neutropenia_event_or_censor_date - date_treatment_1st_line), neutropenia_binary) ~ treatment_type, data = cox_neutropenia_data)

summary(cox_uni_neutropenia_treatment)


#Univariable analysis for IGHV, TP53_ab and del17p_at_treatment
#IGHV
cox_neutropenia_data$IGHV_grouped <- relevel(factor(cox_neutropenia_data$IGHV_grouped), ref = "Unmutated")
cox_uni_neutropenia_IGHV <- coxph(            
  Surv(as.numeric(neutropenia_event_or_censor_date - date_treatment_1st_line), neutropenia_binary) ~ IGHV_grouped, data = cox_neutropenia_data)
summary(cox_uni_neutropenia_IGHV)

#TP53_ab
cox_neutropenia_data$TP53_ab_grouped <- relevel(factor(cox_neutropenia_data$TP53_ab_grouped), ref = "No")
cox_uni_neutropenia_TP53 <- coxph(            
  Surv(as.numeric(neutropenia_event_or_censor_date - date_treatment_1st_line), neutropenia_binary) ~ TP53_ab_grouped, data = cox_neutropenia_data)
summary(cox_uni_neutropenia_TP53)

#del17p_at_treatment
cox_neutropenia_data$del17p_at_treatment_grouped <- relevel(factor(cox_neutropenia_data$del17p_at_treatment_grouped), ref = "No")
cox_uni_neutropenia_del17p <- coxph(            
  Surv(as.numeric(neutropenia_event_or_censor_date - date_treatment_1st_line), neutropenia_binary) ~ del17p_at_treatment_grouped, data = cox_neutropenia_data)
summary(cox_uni_neutropenia_del17p)


#Lymphopenia
cox_lymphopenia_data <- cohort_done_statistics %>%
  filter(!is.na(lymphopenia_event_or_censor_date),
         !is.na(date_treatment_1st_line),
         !is.na(lymphopenia_binary),
         !is.na(treatment_type),
         !is.na(sex))

cox_uni_lymphopenia_sex <- coxph(            #Sex
  Surv(as.numeric(lymphopenia_event_or_censor_date - date_treatment_1st_line), lymphopenia_binary) ~ sex,
  data = cox_lymphopenia_data
)
summary(cox_uni_lymphopenia_sex)

cox_lymphopenia_data$age_group <- factor(cox_lymphopenia_data$age_group)

cox_uni_lymphopenia_age <- coxph(          #Age
  Surv(as.numeric(lymphopenia_event_or_censor_date - date_treatment_1st_line), lymphopenia_binary) ~ age_group, data = cox_lymphopenia_data)

summary(cox_uni_lymphopenia_age)


cox_lymphopenia_data$treatment_type <- factor(cox_lymphopenia_data$treatment_type)  #Treatment

cox_lymphopenia_data$treatment_type <- relevel(factor(cox_lymphopenia_data$treatment_type), ref = "btk-inhibitor")

cox_uni_lymphopenia_treatment <- coxph(      
  Surv(as.numeric(lymphopenia_event_or_censor_date - date_treatment_1st_line), lymphopenia_binary) ~ treatment_type, data = cox_lymphopenia_data)

summary(cox_uni_lymphopenia_treatment)

#Univariable analysis for IGHV, TP53_ab and del17p_at_treatment
#IGHV
cox_lymphopenia_data$IGHV_grouped <- relevel(factor(cox_lymphopenia_data$IGHV_grouped), ref = "Unmutated")
cox_uni_lymphopenia_IGHV <- coxph(            
  Surv(as.numeric(lymphopenia_event_or_censor_date - date_treatment_1st_line), lymphopenia_binary) ~ IGHV_grouped, data = cox_lymphopenia_data)
summary(cox_uni_lymphopenia_IGHV)

#TP53_ab
cox_lymphopenia_data$TP53_ab_grouped <- relevel(factor(cox_lymphopenia_data$TP53_ab_grouped), ref = "No")
cox_uni_lymphopenia_TP53 <- coxph(            
  Surv(as.numeric(lymphopenia_event_or_censor_date - date_treatment_1st_line), lymphopenia_binary) ~ TP53_ab_grouped, data = cox_lymphopenia_data)
summary(cox_uni_lymphopenia_TP53)

#del17p_at_treatment
cox_lymphopenia_data$del17p_at_treatment_grouped <- relevel(factor(cox_lymphopenia_data$del17p_at_treatment_grouped), ref = "No")
cox_uni_lymphopenia_del17p <- coxph(            
  Surv(as.numeric(lymphopenia_event_or_censor_date - date_treatment_1st_line), lymphopenia_binary) ~ del17p_at_treatment_grouped, data = cox_lymphopenia_data)
summary(cox_uni_lymphopenia_del17p)


#Thrombocytopenia
cox_thrombocytopenia_data <- cohort_done_statistics %>%
  filter(!is.na(thrombocytopenia_event_or_censor_date),
         !is.na(date_treatment_1st_line),
         !is.na(thrombocytopenia_binary),
         !is.na(treatment_type),
         !is.na(sex))

cox_uni_thrombocytopenia_sex <- coxph(            #Sex
  Surv(as.numeric(thrombocytopenia_event_or_censor_date - date_treatment_1st_line), thrombocytopenia_binary) ~ sex,
  data = cox_thrombocytopenia_data
)
summary(cox_uni_thrombocytopenia_sex)

cox_thrombocytopenia_data$age_group <- factor(cox_thrombocytopenia_data$age_group)

cox_uni_thrombocytopenia_age <- coxph(          #Age
  Surv(as.numeric(thrombocytopenia_event_or_censor_date - date_treatment_1st_line), thrombocytopenia_binary) ~ age_group, data = cox_thrombocytopenia_data)

summary(cox_uni_thrombocytopenia_age)


cox_thrombocytopenia_data$treatment_type <- factor(cox_thrombocytopenia_data$treatment_type)  #Treatment

cox_thrombocytopenia_data$treatment_type <- relevel(factor(cox_thrombocytopenia_data$treatment_type), ref = "btk-inhibitor")

cox_uni_thrombocytopenia_treatment <- coxph(      
  Surv(as.numeric(thrombocytopenia_event_or_censor_date - date_treatment_1st_line), thrombocytopenia_binary) ~ treatment_type, data = cox_thrombocytopenia_data)

summary(cox_uni_thrombocytopenia_treatment)

#Univariable analysis for IGHV, TP53_ab and del17p_at_treatment
#IGHV
cox_thrombocytopenia_data$IGHV_grouped <- relevel(factor(cox_thrombocytopenia_data$IGHV_grouped), ref = "Unmutated")
cox_uni_thrombocytopenia_IGHV <- coxph(            
  Surv(as.numeric(thrombocytopenia_event_or_censor_date - date_treatment_1st_line), thrombocytopenia_binary) ~ IGHV_grouped, data = cox_thrombocytopenia_data)
summary(cox_uni_thrombocytopenia_IGHV)

#TP53_ab
cox_thrombocytopenia_data$TP53_ab_grouped <- relevel(factor(cox_thrombocytopenia_data$TP53_ab_grouped), ref = "No")
cox_uni_thrombocytopenia_TP53 <- coxph(            
  Surv(as.numeric(thrombocytopenia_event_or_censor_date - date_treatment_1st_line), thrombocytopenia_binary) ~ TP53_ab_grouped, data = cox_thrombocytopenia_data)
summary(cox_uni_thrombocytopenia_TP53)

#del17p_at_treatment
cox_thrombocytopenia_data$del17p_at_treatment_grouped <- relevel(factor(cox_thrombocytopenia_data$del17p_at_treatment_grouped), ref = "No")
cox_uni_thrombocytopenia_del17p <- coxph(            
  Surv(as.numeric(thrombocytopenia_event_or_censor_date - date_treatment_1st_line), thrombocytopenia_binary) ~ del17p_at_treatment_grouped, data = cox_thrombocytopenia_data)
summary(cox_uni_thrombocytopenia_del17p)


#Leukopenia
cox_leukopenia_data <- cohort_done_statistics %>%
  filter(!is.na(leukopenia_event_or_censor_date),
         !is.na(date_treatment_1st_line),
         !is.na(leukopenia_binary),
         !is.na(treatment_type),
         !is.na(sex))

cox_uni_leukopenia_sex <- coxph(            #Sex
  Surv(as.numeric(leukopenia_event_or_censor_date - date_treatment_1st_line), leukopenia_binary) ~ sex,
  data = cox_leukopenia_data
)
summary(cox_uni_leukopenia_sex)

cox_leukopenia_data$age_group <- factor(cox_leukopenia_data$age_group)

cox_uni_leukopenia_age <- coxph(          #Age
  Surv(as.numeric(leukopenia_event_or_censor_date - date_treatment_1st_line), leukopenia_binary) ~ age_group, data = cox_leukopenia_data)

summary(cox_uni_leukopenia_age)


cox_leukopenia_data$treatment_type <- factor(cox_leukopenia_data$treatment_type)  #Treatment

cox_leukopenia_data$treatment_type <- relevel(factor(cox_leukopenia_data$treatment_type), ref = "btk-inhibitor")

cox_uni_leukopenia_treatment <- coxph(      
  Surv(as.numeric(leukopenia_event_or_censor_date - date_treatment_1st_line), leukopenia_binary) ~ treatment_type, data = cox_leukopenia_data)

summary(cox_uni_leukopenia_treatment)

#Univariable analysis for IGHV, TP53_ab and del17p_at_treatment
#IGHV
cox_leukopenia_data$IGHV_grouped <- relevel(factor(cox_leukopenia_data$IGHV_grouped), ref = "Unmutated")
cox_uni_leukopenia_IGHV <- coxph(            
  Surv(as.numeric(leukopenia_event_or_censor_date - date_treatment_1st_line), leukopenia_binary) ~ IGHV_grouped, data = cox_leukopenia_data)
summary(cox_uni_leukopenia_IGHV)

#TP53_ab
cox_leukopenia_data$TP53_ab_grouped <- relevel(factor(cox_leukopenia_data$TP53_ab_grouped), ref = "No")
cox_uni_leukopenia_TP53 <- coxph(            
  Surv(as.numeric(leukopenia_event_or_censor_date - date_treatment_1st_line), leukopenia_binary) ~ TP53_ab_grouped, data = cox_leukopenia_data)
summary(cox_uni_leukopenia_TP53)

#del17p_at_treatment
cox_leukopenia_data$del17p_at_treatment_grouped <- relevel(factor(cox_leukopenia_data$del17p_at_treatment_grouped), ref = "No")
cox_uni_leukopenia_del17p <- coxph(            
  Surv(as.numeric(leukopenia_event_or_censor_date - date_treatment_1st_line), leukopenia_binary) ~ del17p_at_treatment_grouped, data = cox_leukopenia_data)
summary(cox_uni_leukopenia_del17p)


#Creatinine increase (AKI)
cox_creatinine_data <- cohort_done_statistics %>%
  filter(!is.na(creatinine_event_or_censor_date),
         !is.na(date_treatment_1st_line),
         !is.na(creatinine_binary),
         !is.na(treatment_type),
         !is.na(sex))

cox_uni_creatinine_sex <- coxph(            #Sex
  Surv(as.numeric(creatinine_event_or_censor_date - date_treatment_1st_line), creatinine_binary) ~ sex,
  data = cox_creatinine_data
)
summary(cox_uni_creatinine_sex)

cox_creatinine_data$age_group <- factor(cox_creatinine_data$age_group)

cox_uni_creatinine_age <- coxph(          #Age
  Surv(as.numeric(creatinine_event_or_censor_date - date_treatment_1st_line), creatinine_binary) ~ age_group, data = cox_creatinine_data)

summary(cox_uni_creatinine_age)


cox_creatinine_data$treatment_type <- factor(cox_creatinine_data$treatment_type)  #Treatment

cox_creatinine_data$treatment_type <- relevel(factor(cox_creatinine_data$treatment_type), ref = "btk-inhibitor")

cox_uni_creatinine_treatment <- coxph(      
  Surv(as.numeric(creatinine_event_or_censor_date - date_treatment_1st_line), creatinine_binary) ~ treatment_type, data = cox_creatinine_data)

summary(cox_uni_creatinine_treatment)

#Univariable analysis for IGHV, TP53_ab and del17p_at_treatment
#IGHV
cox_creatinine_data$IGHV_grouped <- relevel(factor(cox_creatinine_data$IGHV_grouped), ref = "Unmutated")
cox_uni_creatinine_IGHV <- coxph(            
  Surv(as.numeric(creatinine_event_or_censor_date - date_treatment_1st_line), creatinine_binary) ~ IGHV_grouped, data = cox_creatinine_data)
summary(cox_uni_creatinine_IGHV)

#TP53_ab
cox_creatinine_data$TP53_ab_grouped <- relevel(factor(cox_creatinine_data$TP53_ab_grouped), ref = "No")
cox_uni_creatinine_TP53 <- coxph(            
  Surv(as.numeric(creatinine_event_or_censor_date - date_treatment_1st_line), creatinine_binary) ~ TP53_ab_grouped, data = cox_creatinine_data)
summary(cox_uni_creatinine_TP53)

#del17p_at_treatment
cox_creatinine_data$del17p_at_treatment_grouped <- relevel(factor(cox_creatinine_data$del17p_at_treatment_grouped), ref = "No")
cox_uni_creatinine_del17p <- coxph(            
  Surv(as.numeric(creatinine_event_or_censor_date - date_treatment_1st_line), creatinine_binary) ~ del17p_at_treatment_grouped, data = cox_creatinine_data)
summary(cox_uni_creatinine_del17p)


#### Cox proportional hazard models (multivariable) ####

#Prepping
class(cohort_done_statistics$age_group)
class(cohort_done_statistics$sex)
class(cohort_done_statistics$treatment_type)

cox_multi_data <- cohort_done_statistics %>%       #Changing sex, age_group and treatment_type to factor (necessary for Cox models)
  mutate(
    treatment_type = factor(treatment_type),
    age_group = factor(age_group),
    sex = factor(sex)
  )

#Anemia
cox_anemia_data <- cox_multi_data %>%
  filter(!is.na(anemia_event_or_censor_date),
         !is.na(date_treatment_1st_line),
         !is.na(anemia_binary),
         !is.na(age_group),
         !is.na(sex),
         !is.na(treatment_type))

cox_anemia_data$treatment_type <- relevel(cox_anemia_data$treatment_type, ref = "btk-inhibitor") #Making BTK-inhibitor the reference group for treatment types

cox_multi_anemia <- coxph(
  Surv(as.numeric(anemia_event_or_censor_date - date_treatment_1st_line), anemia_binary) ~
    age_group + sex + treatment_type,          #Adjusting for age, sex and treatment type
  data = cox_anemia_data
)
summary(cox_multi_anemia)


cox_anemia_data$IGHV_grouped <- relevel(factor(cox_anemia_data$IGHV_grouped), ref = "Unmutated")
cox_anemia_data$TP53_ab_grouped <- relevel(factor(cox_anemia_data$TP53_ab_grouped), ref = "No")
cox_anemia_data$del17p_at_treatment_grouped <- relevel(factor(cox_anemia_data$del17p_at_treatment_grouped), ref = "No")

cox_multi_anemia_all <- coxph(
  Surv(as.numeric(anemia_event_or_censor_date - date_treatment_1st_line), anemia_binary) ~
    age_group + sex + treatment_type + TP53_ab_grouped + del17p_at_treatment_grouped + IGHV_grouped, #Adjusting for age, sex, treatment type, IGHV, TP53 and del17p_at_treatment
  data = cox_anemia_data
)
summary(cox_multi_anemia_all)


#Neutropenia
cox_neutropenia_data <- cox_multi_data %>%
  filter(!is.na(neutropenia_event_or_censor_date),
         !is.na(date_treatment_1st_line),
         !is.na(neutropenia_binary),
         !is.na(age_group),
         !is.na(sex),
         !is.na(treatment_type))

cox_neutropenia_data$treatment_type <- relevel(cox_neutropenia_data$treatment_type, ref = "btk-inhibitor") #Making BTK-inhibitor the reference group for treatment types

cox_multi_neutropenia <- coxph(
  Surv(as.numeric(neutropenia_event_or_censor_date - date_treatment_1st_line), neutropenia_binary) ~
    age_group + sex + treatment_type,          #Adjusting for age, sex and treatment type
  data = cox_neutropenia_data
)
summary(cox_multi_neutropenia)

cox_neutropenia_data$IGHV_grouped <- relevel(factor(cox_neutropenia_data$IGHV_grouped), ref = "Unmutated")
cox_neutropenia_data$TP53_ab_grouped <- relevel(factor(cox_neutropenia_data$TP53_ab_grouped), ref = "No")
cox_neutropenia_data$del17p_at_treatment_grouped <- relevel(factor(cox_neutropenia_data$del17p_at_treatment_grouped), ref = "No")

cox_multi_neutropenia_all <- coxph(
  Surv(as.numeric(neutropenia_event_or_censor_date - date_treatment_1st_line), neutropenia_binary) ~
    age_group + sex + treatment_type + TP53_ab_grouped + del17p_at_treatment_grouped + IGHV_grouped, #Adjusting for age, sex, treatment type, IGHV, TP53 and del17p_at_treatment
  data = cox_neutropenia_data
)
summary(cox_multi_neutropenia_all)

#Lymphopenia
cox_lymphopenia_data <- cox_multi_data %>%
  filter(!is.na(lymphopenia_event_or_censor_date),
         !is.na(date_treatment_1st_line),
         !is.na(lymphopenia_binary),
         !is.na(age_group),
         !is.na(sex),
         !is.na(treatment_type))

cox_lymphopenia_data$treatment_type <- relevel(cox_lymphopenia_data$treatment_type, ref = "btk-inhibitor") #Making BTK-inhibitor the reference group for treatment types

cox_multi_lymphopenia <- coxph(
  Surv(as.numeric(lymphopenia_event_or_censor_date - date_treatment_1st_line), lymphopenia_binary) ~
    age_group + sex + treatment_type,          #Adjusting for age, sex and treatment type
  data = cox_lymphopenia_data
)
summary(cox_multi_lymphopenia)


cox_lymphopenia_data$IGHV_grouped <- relevel(factor(cox_lymphopenia_data$IGHV_grouped), ref = "Unmutated")
cox_lymphopenia_data$TP53_ab_grouped <- relevel(factor(cox_lymphopenia_data$TP53_ab_grouped), ref = "No")
cox_lymphopenia_data$del17p_at_treatment_grouped <- relevel(factor(cox_lymphopenia_data$del17p_at_treatment_grouped), ref = "No")

cox_multi_lymphopenia_all <- coxph(
  Surv(as.numeric(lymphopenia_event_or_censor_date - date_treatment_1st_line), lymphopenia_binary) ~
    age_group + sex + treatment_type + TP53_ab_grouped + del17p_at_treatment_grouped + IGHV_grouped, #Adjusting for age, sex, treatment type, IGHV, TP53 and del17p_at_treatment
  data = cox_lymphopenia_data
)
summary(cox_multi_lymphopenia_all)


#Thrombocytopenia
cox_thrombocytopenia_data <- cox_multi_data %>%
  filter(!is.na(thrombocytopenia_event_or_censor_date),
         !is.na(date_treatment_1st_line),
         !is.na(thrombocytopenia_binary),
         !is.na(age_group),
         !is.na(sex),
         !is.na(treatment_type))

cox_thrombocytopenia_data$treatment_type <- relevel(cox_thrombocytopenia_data$treatment_type, ref = "btk-inhibitor") #Making BTK-inhibitor the reference group for treatment types

cox_multi_thrombocytopenia <- coxph(
  Surv(as.numeric(thrombocytopenia_event_or_censor_date - date_treatment_1st_line), thrombocytopenia_binary) ~
    age_group + sex + treatment_type,          #Adjusting for age, sex and treatment type
  data = cox_thrombocytopenia_data
)
summary(cox_multi_thrombocytopenia)


cox_thrombocytopenia_data$IGHV_grouped <- relevel(factor(cox_thrombocytopenia_data$IGHV_grouped), ref = "Unmutated")
cox_thrombocytopenia_data$TP53_ab_grouped <- relevel(factor(cox_thrombocytopenia_data$TP53_ab_grouped), ref = "No")
cox_thrombocytopenia_data$del17p_at_treatment_grouped <- relevel(factor(cox_thrombocytopenia_data$del17p_at_treatment_grouped), ref = "No")

cox_multi_thrombocytopenia_all <- coxph(
  Surv(as.numeric(thrombocytopenia_event_or_censor_date - date_treatment_1st_line), thrombocytopenia_binary) ~
    age_group + sex + treatment_type + TP53_ab_grouped + del17p_at_treatment_grouped + IGHV_grouped, #Adjusting for age, sex, treatment type, IGHV, TP53 and del17p_at_treatment
  data = cox_thrombocytopenia_data
)
summary(cox_multi_thrombocytopenia_all)


#Leukopenia
cox_leukopenia_data <- cox_multi_data %>%
  filter(!is.na(leukopenia_event_or_censor_date),
         !is.na(date_treatment_1st_line),
         !is.na(leukopenia_binary),
         !is.na(age_group),
         !is.na(sex),
         !is.na(treatment_type))

cox_leukopenia_data$treatment_type <- relevel(cox_leukopenia_data$treatment_type, ref = "btk-inhibitor") #Making BTK-inhibitor the reference group for treatment types

cox_multi_leukopenia <- coxph(
  Surv(as.numeric(leukopenia_event_or_censor_date - date_treatment_1st_line), leukopenia_binary) ~
    age_group + sex + treatment_type,          #Adjusting for age, sex and treatment type
  data = cox_leukopenia_data
)
summary(cox_multi_leukopenia)


cox_leukopenia_data$IGHV_grouped <- relevel(factor(cox_leukopenia_data$IGHV_grouped), ref = "Unmutated")
cox_leukopenia_data$TP53_ab_grouped <- relevel(factor(cox_leukopenia_data$TP53_ab_grouped), ref = "No")
cox_leukopenia_data$del17p_at_treatment_grouped <- relevel(factor(cox_leukopenia_data$del17p_at_treatment_grouped), ref = "No")

cox_multi_leukopenia_all <- coxph(
  Surv(as.numeric(leukopenia_event_or_censor_date - date_treatment_1st_line), leukopenia_binary) ~
    age_group + sex + treatment_type + TP53_ab_grouped + del17p_at_treatment_grouped + IGHV_grouped, #Adjusting for age, sex, treatment type, IGHV, TP53 and del17p_at_treatment
  data = cox_leukopenia_data
)
summary(cox_multi_leukopenia_all)


#Creatinine increase (AKI)
cox_creatinine_data <- cox_multi_data %>%
  filter(!is.na(creatinine_event_or_censor_date),
         !is.na(date_treatment_1st_line),
         !is.na(creatinine_binary),
         !is.na(age_group),
         !is.na(sex),
         !is.na(treatment_type))

cox_creatinine_data$treatment_type <- relevel(cox_creatinine_data$treatment_type, ref = "btk-inhibitor") #Making BTK-inhibitor the reference group for treatment types

cox_multi_creatinine <- coxph(
  Surv(as.numeric(creatinine_event_or_censor_date - date_treatment_1st_line), creatinine_binary) ~
    age_group + sex + treatment_type,          #Adjusting for age, sex and treatment type
  data = cox_creatinine_data
)
summary(cox_multi_creatinine)

cox_creatinine_data$IGHV_grouped <- relevel(factor(cox_creatinine_data$IGHV_grouped), ref = "Unmutated")
cox_creatinine_data$TP53_ab_grouped <- relevel(factor(cox_creatinine_data$TP53_ab_grouped), ref = "No")
cox_creatinine_data$del17p_at_treatment_grouped <- relevel(factor(cox_creatinine_data$del17p_at_treatment_grouped), ref = "No")

cox_multi_creatinine_all <- coxph(
  Surv(as.numeric(creatinine_event_or_censor_date - date_treatment_1st_line), creatinine_binary) ~
    age_group + sex + treatment_type + TP53_ab_grouped + del17p_at_treatment_grouped + IGHV_grouped, #Adjusting for age, sex, treatment type, IGHV, TP53 and del17p_at_treatment
  data = cox_creatinine_data
)
summary(cox_multi_creatinine_all)



#### Forest plots (univariable) ####
#Making forest plots showing univariable hazard ratios (HR) with 95% confidence intervals for each adverse event, estimated using seperate cox proportianal hazard models. 
#Using data from above ("Cox proportional hazard models (univariable)")

library(broom)

#Helping function 
make_forest_univ <- function(models_names, title, file_out = NULL) {
  
  df <- imap_dfr(models_named, function(mod, lbl) {
    tidy(mod, exponentiate = TRUE, conf.int = TRUE) |>
      slice(1)|>
      transmute(
        label = lbl, 
        HR = estimate,
        lower = conf.low,
        upper = conf.high
      )
  })|>
    arrange(HR) |>
    mutate(label = factor(label, levels = rev(label)))
  
p <- ggplot(df, aes(x=label, y=HR, ymin=lower, ymax=upper))+
  geom_hline(yintercept = 1, linetype = 2)+
  geom_pointrange(size = 0.9)+
  coord_flip()+
  scale_y_log10()+
  labs(
    title = title, 
    x = NULL,
    y = "Hazard ratio (95% CI, log scale)"
  )+
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.y = element_blank(),
        plot.title.position = "plot")

if(!is.null(file_out)) ggsave(file_out, p, width = 6.5, height = 4.2, dpi = 300)
p
}

#Sex: univariable forest plots for all adverse events (ref. female)
sex_fp_uni_models <- list(
  Anemia = cox_uni_anemia_sex,
  Thrombocytopenia = cox_uni_thrombocytopenia_sex,
  Leukopenia = cox_uni_leukopenia_sex,
  Neutropenia = cox_uni_neutropenia_sex,
  Lymphopenia = cox_uni_lymphopenia_sex,
  "Acute kidney injury" = cox_uni_creatinine_sex
)

sex_df <- imap_dfr(sex_fp_uni_models, ~ tidy(.x, exponentiate = TRUE, conf.int = TRUE) %>%
                     slice(1)%>%
                     transmute(
                       AE = .y, 
                       HR = estimate, 
                       lower = conf.low, 
                       upper = conf.high
                     ))

ggplot(sex_df, aes(x = AE, y = HR, ymin = lower, ymax = upper))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  geom_pointrange(size = 0.9)+
  coord_flip()+
  scale_y_log10()+
  labs(
    title = "Univariable hazard ratios by adverse event",
    subtitle = "Covariate: Sex (male vs. female)",
    x = NULL, y = "Hazard ratio (95% CI)"
  )+
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.y = element_blank(),
        plot.title.position = "plot",
        axis.title = element_text(size=16),
        axis.text = element_text(size=14))

#Age: univariable forest plots for all adverse events (ref. >60)
age_fp_uni_models <- list(
  Anemia = cox_uni_anemia_age,
  Thrombocytopenia = cox_uni_thrombocytopenia_age,
  Leukopenia = cox_uni_leukopenia_age,
  Neutropenia = cox_uni_neutropenia_age,
  Lymphopenia = cox_uni_lymphopenia_age,
  "Acute kidney injury" = cox_uni_creatinine_age
)

age_df <- imap_dfr(age_fp_uni_models, ~ tidy(.x, exponentiate = TRUE, conf.int = TRUE) %>%
                     slice(1)%>%
                     transmute(
                       AE = .y, 
                       HR = estimate, 
                       lower = conf.low, 
                       upper = conf.high
                     ))

ggplot(age_df, aes(x = AE, y = HR, ymin = lower, ymax = upper))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  geom_pointrange(size = 0.9)+
  coord_flip()+
  scale_y_log10()+
  labs(
    title = "Univariable hazard ratios by adverse event",
    subtitle = "Covariate: Age (ref. >60)",
    x = NULL, y = "Hazard ratio (95% CI)"
  )+
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.y = element_blank(),
        plot.title.position = "plot",
        axis.title = element_text(size=16),
        axis.text = element_text(size=14))

#Treatment type: univariable forest plots for all adverse events (ref. BTK-inhibitor)
treatment_fp_uni_models <- list(
  Anemia = cox_uni_anemia_treatment,
  Thrombocytopenia = cox_uni_thrombocytopenia_treatment,
  Leukopenia = cox_uni_leukopenia_treatment,
  Neutropenia = cox_uni_neutropenia_treatment,
  Lymphopenia = cox_uni_lymphopenia_treatment,
  "Acute kidney injury" = cox_uni_creatinine_treatment
)

treatment_df <- imap_dfr(treatment_fp_uni_models, ~ tidy(.x, exponentiate = TRUE, conf.int = TRUE) %>%
                     slice(1)%>%
                     transmute(
                       AE = .y, 
                       HR = estimate, 
                       lower = conf.low, 
                       upper = conf.high
                     ))

ggplot(treatment_df, aes(x = AE, y = HR, ymin = lower, ymax = upper))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  geom_pointrange(size = 0.9)+
  coord_flip()+
  scale_y_log10()+
  labs(
    title = "Univariable hazard ratios by adverse event",
    subtitle = "Covariate: Treatment type (ref. BTK-inhibitor)",
    x = NULL, y = "Hazard ratio (95% CI)"
  )+
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.y = element_blank(),
        plot.title.position = "plot")


#### Forest plots (multivariable) ####
#Making forest plots showing multivariable hazard ratios (HR) with 95% confidence intervals for each adverse event, estimated using seperate cox proportianal hazard models. 
#Using data from above ("Cox proportional hazard models (multivariable)")


#Sex: univariable forest plots for all adverse events (ref. female)
fp_multi_models <- list(
  Anemia = cox_multi_anemia_all,
  Thrombocytopenia = cox_multi_thrombocytopenia_all,
  Leukopenia = cox_multi_leukopenia_all,
  Neutropenia = cox_multi_neutropenia_all,
  Lymphopenia = cox_multi_lymphopenia_all,
  "Acute kidney injury" = cox_multi_lymphopenia_all
)%>%
  discard(is.null)

#Samler al data i ét frame
tidy_all <- imap_dfr(fp_multi_models, ~{
  tidy(.x, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(event = .y)
  })

#Cleaning in rows and labels
tidy_clean <- tidy_all %>%
  filter(term !="(Intercept)")%>%
  mutate(
    term_nice = term %>%
      str_replace_all("_grouped", "") %>%
      str_replace_all("_at_treatment", "") %>%
      str_replace_all("treatment_type", "Treatment") %>%
      str_replace_all("TP53_ab", "TP53 abnormality") %>%
      str_replace_all("sex", "Sex") %>%
      str_replace_all("age_group", "Age Group") %>%
      str_replace_all("age", "Age"),
    event = factor(event, levels = names(fp_multi_models))
  )%>%
  arrange(event, desc(estimate))

tidy_clean <- tidy_clean %>%
  mutate(
    term_nice = factor(
      term_nice, 
      levels = c(
        "Treatment", 
        "Age",
        "Sex",
        "IGHV",
        "TP53 abnormality",
        "del17p"
      )
    )
  )

#Forest plot
ggplot(tidy_clean, aes(y=term_nice, x=estimate, xmin=conf.low, xmax=conf.high))+
  geom_vline(xintercept = 1, linetype = 2)+
  geom_pointrange(size = 1.1)+
  geom_point(size=3)+
  facet_wrap(~event, scales = "free_y")+
  labs(
    title = "Multivariable hazard ratios (95% CI) for adverse events",
    x = "Hazard ratio (log scale)",
    y = NULL
  )+
  scale_x_log10()+
  theme_minimal(base_size = 18)+
  theme(
    plot.title = element_text(size=26, face="bold", hjust=0,5),
    strip.text = element_text(size=20, face="bold"),
    axis.title.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=16),
    panel.spacing = unit(1.2, "lines"),
    plot.margin = margin(20, 20, 20, 20)
  )



ggplot(tidy_all, aes(y=term_nice, x=estimate, xmin=conf.low, xmax=conf.high))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  coord_flip()+
  scale_y_log10()+
  labs(
    title = "Multivariable hazard ratios by adverse event",
    x = NULL, y = "Hazard ratio (95% CI)"
  )+
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.y = element_blank(),
        plot.title.position = "plot")
