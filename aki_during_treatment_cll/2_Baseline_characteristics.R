# This script adds baseline characteristics to the CLL and SLL cohort
# date: "2025-02-10"
# author: christian brieghel og emma siig
# comment: run after 1cohort.R or downloade table: 1cohort_clean.rds

source('/ngc/projects2/dalyca_r/clean_r/load_dalycare_package.R')

cohort_clean <- readRDS("/ngc/projects2/dalyca_r/emmsii_r/sandbox/tables/1cohort_clean.rds")
cohort_riskfactors <- cohort_clean

#### CLL_TREAT - ready ####
load_dataset("CLL_TREAT", cohort_riskfactors$patientid)

#### Smoking history - ready ####
load_dataset("SP_Social_Hx", cohort_riskfactors$patientid)
SP_Social_Hx_subset %>% head2()
SP_Social_Hx_subset$ryger %>% table(exclude = NULL)
SP_Social_Hx_subset %>% nrow_npatients()

#Smoking and drinking history: 
history_smoking_drinking <- SP_Social_Hx_subset %>% select(patientid, timestamp = registringsdato, smoking_history = ryger, drinking_history = drikker) 

# Filtrer mistakes on date
mistake <- history_smoking_drinking %>%
  filter(timestamp < as.Date("2000-01-01")) %>%
  mutate(seconds = as.numeric(difftime(timestamp, as.Date("1970-01-01"), units = "secs")))

# New dates based on seconds
mistake <- mistake %>%
  mutate(date = timestamp + seconds * ddays(1),
         date = as.Date(date))

# Merge back on data
history_smoking_drinking <- history_smoking_drinking %>%
  left_join(mistake %>% select(patientid, timestamp, date), by = c("patientid", "timestamp")) %>%
  mutate(timestamp = as.Date(timestamp))

# Update timestamp
history_smoking_drinking <- history_smoking_drinking %>%
  mutate(
    timestamp = case_when(
      !is.na(date) ~ date,    
      TRUE ~ timestamp       
    )) 

history_smoking_drinking <- history_smoking_drinking %>%
  select(-date)

load_dataset("CLL_TREAT", cohort_riskfactors$patientid)
CLL_TREAT_subset_smoking <- CLL_TREAT_subset %>% 
  select(patientid, Smoking)

CLL_TREAT_subset_smoking <- CLL_TREAT_subset_smoking %>%
  mutate(Smoking = case_when(
    Smoking == 0 ~ "Aldrig",
    Smoking == 1 ~ "Er holdt op",
    Smoking == 2 ~ "Ja",
    TRUE ~ NA_character_  
  ))

#Adds smoking status from CLL_TREAT to history_smoking_drinking, if there isn't already is a datapoint
CLL_TREAT_subset_smoking$patientid <- as.numeric(CLL_TREAT_subset_smoking$patientid) #Makes patientid numeric 
history_smoking_drinking$patientid <- as.numeric(history_smoking_drinking$patientid)

CLL_TREAT_subset_smoking <- CLL_TREAT_subset_smoking %>% #Filters CLL_TREAT for patients there isn't in history_smoking_drinking
  filter(!patientid %in% history_smoking_drinking$patientid)

history_smoking_drinking <- bind_rows(history_smoking_drinking, CLL_TREAT_subset_smoking) #Adds the (extra) patients from CLL_TREAT to history_smoking_drinking

history_smoking_drinking <- history_smoking_drinking %>%
  mutate(ryger = coalesce(smoking_history, Smoking)) %>%  #Mutates smoking_history and Smoking, so values from Smoking is inserted if smoking_history is NA
  select(-Smoking)  #Deletes the column Smoking

history_smoking_drinking <- history_smoking_drinking %>% #Making "Ikke spurgt" into NA's
  mutate(smoking_history = case_when(
    smoking_history == "Ikke spurgt" ~ NA_character_,
    TRUE ~ smoking_history
  ))

history_smoking_drinking <- history_smoking_drinking %>%
  left_join(cohort_riskfactors %>% select(patientid, date_treatment_1st_line), by = "patientid") #Adds date_treatment_1st_line from cohort_riskfactors into history_smoking_drinking

history_smoking_drinking <- history_smoking_drinking %>%
  filter(timestamp <= date_treatment_1st_line + 14) %>% #Only keeps measures before or max 2 weeks after first treatment date
  group_by(patientid) %>%  
  slice_min(abs(timestamp - date_treatment_1st_line), n = 1) %>% #Keeps the date closest to the treatment date
  ungroup()

history_smoking_drinking <- history_smoking_drinking %>% distinct() #Removes duplicates

history_smoking_drinking$diff_days <- as.numeric(history_smoking_drinking$timestamp - history_smoking_drinking$date_treatment_1st_line) #Making sure the code above worked (regarding dates)
hist(history_smoking_drinking$diff_days)
summary(history_smoking_drinking$diff_days)

history_smoking_drinking <- dplyr::rename(history_smoking_drinking, date_smoking_drinking = timestamp) #Renames timestamp to date_smoking_drinking

history_smoking_drinking %>% nrow_npatients()

cohort_riskfactors <- cohort_riskfactors %>%
  left_join(history_smoking_drinking %>% select(patientid, smoking_history, date_smoking_drinking), by = "patientid") #Inserts smoking_history from history_smoking_drinking into our cohort_riskfactors

history_smoking_drinking %>%
  count(smoking_history)

#### Drinking history - not using ####
#Kør rygestatus inden

#SP_Social_Hx_subset$drikker %>% table(exclude = NULL)
#SP_Social_Hx_subset %>% nrow_npatients()

#cohort_riskfactors <- cohort_riskfactors %>%
#  left_join(history_smoking_drinking %>% select(patientid, drinking_history), by = "patientid") #Rykker drinking_history fra history_smoking_drinking ind i cohort

#### BMI - ready ####
load_dataset("SP_VitaleVaerdier", cohort_riskfactors$patientid)
SP_VitaleVaerdier_subset %>%  head2()
SP_VitaleVaerdier_subset$displayname %>% table()
BMI_data <- SP_VitaleVaerdier_subset %>% BMI()
BMI_data %>% nrow_npatients()

BMI_data <- BMI_data %>%
  left_join(cohort_riskfactors %>% select(patientid, date_treatment_1st_line), by = "patientid") #Inserts date_treatment_1st_line from cohort_riskfactors into BMI_data

BMI_data <- BMI_data %>% #Makes it into date format
  mutate(
    date = as.Date(date),
    date_treatment_1st_line = as.Date(date_treatment_1st_line)
  )

BMI_data <- BMI_data %>%
  filter(date >= date_treatment_1st_line - 365 & date <= date_treatment_1st_line + 14) %>% #Only keeps measurements 1 year before or max 2 weeks after first date of treatment
  group_by(patientid) %>%
  slice_min(abs(difftime(date, date_treatment_1st_line, units = "days")), n = 1, with_ties = FALSE) %>% #Chooses the date closest to the first treatment date
  ungroup()

BMI_data %>% nrow_npatients()

BMI_data$diff_days <- as.numeric(BMI_data$date - BMI_data$date_treatment_1st_line) #Makes sure the codes worked
hist(BMI_data$diff_days)
summary(BMI_data$diff_days) #YES :D

BMI_data <- BMI_data %>% distinct() #Removes duplicates

BMI_data <- dplyr::rename(BMI_data, date_BMI = date) #Renames date to date_BMI

BMI_data %>% nrow_npatients()

cohort_riskfactors <- cohort_riskfactors %>%
  left_join(BMI_data %>% select(patientid, BMI, date_BMI), by = "patientid") #Inserts BMI and date_BMI from BMI_data into cohort_riskfactors

hist(cohort_riskfactors$BMI) #One patient with BMI on 144? 
summary(cohort_riskfactors$BMI)
cohort_riskfactors <- cohort_riskfactors %>% #Removes data on BMI over 100 
  mutate(BMI = ifelse(BMI > 100, NA, BMI))

#### Comorbidities ####
load_dataset()
load_dataset("diagnoses_all", cohort_riskfactors$patientid)
diagnoses_all_subset %>%  head2()
diagnoses_all_subset$diagnosis %>% table()

CCI_data <- diagnoses_all_subset %>% 
  left_join(cohort_riskfactors %>% select(patientid, date_treatment_CLL_SLL = date_treatment_1st_line)) %>% 
  filter(date_diagnosis < date_treatment_CLL_SLL)%>% #Filters all diagnoses, so we only keep diagnoses before the first treatment of CLL/SLL
  CCI(exclude_CLL_score = TRUE) #Excludes the CLL-score in the overall CCI-score

ggplot(CCI_data, aes(x = CCI.2011.update)) + #CCI score without CLL-score
  geom_bar() +
  labs(x = "CCI.2011.update", y = "Patients") +
  theme_minimal()

CCI_data %>% nrow_npatients()

#table(CCI_data$CCI.2011.update, exclude = if(CLL = No))

CCI_data <- CCI_data %>% #categorizes CCI.2011.update score so the highest category is grouped as 3+
  mutate(
    CCI.2011.update_cat = case_when(
      CCI.2011.update == 0 ~ "0",
      CCI.2011.update == 1 ~ "1",
      CCI.2011.update == 2 ~ "2",
      CCI.2011.update >= 3 ~ "3+",
      TRUE ~ NA_character_
    ),
    CCI.2011.update_cat = factor(CCI.2011.update_cat, levels = c("0", "1", "2", "3+"), ordered = TRUE) #Makes it a factor instead of a character 
  )

ggplot(CCI_data, aes(x = CCI.2011.update_cat)) + #Visuals of patients and CCI score
  geom_bar(fill = "#2C7BB6") +
  labs(
    title = "Fordeling af patienter efter CCI-score",
    x = "CCI-score (kategoriseret)",
    y = "Antal patienter"
  ) +
  theme_minimal()

#Ready to be used, when approved by Emelie
#cohort_riskfactors <- cohort_riskfactors %>% #Inserts it in cohort
#  left_join(
#    CCI_data %>% select(patientid, CCI.2011.update_cat), 
#    by = "patientid"
#    )


#### FISH - ready ####
RKKP_CLL %>% names()
RKKP_CLL$Beh_FISH_TP53
RKKP_CLL$Beh_Del17p
utable(~del17p+del13q+del11q+tri12,CLL_clean) #at time of diagnosis
utable(~del17p_at_treatment+TP53_at_treatment,CLL_clean) #at time of treatment

CLL_clean %>% names()
CLL_clean$IPI_score_del17p_only

#Making FISH subset
CLL_FISH_subset <- CLL_clean %>%
  semi_join(cohort_riskfactors, by = "patientid") %>%
  select(patientid, del17p_at_treatment, TP53_at_treatment, del13q, del11q, del17p, tri12) #Variables from CLL_clean

CLL_FISH_subset <- CLL_FISH_subset %>%
  left_join(CLL_TREAT_subset %>% select(patientid, FISH_diag), by = "patientid") #Variable from CCL_TREAT

CLL_FISH_subset <- CLL_FISH_subset %>%
  left_join(RKKP_CLL_CLEAN %>% select(patientid, FISH.status), by = "patientid") #Variable from RKKP

nrow_npatients(CLL_FISH_subset) #Checking number of rows and patients
summary(CLL_FISH_subset)
CLL_FISH_subset %>%
  count(FISH.status, sort = TRUE)

#Making FISH_status
CLL_FISH_subset <- dplyr::rename(CLL_FISH_subset, FISH_status = `FISH.status`) #Changing name from FISH.status to FISH_status

CLL_FISH_subset_status <- CLL_FISH_subset %>% #Inserts info on del13q, del11q, del17p and tri12 if FISH_status=NA
  select(patientid, FISH_status, del13q, del11q, del17p, tri12, FISH_diag) %>% 
  mutate(
      FISH_status = if_else(
      is.na(FISH_status),
      paste(
        ifelse(del13q == "Yes", "Del13q", NA),
        ifelse(del11q == "Yes", "Del11q", NA),
        ifelse(del17p == "Yes", "Del17p", NA),
        ifelse(tri12  == "Yes", "Tri12",  NA),
        sep = "-"
      ),
      FISH_status
    )
  ) %>% 
  mutate( #Removes excess - and " " and more
    FISH_status = gsub("NA-|NA|-$| |- ", "", FISH_status), 
    FISH_status = gsub("NA-|NA|-$| |- ", "", FISH_status),
    FISH_status = if_else(FISH_status == "", NA, FISH_status)
  )

CLL_FISH_subset_status %>%          #Checking the code above worked and checking the new FISH_status groups (eg. Del13q-Del11q)
  count(FISH_status, sort = TRUE)

CLL_FISH_subset_status <- CLL_FISH_subset_status %>% #Simplifies FISH_status variable when multiple mutations occurs - by only keeping the highest mutation in the Döner hierarchy
  mutate(
    FISH_status = case_when(
      str_detect(FISH_status, "Del17p") ~ "Del17p",
      str_detect(FISH_status, "Del11q") ~ "Del11q",
      str_detect(FISH_status, "Tri12")  ~ "Tri12",
      str_detect(FISH_status, "Del13q") ~ "Del13q",
      TRUE ~ FISH_status
    )
  )

CLL_FISH_subset_status <- CLL_FISH_subset_status %>%  #Inserts "Normal" if del13q, del11q, del17p and tri12 all = "No"
  mutate(
    FISH_status = if_else(
      is.na(FISH_status) &
        del13q == "No" & del11q == "No" & del17p == "No" & tri12 == "No",
      "Normal",
      FISH_status
    )
  )

#Inserts info from FISH_diag (from CLL_TREAT) when FISH_status = NA
CLL_FISH_subset_status %>%          #Checking FISH_diag
  count(FISH_diag, sort = TRUE)

CLL_FISH_subset_status <- CLL_FISH_subset_status %>% #Cleaning FISH_diag
  mutate(
    FISH_diag_clean = case_when(
      str_detect(FISH_diag, "Missing") ~ NA,
      str_detect(FISH_diag, "Udført,\\s*[Nn]ormal") ~ "Normal",
      TRUE ~ str_remove(FISH_diag, "Udført,\\s*")
    ),
    
    FISH_status = if_else( #Insert info from FISH_diag when FISH_status = NA
      is.na(FISH_status) & !is.na(FISH_diag_clean),
      FISH_diag_clean,
      FISH_status
    )
  ) %>%
  select(-FISH_diag_clean) #Deletes extra variable

CLL_FISH_subset_status %>%          #Checking FISH_status
  count(FISH_status, sort = TRUE)
#FISH_status is done now

#Inserts FISH_status as FISH_status_all in CLL_FISH_subset
CLL_FISH_subset <- CLL_FISH_subset %>%
  left_join(CLL_FISH_subset_status %>% select(patientid, FISH_status) %>%
              dplyr::rename(FISH_status_all = FISH_status),
            by = "patientid")

CLL_FISH_subset %>%          #Checking del17p_at_treatment
  count(del17p_at_treatment, sort = TRUE)

CLL_FISH_subset <- CLL_FISH_subset %>%      #Makes UNK in del17p_at_treatment into NA
  mutate(del17p_at_treatment = ifelse(del17p_at_treatment == "UNK", NA_character_, as.character(del17p_at_treatment)))


CLL_FISH_subset <- CLL_FISH_subset %>% #Inserts "Yes" if del17p_at_treatment is NA and FISH_status_all has "Del17" in it
  mutate(
    del17p_at_treatment = if_else(
      is.na(del17p_at_treatment) & str_detect(FISH_status_all, "Del17p"),
      "Yes",
      del17p_at_treatment
    )
  )

#Insert in cohort
cohort_riskfactors <- cohort_riskfactors %>%
  mutate(patientid = as.numeric(patientid))

CLL_FISH_subset <- CLL_FISH_subset %>%
  mutate(patientid = as.numeric(patientid))

cohort_riskfactors <- cohort_riskfactors %>%
  left_join(CLL_FISH_subset %>% select(patientid, del17p_at_treatment, FISH_status_all), by = "patientid") #Inserts del17p_at_treatment and FISH_status_all from CLL_FISH_subset into cohort_riskfactors


#### IGHV - ready ####
cohort_riskfactors %>% summary 
count(CLL_TREAT_subset_IGHV, IGHV)
count(cohort_riskfactors, IGHV)
class(CLL_TREAT_subset$patientid)
class(cohort_riskfactors$patientid)
summary(CLL_TREAT_subset$IGHV)

CLL_TREAT_subset$patientid <- as.numeric(CLL_TREAT_subset$patientid)

CLL_TREAT_subset_IGHV <- CLL_TREAT_subset %>% #Renames "missing" to NA
  mutate(IGHV = na_if(IGHV, "Missing"))

cohort_riskfactors <- merge(cohort_riskfactors, CLL_TREAT_subset_IGHV[, c("patientid", "IGHV")], #Merge IGHV from CLL_TREAT into cohort by patientid
                by = "patientid", all.x = TRUE, suffixes = c("", "_cll_treat"))

cohort_riskfactors$IGHV[is.na(cohort_riskfactors$IGHV) & !is.na(cohort_riskfactors$IGHV_cll)] <- cohort_riskfactors$IGHV_cll[is.na(cohort_riskfactors$IGHV) & !is.na(cohort_riskfactors$IGHV_cll)] #Updates IGHV in cohort, if it is NA and a value is in CLL_TREAT

cohort_riskfactors <- cohort_riskfactors %>% #Deletes extra variable 
  select(-IGHV_cll_treat)

#### LSR - not using ####
SDS_DATASETS
load_dataset("SDS_epikur", cohort_riskfactors$patientid)
head2(SDS_epikur_subset)
hypertension_data <- SDS_epikur_subset %>% ATC_hypertensives() #hypertension
ab_data <- SDS_epikur_subset %>% ATC_AB() #antibiotics


#### Solving date mistake in SP_Social ####
#This code also runs in the "Smoking history - ready"-section

#Smoking and drinking history: 
load_dataset("SP_Social_Hx", feature_tabel$patientid)

history_smoking_drinking <- SP_Social_Hx_subset %>% select(patientid, timestamp = registringsdato, smoking_history = ryger, drinking_history = drikker) 

# Filtrer mistakes on date
mistake <- history_smoking_drinking %>%
  filter(timestamp < as.Date("2000-01-01")) %>%
  mutate(seconds = as.numeric(difftime(timestamp, as.Date("1970-01-01"), units = "secs")))

# New dates based on seconds
mistake <- mistake %>%
  mutate(date = timestamp + seconds * ddays(1),
         date = as.Date(date))


# Merge back on data
history_smoking_drinking <- history_smoking_drinking %>%
  left_join(mistake %>% select(patientid, timestamp, date), by = c("patientid", "timestamp")) %>%
  mutate(timestamp = as.Date(timestamp))

# Update timestamp and slice for latest value
history_smoking_drinking <- history_smoking_drinking %>%
  mutate(
    timestamp = case_when(
      !is.na(date) ~ date,    
      TRUE ~ timestamp       
    )) %>%
  select(-date) %>%
  left_join(RKKP_dates %>% select(patientid, date_flt), by = "patientid") %>%
  group_by(patientid) %>%
  arrange(patientid, timestamp) %>%
  filter(timestamp < date_flt) %>%
  slice_max(timestamp) %>%
  slice(1) %>%
  mutate(drinking_history = case_when(
    drinking_history == "Ja" ~ 1,         
    drinking_history == "Nej" ~ 0,        
    drinking_history == "Aldrig" ~ 0,     
    drinking_history == "Ikke aktuelt" ~ 0,  
    drinking_history == "Ikke spurgt" ~ NA,
    drinking_history == "Udskyd" ~ NA,
    TRUE ~ NA                           
  )) %>%
  mutate(smoking_history = case_when(
    smoking_history == "Aldrig" ~ 0,         
    smoking_history == "Er holdt op" ~ 1,   #Earlier and passive smoking in same 
    smoking_history == "Ikke spurgt" ~ 0,    
    smoking_history == "Passive" ~ 1,        
    smoking_history == "Ja" ~ 2,             
    TRUE ~ NA                               
  ))


#### SAVE TABLES ####
saveRDS(cohort_riskfactors, "/ngc/projects2/dalyca_r/emmsii_r/sandbox/tables/2cohort_riskfactors.rds")

