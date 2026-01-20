# This script identifies patients with CLL/SLL, and adds characteristics
# date: "2025-02-10"
# author: christian brieghel and emma siig

source('/ngc/projects2/dalyca_r/clean_r/load_dalycare_package.R')
load_dataset(c('patient', 't_dalycare_diagnoses', 'RKKP_CLL', 'RKKP_CLL_CLEAN', 'RKKP_LYFO'))

# CLL = t_dalycare_diagnoses %>%
#   filter_first_diagnosis(ICD10.CLL, str_contains = FALSE)

CLL_clean = RKKP_CLL %>% 
  clean_RKKP_CLL()
SLL_clean = RKKP_LYFO %>% 
  clean_RKKP_LYFO() %>% 
  filter(subtype == 'SLL')

CLL_clean$treatment_type
CLL_clean %>% nrow_npatients()
CLL_clean$TP53_ab
utable(~age+binet+B2M+IGHV+TP53_ab, CLL_clean) 
SLL_clean %>% nrow_npatients()
SLL_clean$date_treatment_1st_line %>% summary
SLL_clean$B2M_diagnosis
SLL_clean$AA_stage_diagnosis
utable(~age_diagnosis+AA_stage_diagnosis+B2M_diagnosis, SLL_clean) 
SLL_clean$chemo_regime_1_type_1st_line %>% table
SLL_clean$immunotherapy_type_1st_line
CLL_clean$treatment_type

#### The making of the cohort ####
#Making a new dataset (cohort), that combines data from CLL_clean and SLL_clean
cohort = bind_rows(CLL_clean %>% transmute(patientid, date_diagnosis, hospital_id, date_treatment_1st_line, treatment_type,
                                           age, stage = binet, B2M, IGHV, TP53_ab, type = 'CLL'),
                   SLL_clean %>% transmute(patientid, date_diagnosis, hospital_id, date_treatment_1st_line, 
                                           treatment_type = paste0(chemo_regime_1_type_1st_line, '-', immunotherapy_type_1st_line),
                                           age = age_diagnosis,
                                           stage = factor(AA_stage_diagnosis, levels = c(1:4)), 
                                           B2M = ifelse(B2M_diagnosis < 4.0, '<4.0 mg/L', '>4.0 mg/L'),
                                           IGHV = NA,
                                           TP53_ab = NA,
                                           type = 'SLL')) %>% 
  filter(!is.na(date_treatment_1st_line)) %>% 
  group_by(patientid) %>% 
  arrange(date_diagnosis) %>% 
  slice(1) %>% 
  ungroup() %>% 
  left_join(patient)
cohort$treatment_type %>% table # needs filtering!!
cohort %>% nrow_npatients()

load_dataset('LAB_IGHVIMGT', cohort$patientid)
load_npu_common()
B2M_df = load_biochemistry(c(NPU.B2M)) %>% 
  filter(patientid %in% cohort$patientid) %>% 
  clean_lab_values()

#Making cohort_2
cohort_2 = cohort %>%
  left_join(LAB_IGHVIMGT_subset %>% select(patientid, date_ighv = date_sample, ighv_lab = ighv), 'patientid') %>% 
  mutate(time_ighv = diff_days(date_diagnosis, date_ighv),
         time_ighv = ifelse(time_ighv < 90 ,time_ighv, NA),
         ighv_lab = factor(ifelse(is.na(time_ighv), NA, ighv_lab), c('Mutated', 'Unmutated')),
         IGHV = if_else(is.na(IGHV), ighv_lab, IGHV)) %>% 
  left_join(B2M_df %>% 
              transmute(patientid,
                        date_B2M = as_date(samplingdate), 
                        value =  as.numeric(value2))) %>%
  slice_closest_value(date_baseline = date_diagnosis, date_value = date_B2M, value = value) %>% 
  mutate(B2M_value = ifelse(value < 4.0, '<4.0 mg/L', '>4.0 mg/L'),
         B2M = if_else(is.na(B2M), B2M_value, B2M),
         B2M = factor(B2M, c('<4.0 mg/L', '>4.0 mg/L')))

cohort_2 %>% nrow_npatients() # OK!
write_csv2(cohort_2, '/ngc/projects2/dalyca_r/chribr_r/EMELIE/CLL_AE_EMMA/data/cohort.csv')

#### Time from diagnosis to first treatment ####
cohort$date_diagnosis <- as.Date(cohort$date_diagnosis) #converts to date-format 
cohort$date_treatment_1st_line <- as.Date(cohort$date_treatment_1st_line) #converts to date-format
cohort$time_to_treatment_days <- as.numeric(cohort$date_treatment_1st_line - cohort$date_diagnosis) #converts to number of days and finds the difference
mean_time <- mean(cohort$time_to_treatment_days, na.rm = TRUE) #calculates mean time (days) from diagnosis to treatment
print(mean_time/365) #mean time in year
median_time <- median(cohort$time_to_treatment_days, na.rm = TRUE) #calculates median time from diagnosis to treatment
print(median_time/365) #median (years)

#### Mean and median age for first treatment ####
cohort$age_at_treatment <- as.numeric(difftime(cohort$date_treatment_1st_line, 
                                                 cohort$date_birth, 
                                                 units = "days")) / 365.25 #Calculates age at treatment
cohort <- cohort %>%
  mutate(age_at_treatment = floor(age_at_treatment))
mean_age <- mean(cohort$age_at_treatment, na.rm = TRUE) #Mean age
print(mean_age)
utable(~Q(age_at_treatment), cohort) #Median age at first treatment
hist(cohort$age_at_treatment)

#### Cleaning of treatment type ####
sum(is.na(cohort$treatment_type)) #Number of NA's in treatment_type
mean(is.na(cohort$treatment_type)) * 100 #%NA's

cohort$treatment_type %>% table #Needs cleaning - same type is written in different ways
unique(cohort$treatment_type) #97 different treatment types
cohort_clean <- cohort %>%
  mutate(treatment_type = trimws(tolower(treatment_type)))  # Removes space "" and makes it only with small letters

cohort_clean %>% 
  mutate(
    treatment_type = gsub("-na|na-|na|no_chemo-|-no_chemo|no_chemo", "", treatment_type),
    treatment_type = if_else(treatment_type == '', NA, treatment_type)
  ) 
  
cohort_clean <- cohort_clean %>%
  mutate(treatment_type = gsub("^no_chemo[-]*|[-]na[-]*|na", "", treatment_type),
         treatment_type = gsub("cvp", "cop", treatment_type),
         treatment_type = ifelse(treatment_type == "", NA, treatment_type))

#Adjusts spelling
   cohort_clean <- cohort_clean %>%
  mutate(treatment_type = case_when(
    treatment_type %in% c("bendamustine-rituximab", "bendamustin-rituximab") ~ "bendamustine-rituximab",
    treatment_type %in% c("fludara-rituximab", "fludarabin-rituximab") ~ "fludarabin-rituximab",
    treatment_type %in% c("rituximab", "-rituximab") ~ "rituximab",
    treatment_type %in% c("bendamustin", "bendamustine") ~ "bendamustine",
      TRUE ~ treatment_type))
  
   unique(cohort_clean$treatment_type) #83 different
   cohort_clean$treatment_type %>% table 
   
   cohort_clean <- cohort_clean %>%   #Changes "other" to "other_imm" except when "other_chemo"
     mutate(treatment_type = gsub("\\bother\\b", "other_imm", treatment_type))

   cohort_clean <- cohort_clean %>% 
     mutate(treatment_type = case_when(
       treatment_type %in% c("rituximab", "obinutuzumab", "ofatumumab") ~ "CD20anti",
        TRUE ~ treatment_type))
  
   cohort_clean <- cohort_clean %>% 
     mutate(treatment_type = case_when(
       treatment_type %in% c("acalabrutinib", "ibrutinib", "other_chemo-other_imm-ibrutinib", "other_chemo-ibrutinib", "other_imm-ibrutinib", "rituximab-acalabrutinib", "rituximab-ibrutinib", "other_chemo-rituximab-ibrutinib", "obinutuzumab-ibrutinib", "obinutuzumab-acalabrutinib") ~ "btk-inhibitor",
       treatment_type %in% c("other_imm-ibrutinib-venetoclax", "other_chemo-ibrutinib-venetoclax", "rituximab-ibrutinib-venetoclax", "obinutuzumab-ibrutinib-venetoclax", "ofatumumab-ibrutinib-venetoclax") ~ "ibrutinib-venetoclax",
       treatment_type %in% c("other_chemo", "other_chemo-other_imm", "other_imm", "alemtuzumab", "idelalisib", "velcade-rituximab", "fludara-rituximab-idelalisib", "rituximab-idelalisib","cop-rituximab", "cop", "ccop", "chop", "chop-rituximab", "cyclofosfamid", "cyclofosfamid-rituximab", "minichop-rituximab", "oncovin-rituximab") ~ "other",
       treatment_type %in% c("other_chemo-obinutuzumab", "other_chemo-ofatumumab", "other_chemo-rituximab", "other_chemo-alemtuzumab") ~ "CD20anti",
       treatment_type %in% c("other_chemo-other_imm-venetoclax", "other_chemo-venetoclax", "other_chemo-rituximab-venetoclax", "other_chemo-obinutuzumab-venetoclax", "obinutuzumab-venetoclax", "venetoclax-obinutuzumab", "venetoclax-rituximab", "rituximab-venetoclax") ~ "venetoclax",
       treatment_type %in% c("fludarabin-rituximab", "fludara", "fludara-other_imm", "fcd", "fcd-rituximab", "fludara-alemtuzumab", "fludara-rituximab-ibrutinib", "fludara-ibrutinib") ~ "fludarabin",
       treatment_type %in% c("bendamustine-obinutuzumab", "bendamustine-ofatumumab", "bendamustine-other_imm", "bendamustine-rituximab", "bendamustine-alemtuzumab", "bendamustine-other_imm-ibrutinib", "bendamustine-ibrutinib", "bendamustine-rituximab-venetoclax", "bendamustine-venetoclax") ~ "bendamustine",
       treatment_type %in% c("chlorambucil-obinutuzumab", "chlorambucil-ofatumumab", "chlorambucil-other_imm", "chlorambucil-rituximab", "chlorambucil-alemtuzumab", "chlorambucil-rituximab-ibrutinib", "chlorambucil-rituximab-venetoclax") ~ "chlorambucil",
       treatment_type %in% c("abvd-rituximab") ~ "abvd",
       TRUE ~ treatment_type))

#Exclusion of patients with certain treatment types
   cohort_clean <- cohort_clean %>%
     filter(!treatment_type %in% c("abvd", "abvd-rituximab","cladribin", "cladribin-rituximab", "cns_ielsg", "cns_mvp-rituximab", "hdmtx"))

   sum(is.na(cohort_clean$treatment_type)) #NA's
   mean(is.na(cohort_clean$treatment_type)) * 100 #% NA's
   
   unique(cohort_clean$treatment_type) #8 different
   cohort_clean$treatment_type %>% table 
   
 
#### SAVE TABLES ####
  saveRDS(cohort_clean, "/ngc/projects2/dalyca_r/emmsii_r/sandbox/tables/1cohort_clean.rds")
  
  