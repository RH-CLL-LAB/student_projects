## This script finds biochemical AEs and adds them to the CLL and SLL cohort
# date: 22/4-2025
# Emma Siig and Christian Brieghel
# Rum 1cohort.R and 2Baseline_characteristics.R first or downloade table: 2cohort_riskfactors


source('/ngc/projects2/dalyca_r/clean_r/load_dalycare_package.R')


cohort_riskfactors <- readRDS("/ngc/projects2/dalyca_r/emmsii_r/sandbox/tables/2cohort_riskfactors.rds")
cohort_riskfactors_AE <- cohort_riskfactors

#set.seed(1234)
#cohort_test  = sample(RKKP_CLL$patientid, 200)
#load_dataset('laboratorymeasurements', cohort_test)


load_dataset('laboratorymeasurements', cohort_riskfactors_AE$patientid)
laboratorymeasurements_subset %>% n_patients()
laboratorymeasurements_subset %>% head2
load_npu_common()

#### Looking into patients in cohort without lab-measurements ####
#Finding patients with no lab measures
patients_no_labs <- cohort_riskfactors_AE %>%                          #Making cohort of patients with no lab values
  filter(!(patientid %in% laboratorymeasurements_subset$patientid))

patients_no_labs %>%                                              #Looking into where the patients come from (hospital id)
  count(hospital_id) %>%
  arrange(desc(n))

patients_no_labs %>%
  mutate(year_diagnosis = lubridate::year(date_diagnosis)) %>%    #Looking into date of diagnosis
  count(year_diagnosis) %>%
  arrange(year_diagnosis)


no_lab_summary <- cohort_riskfactors_AE %>%                     #Making no_lab_summary including date_diagnosis and hospital_id
  filter(!(patientid %in% laboratorymeasurements_subset$patientid)) %>%
  mutate(year_diagnosis = year(date_diagnosis)) %>%
  count(hospital_id, year_diagnosis)

ggplot(no_lab_summary, aes(x = factor(year_diagnosis), y = n, fill = factor(hospital_id))) +    #Makes plot of it
  geom_col() +
  labs(
    title = "Antal patienter uden labmålinger fordelt på hospital og diagnoseår",
    x = "Diagnoseår",
    y = "Antal patienter",
    fill = "Hospital ID"
  ) +
  theme_minimal()


#Treatment year
patients_no_labs_2 <- patients_no_labs%>%
  mutate(date_treatment_1st_line = lubridate::year(date_treatment_1st_line)) %>%    #Looking into date of treatment
  count(date_treatment_1st_line) %>%
  arrange(date_treatment_1st_line)

no_lab_summary_2 <- cohort_riskfactors_AE %>%                     #Making no_lab_summary including date_diagnosis and hospital_id
  filter(!(patientid %in% laboratorymeasurements_subset$patientid)) %>%
  mutate(year_treatment = year(date_treatment_1st_line)) %>%
  count(hospital_id, year_treatment)

ggplot(no_lab_summary_2, aes(x = factor(year_treatment), y = n, fill = factor(hospital_id))) +    #Makes plot of it
  geom_col() +
  labs(
    title = "Antal patienter uden labmålinger fordelt på hospital og behandlingsår",
    x = "Behandlingsår",
    y = "Antal patienter",
    fill = "Hospital ID"
  ) +
  theme_minimal()





#Sammenligner med dem, som har lab_data
lab_data_coverage <- cohort_riskfactors_AE %>%              # 1 = Patients WITH lab-data, 0 = patients WITHOUT lab-data
  mutate(has_labdata = if_else(patientid %in% laboratorymeasurements_subset$patientid, 1, 0))

lab_data_coverage <- lab_data_coverage %>%
  mutate(year_diagnosis = lubridate::year(date_diagnosis))

#overview_lab_data <- lab_data_coverage %>%
  group_by(hospital_id, year_diagnosis, has_labdata) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = has_labdata,         #Making with lab-results and without lab_result into the same row
                     values_from = n,
                     names_prefix = "lab_",
                     values_fill = 0) %>%
  mutate(summary = paste0("lab_yes: ", lab_1, " / lab_no: ", lab_0))   #Making a summary for each row


lab_data_coverage <- lab_data_coverage %>%                   #Cleaning
    mutate(hospital_id = case_when(
      hospital_id %in% c("NAE", "NAESTVED") ~ "NAESTVED",
      hospital_id %in% c("SHS", "SHS / AABENRAA") ~ "SHS/AABENRAA",
      TRUE ~ hospital_id
    ))
  
overview_lab_data <- lab_data_coverage %>%                  #Making variables for % with no lab_results
  group_by(hospital_id, year_diagnosis, has_labdata) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(
      names_from = has_labdata,
      values_from = n,
      names_prefix = "lab_",
      values_fill = 0
    ) %>%
    mutate(
      total = lab_0 + lab_1,
      pct_no_lab = round(100 * lab_0 / total, 1),
      summary = paste0("(", lab_0, "/", total, ") ", pct_no_lab, "%")
    )
  
overview_lab_data <- overview_lab_data %>% filter(lab_0 > 0)        #Filters out result with full coverage

#Plot with multiple bars for each hospital
ggplot(overview_lab_data, aes(x = factor(year_diagnosis), y = lab_0, fill = hospital_id)) +
  geom_col(position = "dodge") +
  geom_text(aes(label =  paste0(summary)),
            position = position_dodge(width = 0.9),
            angle = 90, size = 3) +
  labs(
    title = "Proportion of patients without laboratory measurements across hospitals and diagnosis years",
    x = "Year of diagnosis",
    y = "Number of patients without laboratory measurements",
    fill = "Hospital"
  ) +
  theme_minimal()

#Plot with one bar for each year - stacked (still divided into hospitals and years)
ggplot(overview_lab_data, aes(x = factor(year_diagnosis), y = lab_0, fill = hospital_id)) +
  geom_col() +
  geom_text(aes(label = summary), 
            position = position_stack(vjust = 0.3),     
            size = 3
  )+
  labs(
    title = "Proportion of patients without laboratory measurements across hospitals and diagnosis years",
    x = "Year of diagnosis",
    y = "Number of patients without laboratory measurements",
    fill = "Hospital"
  ) +
  theme_minimal()





#By treatment year
#Sammenligner med dem, som har lab_data
lab_data_coverage_2 <- cohort_riskfactors_AE %>%              # 1 = Patients WITH lab-data, 0 = patients WITHOUT lab-data
  mutate(has_labdata_2 = if_else(patientid %in% laboratorymeasurements_subset$patientid, 1, 0))

lab_data_coverage_2 <- lab_data_coverage_2 %>%
  mutate(year_treatment = lubridate::year(date_treatment_1st_line))

lab_data_coverage_2 <- lab_data_coverage_2 %>%                   #Cleaning
  mutate(hospital_id = case_when(
    hospital_id %in% c("NAE", "NAESTVED") ~ "NAESTVED",
    hospital_id %in% c("SHS", "SHS / AABENRAA") ~ "SHS/AABENRAA",
    TRUE ~ hospital_id
  ))

overview_lab_data_2 <- lab_data_coverage_2 %>%                  #Making variables for % with no lab_results
  group_by(hospital_id, year_treatment, has_labdata_2) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = has_labdata_2,
    values_from = n,
    names_prefix = "lab_",
    values_fill = 0
  ) %>%
  mutate(
    total = lab_0 + lab_1,
    pct_no_lab = round(100 * lab_0 / total, 1),
    summary = paste0("(", lab_0, "/", total, ") ", pct_no_lab, "%")
  )

overview_lab_data_2 <- overview_lab_data_2 %>% filter(lab_0 > 0)        #Filters out result with full coverage

#Plot with one bar for each year - stacked (still divided into hospitals and years)
ggplot(overview_lab_data_2, aes(x = factor(year_treatment), y = lab_0, fill = hospital_id)) +
  geom_col() +
  geom_text(aes(label = summary), 
            position = position_stack(vjust = 0.3),     
            size = 3
  )+
  labs(
    title = "Proportion of patients without laboratory measurements across hospitals and treatment years",
    x = "Year of 1st line treatment",
    y = "Number of patients without laboratory measurements",
    fill = "Hospital"
  ) +
  theme_minimal()



#### Anemia - ready ####
docstring(AE_anemia)
AE_anemia_data = laboratorymeasurements_subset %>%
  AE_anemia(days_grade_5 = 2) %>%                                   #Grade 5 is if a patient dies within 2 days
  left_join(cohort_riskfactors_AE %>% select(patientid, date_treatment_1st_line), by = "patientid") %>%  #inserts treatment date from cohort
  filter(
    !is.na(date_treatment_1st_line),                                 #Removes patients with no treatment start
    samplingdate >= date_treatment_1st_line,                         #Sampling date must be after/or start of treatment
    samplingdate <= date_treatment_1st_line + 180                    #Sampling date must be under/or 180 days after start of treatment
  ) %>%
  group_by(patientid) %>%
  arrange(patientid, desc(ae_anemia)) %>%                         #Sorts every patient by patientid and grade (highest first)
  slice(1) %>%                                                     #Chooses the worst anemia grade pr. patient
  ungroup() 
  #filter(ae_anemia != 0)                                           #Removes patients without anemia

AE_anemia_data %>% nrow_npatients()

utable(~ae_anemia, AE_anemia_data)
AE_anemia_data$ae_anemia %>% class
AE_anemia_data %>% head()

AE_anemia_data <- AE_anemia_data %>% #Renames samplingdate to samplingdate_anemia
  dplyr::rename(samplingdate_anemia = samplingdate)

cohort_riskfactors_AE1 <- cohort_riskfactors_AE

cohort_riskfactors_AE1 <- left_join( #Inserts samplingdate_anemia and ae_anemia in our cohort (cohort_riskfactors_AE1)
  cohort_riskfactors_AE1,
  AE_anemia_data %>% select(patientid, samplingdate_anemia, ae_anemia),
  by = "patientid"
)

#Binary data on anemia
str(AE_anemia_data$ae_anemia)

AE_anemia_data_bin = laboratorymeasurements_subset %>%
  AE_anemia(days_grade_5 = 2) %>%                                   #Grade 5 is if a patient dies within 2 days
  left_join(cohort_riskfactors_AE %>% select(patientid, date_treatment_1st_line), by = "patientid") %>%  #inserts treatment date from cohort
  filter(
    !is.na(date_treatment_1st_line),                                 #Removes patients with no treatment start
    samplingdate >= date_treatment_1st_line,                         #Sampling date must be after/or start of treatment
    samplingdate <= date_treatment_1st_line + 180                    #Sampling date must be under/or 180 days after start of treatment
  ) 

AE_anemia_data_bin <- AE_anemia_data_bin %>%
  mutate(anemia_binary = if_else(ae_anemia == "3-4", 1, 0))    #Making a binary variable. Grade 3-4 = anemia. Grade 0, 1 and 2 = no anemia

AE_anemia_data_bin <- AE_anemia_data_bin %>%
  filter(anemia_binary == 1) %>%                  #Only keeps measures where there is an event (anemia_binary = 1) 
  group_by(patientid) %>%                         #Groups by patientid
  arrange(samplingdate) %>%                       #Sorts dates
  slice(1) %>%                                    #Using only the first observation (first time anemia is measured)
  ungroup() %>%
  select(patientid, first_anemia_date = samplingdate, anemia_binary)  #Only keeps relevant variables and renames the date to first_anemia_date

AE_anemia_data_bin %>% nrow_npatients()

cohort_riskfactors_AE1 <- left_join( #Inserts first_anemia_date and anemia_binary in our cohort (cohort_riskfactors_AE1)
  cohort_riskfactors_AE1,
  AE_anemia_data_bin %>% select(patientid, first_anemia_date, anemia_binary),
  by = "patientid"
)

utable(~anemia_binary, cohort_riskfactors_AE1)

cohort_riskfactors_AE1 <- cohort_riskfactors_AE1 %>%  #Inserts 0 in anemia_binary for patients with a measurement/grade in ae_anemia and without an adverse event (anemia_binary = NA)
  mutate(anemia_binary = case_when(
    !is.na(ae_anemia) & is.na(anemia_binary) ~ 0,
    TRUE ~ anemia_binary                                   #If ae_anemia = NA, anemia_binary shall continue being NA
  ))

dead_without_anemia_in_followup <- cohort_riskfactors_AE1 %>%    #Finding patients who died within the follow-up period and didn't had anemia
  filter(
    anemia_binary == 0,
    status == 1,
    date_death_fu <= date_treatment_1st_line + 180
  )
nrow(dead_without_anemia_in_followup)                            #73 patients died


cohort_riskfactors_AE1 <- cohort_riskfactors_AE1 %>%         #Makes new variable: anemia_event_or_censor_date
  mutate(anemia_event_or_censor_date = case_when(
    anemia_binary == 1 ~ first_anemia_date,     #When a patient has anemia (=1), the first date of anemia is inserted
    anemia_binary == 0 & status == 1 & date_death_fu <= date_treatment_1st_line + 180 ~ date_death_fu, #When a patient is dead within the follow-up period without anemia, death date is inserted
    anemia_binary == 0 ~ date_treatment_1st_line + 180, #When a patient don't has anemia (=0), the last day of follow-up is inserted (treatment date + 180 days)
    TRUE ~ as.Date(NA) #NA if they both are NA
  ))



#### Neutropenia - ready ####
docstring(AE_neutrophil_count_decreased)
AE_neu_data$ae_neutrophil_count_decreased

AE_neu_data = laboratorymeasurements_subset %>% 
  AE_neutrophil_count_decreased() %>% 
  left_join(cohort_riskfactors_AE %>% select(patientid, date_treatment_1st_line), by = "patientid") %>%  #inserts treatment date from cohort
  filter(
    samplingdate >= date_treatment_1st_line,                         #Sampling date must be after/or start of treatment
    samplingdate <= date_treatment_1st_line + 180                    #Sampling date must be under/or 180 days after start of treatment
  ) %>%
  group_by(patientid) %>% 
  arrange(patientid, desc(ae_neutrophil_count_decreased)) %>%       #Sorts every patient by patientid and grade (highest first)
  slice(1) %>%                                                      #Chooses the worst white bloodcell count pr. patient
  ungroup()
# filter(ae_neutrophil_count_decreased != 0)                        #Removes patients without Neutropenia

AE_neu_data %>% nrow_npatients()

utable(~ae_neutrophil_count_decreased, AE_neu_data)

AE_neu_data <- AE_neu_data %>% #Renames samplingdate to samplingdate_neu
  dplyr::rename(samplingdate_neu = samplingdate)

cohort_riskfactors_AE1 <- left_join( #Inserts samplingdate_neu and ae_neutrophil_count_decreased in our cohort (cohort_riskfactors_AE1)
  cohort_riskfactors_AE1,
  AE_neu_data %>% select(patientid, samplingdate_neu, ae_neutrophil_count_decreased),
  by = "patientid"
)


#Binary data on neutropenia
str(AE_neu_data$ae_neutrophil_count_decreased)

AE_neu_data_bin = laboratorymeasurements_subset %>% 
  AE_neutrophil_count_decreased() %>% 
  left_join(cohort_riskfactors_AE %>% select(patientid, date_treatment_1st_line), by = "patientid") %>%  #inserts treatment date from cohort
  filter(
    samplingdate >= date_treatment_1st_line,                         #Sampling date must be after/or start of treatment
    samplingdate <= date_treatment_1st_line + 180                    #Sampling date must be under/or 180 days after start of treatment
  ) 

AE_neu_data_bin <- AE_neu_data_bin %>%                 #Making a binary variable. Grade 3 and 4 is 1. 0 = no neutropenia, 1 = neutropenia
  mutate(neutropenia_binary = case_when(
    as.character(ae_neutrophil_count_decreased) %in% c("3", "4") ~ 1,
    as.character(ae_neutrophil_count_decreased) %in% c("0", "1", "2") ~ 0,
    TRUE ~ NA_real_
  ))

AE_neu_data_bin <- AE_neu_data_bin %>%
  filter(neutropenia_binary == 1) %>%             #Only keeps measures where there is an event (neutropenia_binary = 1) 
  group_by(patientid) %>%                         #Groups by patientid
  arrange(samplingdate) %>%                       #Sorts dates
  slice(1) %>%                                    #Using only the first observation (first time neutropenia is measured)
  ungroup() %>%
  select(patientid, first_neutropenia_date = samplingdate, neutropenia_binary)  #Only keeps relevant variables and renames the date to first_neutropenia_date

AE_neu_data_bin %>% nrow_npatients()

cohort_riskfactors_AE1 <- left_join( #Inserts first_neutropenia_date and neutropenia_binary in our cohort (cohort_riskfactors_AE1)
  cohort_riskfactors_AE1,
  AE_neu_data_bin %>% select(patientid, first_neutropenia_date, neutropenia_binary),
  by = "patientid"
)

utable(~neutropenia_binary, cohort_riskfactors_AE1)

cohort_riskfactors_AE1 <- cohort_riskfactors_AE1 %>%  #Inserts 0 in neutropenia_binary for patients with a measurement/grade in ae_neutrophil_count_decreased and without an adverse event (neutropenia_binary = NA)
  mutate(neutropenia_binary = case_when(
    !is.na(ae_neutrophil_count_decreased) & is.na(neutropenia_binary) ~ 0,
    TRUE ~ neutropenia_binary                                   #If ae_neutrophil_count_decreased = NA, neutropenia_binary will still be NA
  ))

dead_without_neutropenia_in_followup <- cohort_riskfactors_AE1 %>%    #Finding patients who died within the follow-up period and didn't had neutropenia
  filter(
    neutropenia_binary == 0,
    status == 1,
    date_death_fu <= date_treatment_1st_line + 180
  )
nrow(dead_without_neutropenia_in_followup)                            #75 patients died


cohort_riskfactors_AE1 <- cohort_riskfactors_AE1 %>%         #Makes new variable: neutropenia_event_or_censor_date
  mutate(neutropenia_event_or_censor_date = case_when(
    neutropenia_binary == 1 ~ first_neutropenia_date,     #When a patient has neutropenia (=1), the first date of neutropenia is inserted
    neutropenia_binary == 0 & status == 1 & date_death_fu <= date_treatment_1st_line + 180 ~ date_death_fu, #When a patient is dead within the follow-up period without neutropenia, death date is insertet
    neutropenia_binary == 0 ~ date_treatment_1st_line + 180, #When a patient don't has neutropenia (=0), the last day of follow-up is inserted (treatment date + 180 days)
    TRUE ~ as.Date(NA) #NA if they both are NA
  ))


#### Lymphopenia - ready ####
docstring(AE_lymphocyte_count_decreased)

# AE_data_ALC$samplingdate

AE_lympho_data = laboratorymeasurements_subset %>% 
  AE_lymphocyte_count_decreased() %>% 
  left_join(cohort_riskfactors_AE %>% select(patientid, date_treatment_1st_line), by = "patientid") %>%  #inserts treatment date from cohort
  filter(
    samplingdate >= date_treatment_1st_line,                         #Sampling date must be after/or start of treatment
    samplingdate <= date_treatment_1st_line + 180                    #Sampling date must be under/or 180 days after start of treatment
  ) %>%
  group_by(patientid) %>% 
  arrange(patientid, desc(ae_lymphocyte_count_decreased)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  mutate(ae_ALC_G3_4 = recode_factor(ae_lymphocyte_count_decreased,
                                                       `0` = '0',
                                                       `1` = '1',
                                                       `2` = '2',
                                                       `3` = '3-4',
                                                       `4` = '3-4'))


AE_lympho_data %>% nrow_npatients()

utable(~ae_lymphocyte_count_decreased, AE_lympho_data)
utable(~ae_ALC_G3_4, AE_lympho_data)

AE_lympho_data <- AE_lympho_data %>% #Renames samplingdate to samplingdate_lympho
  dplyr::rename(samplingdate_lympho = samplingdate)

cohort_riskfactors_AE1 <- left_join( #Inserts samplingdate_lympho and ae_lymphocyte_count_decreased in our cohort (cohort_riskfactors_AE1)
  cohort_riskfactors_AE1,
  AE_lympho_data %>% select(patientid, samplingdate_lympho, ae_lymphocyte_count_decreased),
  by = "patientid"
)

#Binary data on lymphopenia
AE_lympho_data_bin = laboratorymeasurements_subset %>% 
  AE_lymphocyte_count_decreased() %>% 
  left_join(cohort_riskfactors_AE %>% select(patientid, date_treatment_1st_line), by = "patientid") %>%  #inserts treatment date from cohort
  filter(
    samplingdate >= date_treatment_1st_line,                         #Sampling date must be after/or start of treatment
    samplingdate <= date_treatment_1st_line + 180                    #Sampling date must be under/or 180 days after start of treatment
  )

AE_lympho_data_bin <- AE_lympho_data_bin %>%                 #Making a binary variable. Grade 3 and 4 is 1. 0 = no lymphopenia, 1 = lymphopenia
  mutate(lymphopenia_binary = case_when(
    as.character(ae_lymphocyte_count_decreased) %in% c("3", "4") ~ 1,
    as.character(ae_lymphocyte_count_decreased) %in% c("0", "1", "2") ~ 0,
    TRUE ~ NA_real_
  ))

AE_lympho_data_bin <- AE_lympho_data_bin %>%
  filter(lymphopenia_binary == 1) %>%             #Only keeps measures where there is an event (lymphopenia_binary = 1) 
  group_by(patientid) %>%                         #Groups by patientid
  arrange(samplingdate) %>%                       #Sorts dates
  slice(1) %>%                                    #Using only the first observation (first time lymphopenia is measured)
  ungroup() %>%
  select(patientid, first_lymphopenia_date = samplingdate, lymphopenia_binary)  #Only keeps relevant variables and renames the date to first_lymphopenia_date

AE_lympho_data_bin %>% nrow_npatients()

cohort_riskfactors_AE1 <- left_join( #Inserts first_lymphopenia_date and lymphopenia_binary in our cohort (cohort_riskfactors_AE1)
  cohort_riskfactors_AE1,
  AE_lympho_data_bin %>% select(patientid, first_lymphopenia_date, lymphopenia_binary),
  by = "patientid"
)

utable(~lymphopenia_binary, cohort_riskfactors_AE1)

cohort_riskfactors_AE1 <- cohort_riskfactors_AE1 %>%  #Inserts 0 in lymphopenia_binary for patients with a measurement/grade in ae_lymphocyte_count_decreased and without an adverse event (lymphopenia_binary = NA)
  mutate(lymphopenia_binary = case_when(
    !is.na(ae_lymphocyte_count_decreased) & is.na(lymphopenia_binary) ~ 0,
    TRUE ~ lymphopenia_binary                                   #If ae_lymphocyte_count_decreased = NA, lymphopenia_binary will still be NA
  ))

dead_without_lymphopenia_in_followup <- cohort_riskfactors_AE1 %>%    #Finding patients who died within the follow-up period and didn't had lymphopenia
  filter(
    lymphopenia_binary == 0,
    status == 1,
    date_death_fu <= date_treatment_1st_line + 180
  )
nrow(dead_without_lymphopenia_in_followup)                            #76 patients died

cohort_riskfactors_AE1 <- cohort_riskfactors_AE1 %>%         #Makes new variable: lymphopenia_event_or_censor_date
  mutate(lymphopenia_event_or_censor_date = case_when(
    lymphopenia_binary == 1 ~ first_lymphopenia_date,     #When a patient has lymphopenia (=1), the first date of lymphopenia is inserted
    lymphopenia_binary == 0 & status == 1 & date_death_fu <= date_treatment_1st_line + 180 ~ date_death_fu, #When a patient is dead within the follow-up period without lymphopenia, death date is inserted
    lymphopenia_binary == 0 ~ date_treatment_1st_line + 180, #When a patient don't has lymphopenia (=0), the last day of follow-up is inserted (treatment date + 180 days)
    TRUE ~ as.Date(NA) #NA if they both are NA
  ))

#### Creatinine - ready ####
AE_creatinine_increased()
AE_crea_data = laboratorymeasurements_subset %>% 
  AE_creatinine_increased() %>% 
  left_join(cohort_riskfactors_AE %>% select(patientid, date_treatment_1st_line), by = "patientid") %>%  #inserts treatment date from cohort
  filter(
    samplingdate >= date_treatment_1st_line,                         #Sampling date must be after/or start of treatment
    samplingdate <= date_treatment_1st_line + 180                    #Sampling date must be under/or 180 days after start of treatment
  ) %>%
  group_by(patientid) %>% 
  arrange(patientid, desc(ae_creatinine_increased)) %>% 
  slice(1) %>% 
  ungroup() 

AE_crea_data %>% nrow_npatients()

utable(~ae_creatinine_increased, AE_crea_data)

AE_crea_data <- AE_crea_data %>% #renames samplingdate to samplingdate_crea
  dplyr::rename(samplingdate_crea = samplingdate)

cohort_riskfactors_AE1 <- left_join( #Inserts samplingdate_crea and ae_creatinine_increased in our cohort (cohort_riskfactors_AE1)
  cohort_riskfactors_AE1,
  AE_crea_data %>% select(patientid, samplingdate_crea, ae_creatinine_increased),
  by = "patientid"
)

#Binary data on creatinine
AE_crea_data_bin = laboratorymeasurements_subset %>% 
  AE_creatinine_increased() %>% 
  left_join(cohort_riskfactors_AE %>% select(patientid, date_treatment_1st_line), by = "patientid") %>%  #inserts treatment date from cohort
  filter(
    samplingdate >= date_treatment_1st_line,                         #Sampling date must be after/or start of treatment
    samplingdate <= date_treatment_1st_line + 180                    #Sampling date must be under/or 180 days after start of treatment
  ) 

AE_crea_data_bin <- AE_crea_data_bin %>%                 #Making a binary variable. Grade 2, 3 and 4 is 1. 0 = no creatinine increase (no AKI), 1 = creatinine increase (AKI)
  mutate(creatinine_binary = case_when(
    as.character(ae_creatinine_increased) %in% c("2", "3", "4") ~ 1,
    as.character(ae_creatinine_increased) %in% c("0", "1") ~ 0,
    TRUE ~ NA_real_
  ))

AE_crea_data_bin <- AE_crea_data_bin %>%
  filter(creatinine_binary == 1) %>%              #Only keeps measures where there is an event (creatinine_binary = 1) 
  group_by(patientid) %>%                         #Groups by patientid
  arrange(samplingdate) %>%                       #Sorts dates
  slice(1) %>%                                    #Using only the first observation (first time creatinine increase is measured)
  ungroup() %>%
  select(patientid, first_creatinine_date = samplingdate, creatinine_binary)  #Only keeps relevant variables and renames the date to first_creatinine_date

AE_crea_data_bin %>% nrow_npatients()

cohort_riskfactors_AE1 <- left_join( #Inserts first_creatinine_date and creatinine_binary in our cohort (cohort_riskfactors_AE1)
  cohort_riskfactors_AE1,
  AE_crea_data_bin %>% select(patientid, first_creatinine_date, creatinine_binary),
  by = "patientid"
)

utable(~creatinine_binary, cohort_riskfactors_AE1)

cohort_riskfactors_AE1 <- cohort_riskfactors_AE1 %>%  #Inserts 0 in creatinine_binary for patients with a measurement/grade in ae_creatinine_increased and without an adverse event (creatinine_binary = NA)
  mutate(creatinine_binary = case_when(
    !is.na(ae_creatinine_increased) & is.na(creatinine_binary) ~ 0,
    TRUE ~ creatinine_binary                                   #If ae_creatinine_increased = NA, creatinine_binary will still be NA
  ))

dead_without_creatinine_in_followup <- cohort_riskfactors_AE1 %>%    #Finding patients who died within the follow-up period and didn't had creatinine increase (AKI)
  filter(
    creatinine_binary == 0,
    status == 1,
    date_death_fu <= date_treatment_1st_line + 180
  )
nrow(dead_without_creatinine_in_followup)                            #82 patients died

cohort_riskfactors_AE1 <- cohort_riskfactors_AE1 %>%         #Makes new variable: creatinine_event_or_censor_date
  mutate(creatinine_event_or_censor_date = case_when(
    creatinine_binary == 1 ~ first_creatinine_date,     #When a patient has creatinine (=1), the first date of creatinine increase (AKI) is inserted
    creatinine_binary == 0 & status == 1 & date_death_fu <= date_treatment_1st_line + 180 ~ date_death_fu, #When a patient is dead within the follow-up period without creatinine increase, death date is inserted
    creatinine_binary == 0 ~ date_treatment_1st_line + 180, #When a patient don't has creatinine increase (=0), the last day of follow-up is inserted (treatment date + 180 days)
    TRUE ~ as.Date(NA) #NA if they both are NA
  ))

#### Thrombocytopenia (platelets) - ready ####
AE_platelets_decreased() 
docstring(AE_platelets_decreased)

AE_thrombo_data = laboratorymeasurements_subset %>% 
  AE_platelets_decreased() %>% 
  left_join(cohort_riskfactors_AE %>% select(patientid, date_treatment_1st_line), by = "patientid") %>%  #inserts treatment date from cohort
  filter(
    samplingdate >= date_treatment_1st_line,                         #Sampling date must be after/or start of treatment
    samplingdate <= date_treatment_1st_line + 180                    #Sampling date must be under/or 180 days after start of treatment
  ) %>%
  group_by(patientid) %>% 
  arrange(patientid, desc(ae_platelets_decreased)) %>% 
  slice(1) %>% 
  ungroup() 

AE_thrombo_data %>% nrow_npatients()

utable(~ae_platelets_decreased, AE_thrombo_data)

AE_thrombo_data <- AE_thrombo_data %>% #renames samplingdate to samplingdate_thrombo
  dplyr::rename(samplingdate_thrombo = samplingdate)

cohort_riskfactors_AE1 <- left_join( #Inserts samplingdate_thrombo and ae_platelets_decreased in our cohort (cohort_riskfactors_AE1)
  cohort_riskfactors_AE1,
  AE_thrombo_data %>% select(patientid, samplingdate_thrombo, ae_platelets_decreased),
  by = "patientid"
)

#Binary data on thrombocytopenia
AE_thrombo_data_bin = laboratorymeasurements_subset %>% 
  AE_platelets_decreased() %>% 
  left_join(cohort_riskfactors_AE %>% select(patientid, date_treatment_1st_line), by = "patientid") %>%  #inserts treatment date from cohort
  filter(
    samplingdate >= date_treatment_1st_line,                         #Sampling date must be after/or start of treatment
    samplingdate <= date_treatment_1st_line + 180                    #Sampling date must be under/or 180 days after start of treatment
  )

AE_thrombo_data_bin <- AE_thrombo_data_bin %>%                 #Making a binary variable. Grade 3 and 4 is 1. 0 = no thrombocytopenia, 1 = thrombocytopenia
  mutate(thrombocytopenia_binary = case_when(
    as.character(ae_platelets_decreased) %in% c("3", "4") ~ 1,
    as.character(ae_platelets_decreased) %in% c("0", "1", "2") ~ 0,
    TRUE ~ NA_real_
  ))

AE_thrombo_data_bin <- AE_thrombo_data_bin %>%
  filter(thrombocytopenia_binary == 1) %>%        #Only keeps measures where there is an event (thrombocytopenia_binary = 1) 
  group_by(patientid) %>%                         #Groups by patientid
  arrange(samplingdate) %>%                       #Sorts dates
  slice(1) %>%                                    #Using only the first observation (first time thrombocytopenia (platelets decreased) is measured)
  ungroup() %>%
  select(patientid, first_thrombocytopenia_date = samplingdate, thrombocytopenia_binary)  #Only keeps relevant variables and renames the date to first_thrombocytopenia_date

AE_thrombo_data_bin %>% nrow_npatients()

cohort_riskfactors_AE1 <- left_join( #Inserts first_thrombocytopenia_date and thrombocytopenia_binary in our cohort (cohort_riskfactors_AE1)
  cohort_riskfactors_AE1,
  AE_thrombo_data_bin %>% select(patientid, first_thrombocytopenia_date, thrombocytopenia_binary),
  by = "patientid"
)

utable(~thrombocytopenia_binary, cohort_riskfactors_AE1)

cohort_riskfactors_AE1 <- cohort_riskfactors_AE1 %>%  #Inserts 0 in thrombocytopenia_binary for patients with a measurement/grade in ae_platelets_decreased and without an adverse event (thrombocytopenia_binary = NA)
  mutate(thrombocytopenia_binary = case_when(
    !is.na(ae_platelets_decreased) & is.na(thrombocytopenia_binary) ~ 0,
    TRUE ~ thrombocytopenia_binary                                   #If ae_platelets_decreased = NA, thrombocytopenia_binary will still be NA
  ))

dead_without_thrombocytopenia_in_followup <- cohort_riskfactors_AE1 %>%    #Finding patients who died within the follow-up period and didn't had thrombocytopenia
  filter(
    thrombocytopenia_binary == 0,
    status == 1,
    date_death_fu <= date_treatment_1st_line + 180
  )
nrow(dead_without_thrombocytopenia_in_followup)                            #89 patients died

cohort_riskfactors_AE1 <- cohort_riskfactors_AE1 %>%         #Makes new variable: thrombocytopenia_event_or_censor_date
  mutate(thrombocytopenia_event_or_censor_date = case_when(
    thrombocytopenia_binary == 1 ~ first_thrombocytopenia_date,     #When a patient has thrombocytopenia (=1), the first date of thrombocytopenia is inserted
    thrombocytopenia_binary == 0 & status == 1 & date_death_fu <= date_treatment_1st_line + 180 ~ date_death_fu, #When a patient is dead within the follow-up period without thrombocytopenia, death date is inserted
    thrombocytopenia_binary == 0 ~ date_treatment_1st_line + 180, #When a patient don't has thrombocytopenia (=0), the last day of follow-up is inserted (treatment date + 180 days)
    TRUE ~ as.Date(NA) #NA if they both are NA
  ))

#### Leukopenia (WBC) - ready ####
AE_WBC_decreased()
AE_leuko_data$WBC_decreased
docstring(AE_WBC_decreased)

AE_leuko_data = laboratorymeasurements_subset %>%
  AE_WBC_decreased() %>%  
  left_join(cohort_riskfactors_AE %>% select(patientid, date_treatment_1st_line), by = "patientid") %>%  #inserts treatment date from cohort
  filter(
    samplingdate >= date_treatment_1st_line,                         #Sampling date must be after/or start of treatment
    samplingdate <= date_treatment_1st_line + 180                    #Sampling date must be under/or 180 days after start of treatment
  ) %>%
  group_by(patientid) %>% 
  arrange(patientid, desc(ae_WBC_decreased)) %>% 
  slice(1) %>% 
  ungroup() 

AE_leuko_data %>% nrow_npatients()

utable(~ae_WBC_decreased, AE_leuko_data)

AE_leuko_data <- AE_leuko_data %>% #renames samplingdate to samplingdate_leuko
  dplyr::rename(samplingdate_leuko = samplingdate)

cohort_riskfactors_AE1 <- left_join( #Inserts samplingdate_leuko and ae_WBC_decreased in our cohort (cohort_riskfactors_AE1)
  cohort_riskfactors_AE1,
  AE_leuko_data %>% select(patientid, samplingdate_leuko, ae_WBC_decreased),
  by = "patientid"
)

#Binary data on leukopenia
AE_leuko_data_bin = laboratorymeasurements_subset %>%
  AE_WBC_decreased() %>%  
  left_join(cohort_riskfactors_AE %>% select(patientid, date_treatment_1st_line), by = "patientid") %>%  #inserts treatment date from cohort
  filter(
    samplingdate >= date_treatment_1st_line,                         #Sampling date must be after/or start of treatment
    samplingdate <= date_treatment_1st_line + 180                    #Sampling date must be under/or 180 days after start of treatment
  )

AE_leuko_data_bin <- AE_leuko_data_bin %>%                 #Making a binary variable. Grade 3 and 4 is 1. 0 = no leukopenia, 1 = leukopenia
  mutate(leukopenia_binary = case_when(
    as.character(ae_WBC_decreased) %in% c("3", "4") ~ 1,
    as.character(ae_WBC_decreased) %in% c("0", "1", "2") ~ 0,
    TRUE ~ NA_real_
  ))

AE_leuko_data_bin <- AE_leuko_data_bin %>%
  filter(leukopenia_binary == 1) %>%              #Only keeps measures where there is an event (leukopenia_binary = 1) 
  group_by(patientid) %>%                         #Groups by patientid
  arrange(samplingdate) %>%                       #Sorts dates
  slice(1) %>%                                    #Using only the first observation (first time leukopenia (WBC decreased) is measured)
  ungroup() %>%
  select(patientid, first_leukopenia_date = samplingdate, leukopenia_binary)  #Only keeps relevant variables and renames the date to first_leukopenia_date

AE_leuko_data_bin %>% nrow_npatients()

cohort_riskfactors_AE1 <- left_join( #Inserts first_leukopenia_date and leukopenia_binary in our cohort (cohort_riskfactors_AE1)
  cohort_riskfactors_AE1,
  AE_leuko_data_bin %>% select(patientid, first_leukopenia_date, leukopenia_binary),
  by = "patientid"
)

utable(~leukopenia_binary, cohort_riskfactors_AE1)

cohort_riskfactors_AE1 <- cohort_riskfactors_AE1 %>%  #Inserts 0 in leukopenia_binary for patients with a measurement/grade in ae_WBC_decreased and without an adverse event (leukopenia_binary = NA)
  mutate(leukopenia_binary = case_when(
    !is.na(ae_WBC_decreased) & is.na(leukopenia_binary) ~ 0,
    TRUE ~ leukopenia_binary                                   #If ae_WBC_decreased = NA, leukopenia_binary will still be NA
  ))

#cohort_riskfactors_AE1 <- cohort_riskfactors_AE1 %>%                    #Makes new variable: leukopenia_event_or_censor_date
#  mutate(
#    leukopenia_event_or_censor_date = case_when(
#      leukopenia_binary == 1 ~ first_leukopenia_date,             #When a patient has leukopenia (=1), the first date of leukopenia is inserted
#      leukopenia_binary == 0 ~ date_treatment_1st_line + 180,     #When a patient don't has leukopenia (=0), the last day of follow-up is inserted (treatment date + 180 days)  
#      TRUE ~ NA_Date_                                             #NA if they both are NA
#    ))


dead_without_leukopenia_in_followup <- cohort_riskfactors_AE1 %>%    #Finding patients who died within the follow-up period and didn't had leukopenia
  filter(
    leukopenia_binary == 0,
    status == 1,
    date_death_fu <= date_treatment_1st_line + 180
  )
nrow(dead_without_leukopenia_in_followup)                            #99 patients died

cohort_riskfactors_AE1 <- cohort_riskfactors_AE1 %>%         #Makes new variable: leukopenia_event_or_censor_date
  mutate(leukopenia_event_or_censor_date = case_when(
    leukopenia_binary == 1 ~ first_leukopenia_date,     #When a patient has leukopenia (=1), the first date of leukopenia is inserted
    leukopenia_binary == 0 & status == 1 & date_death_fu <= date_treatment_1st_line + 180 ~ date_death_fu, #When a patient is dead within the follow-up period without leukopenia, death date is inserted
    leukopenia_binary == 0 ~ date_treatment_1st_line + 180, #When a patient don't has leukopenia (=0), the last day of follow-up is inserted (treatment date + 180 days)
    TRUE ~ as.Date(NA) #NA if they both are NA
  ))

#### AE_hospitalization - not using ####
#load_dataset('SP_ADT_haendelser', cohort_test)
#SP_ADT = SP_ADT_haendelser_subset %>% 
#  AE_hospitalization()

#### AE_infections - not using ####
#load_dataset ('patient')
#load_dataset('SP_AdministreretMedicin', cohort_test)
#SP_AB = SP_AdministreretMedicin_subset %>% ATC_AB() 
#SP_infections = SP_AB %>% AE_infection()


#### AKI - not using ####
# AE_AKI() # watch out
#load_dataset('SDS_lab_forsker', cohort_test) #loads creatinine # with laboratorymeasurements?
#SDS_lab_forsker_subset %>% head2
#CREATININE_clean = SDS_lab_forsker_subset %>% 
#  filter(analysiscode %in% NPU.KREA) %>% 
#  clean_lab_values() #cleans creatinine values
#AKI = CREATININE_clean %>%  AE_AKI(value = value2)



#### Summary ####
summary(cohort_riskfactors_AE1)


#Extra code for binary variables (NOT USING)
#cohort_riskfactors_AE1 <- cohort_riskfactors_AE1 %>%                    #Makes new variable: thrombocytopenia_event_or_censor_date
#  mutate(
#    thrombocytopenia_event_or_censor_date = case_when(
#      thrombocytopenia_binary == 1 ~ first_thrombocytopenia_date,       #When a patient has thrombocytopenia (=1), the first date of thrombocytopenia is inserted
#      thrombocytopenia_binary == 0 ~ date_treatment_1st_line + 180,     #When a patient don't has thrombocytopenia (=0), the last day of follow-up is inserted (treatment date + 180 days)  
#      TRUE ~ NA_Date_                                                   #NA if they both are NA
#    ))


#### SAVE TABLES ####
saveRDS(cohort_riskfactors_AE1, "/ngc/projects2/dalyca_r/emmsii_r/sandbox/tables/3cohort_riskfactors_AE1.rds")

