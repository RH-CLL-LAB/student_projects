# title: "CD20-data"
# output: html_document
# date: "2025-11-19"
# Author: Christian Brieghel 
#---

library(data.table)

source('/ngc/projects2/dalyca_r/clean_r/load_dalycare_package.R') 
setwd('/ngc/projects2/dalyca_r/sangho_r/data/CD20_final_data_SG')
study_population <- fread("CD20_df.txt") %>% 
  select(patientid, date_diagnosis, date_treatment_2nd_line)

#### load data ####
load_dataset(c('SDS_ekokur', 'SDS_epikur'), study_population$patientid)

Medicine_all = bind_rows(SDS_epikur_subset %>% transmute(patientid,
                                                         date = eksd,
                                                         atc,
                                                         source = 'LSR'),
                         SDS_ekokur_subset %>% transmute(patientid,
                                                         date = eksd,
                                                         atc,
                                                         source = 'LSR')) %>% 
  mutate(atc = gsub(':|_|\u008f\u008f|ÅÅ|XXX', '', atc),
         atc = gsub('Å', '', atc)) %>% 
  filter(! atc %in% c("N/A", "NUL", '', 'X')) # Time-lapse mins


date_max_atc = Medicine_all$date %>% max
study_population$date_treatment_2nd_line
study_population$date_diagnosis
study_population_polyrx = study_population %>%  
  mutate(date_treatment_2nd_line < date_max_atc) # patients with prescription data
study_population_polyrx %>% nrow_npatients()
cohort_1 = study_population_polyrx$patientid 

Medicine_all %>% nrow_npatients() #Mister en pt her?
#diagnoses_all_subset %>% nrow_npatients(), der er en error men bliver så heller ikke brugt længere nede??

# DX.FIRST = study_population_polyrx

#### N ATC ####
#### PolyRX ####
Medicine_all.1 = Medicine_all %>%
  filter(patientid %in% study_population$patientid) %>% 
  distinct() %>%
  group_by(patientid) %>% 
  mutate(N = n()) %>% 
  slice(1) %>% 
  ungroup() %>% 
  right_join(study_population %>% select(patientid), by = 'patientid') %>% 
  mutate(N = ifelse(is.na(N), 0, N))

## ANY ATC level
ggplot(Medicine_all) +
  geom_histogram(aes(date)) # atc coverage
study_population$date_treatment_2nd_line
DX.FIRST.ATC = study_population_polyrx %>% 
  select(patientid, date_treatment_2nd_line) %>% 
  filter(date_treatment_2nd_line >= ymd('2002-01-01'), # t_dalycare_diagnosis first date
         date_treatment_2nd_line <= as.Date(date_max_atc)) %>%  #LSR max date
  left_join(Medicine_all %>% distinct %>% select(patientid, date_atc = date, atc), 'patientid') %>% 
  distinct()


DX.FIRST.ATC %>% nrow # REPORT in Suppl Info #n.Rx
DX.FIRST.ATC$atc %>% n_distinct() #n.ATC among Rx
DX.FIRST.ATC %>% n_patients()  == n_distinct(cohort_1)
DX.FIRST.ATC %>% 
  left_join(study_population %>% select(patientid, date_treatment_2nd_line), 'patientid') %>% 
  mutate(SAME = date_treatment_2nd_line.x == date_treatment_2nd_line.y) %>% 
  pull(SAME) %>% 
  table # must be TRUE only

# LSR before first DX
DX.FIRST.ATC2 = DX.FIRST.ATC %>% 
  left_join(study_population %>% select(patientid, date_treatment_2nd_line), 'patientid') %>% 
  mutate(Time = diff_days(date_treatment_2nd_line.y, date_atc)) %>%
  filter(Time <= 0,
         Time > -365.25) %>% #0 + 365.25 
  left_join(study_population_polyrx %>% select(patientid, date = date_treatment_2nd_line) %>% 
              filter(date <= as.Date(date_max_atc)), 'patientid') 

exists("DX.FIRST.ATC2")


library(dplyr)
df <- study_population %>% 
  inner_join(DX.FIRST.ATC2, by= "patientid") %>% 
  filter(date<= as.Date(date_diagnosis)) 

  DX.FIRST.POLY <- df %>% 
  group_by(patientid, atc) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(patientid) %>% 
  summarise(n.ATC = n(), .groups="drop") %>% 
    right_join(study_population_polyrx %>% 
                select(patientid, date_diagnosis) %>% 
                 filter(date_diagnosis  <= as.Date(date_max_atc)),
               by= "patientid") %>% 
  mutate(n.ATC = replace_na(n.ATC, 0L), 
         ATC_group = ifelse(n.ATC < 5, '<5', '>5'))

DX.FIRST.POLY %>%  head() # This is the one

getwd()
#saveRDS(DX.FIRST.POLY, file = "DX_FIRST_POLY_SG_20251119.rds")
