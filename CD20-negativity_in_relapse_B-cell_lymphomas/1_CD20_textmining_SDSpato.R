# title: "CD20-data_Reading_pathology_data"
# output: html_document
# date: "2025-05-07"
# Author: Sanaz M. Gholy, Jojo Biel-Nielsen Ditez, Christian Brieghel 
#---


### Loads lymphoma patients and their pathology report to identify CD20 positivity
library(data.table)
library(tidyverse)
library(dplyr)
source('/ngc/projects2/dalyca_r/clean_r/load_dalycare_package.R')
load_dataset(c('t_dalycare_diagnoses', 'patient', 'RKKP_LYFO'))
setwd("/ngc/projects2/dalyca_r/sangho_r")

LYFO_clean = RKKP_LYFO %>% 
  clean_RKKP_LYFO() 
LYFO_clean$subtype %>% table

#Only patients with the three lymphoma subtypes who recived rituximab as part of firstline therapy
RR_CD20_cohort = LYFO_clean %>% 
  filter(subtype %in% c('DLBCL', 'FL', 'MCL'),
         immunotherapy_type_1st_line =='rituximab',
         !is.na(date_treatment_2nd_line))

RR_CD20_cohort %>% nrow_npatients()

#Reading pathology data
load_dataset('SDS_pato', RR_CD20_cohort$patientid)
SDS_pato_subset$patientid
SDS_pato_subset$k_rekvnr # patientid links to k_rekvnr in SDS_pato
SDS_pato_subset$c_snomedkode #so we can identify for only DLBCL, FL and MCL
SDS_pato_subset$d_rekvdato

#Loading pato-text, specifically for k_rekvnr in CD20_cohort
load_dataset('SDS_t_mikro_ny', value = SDS_pato_subset$k_rekvnr, column = 'k_rekvnr')  
SDS_t_mikro_ny_subset %>% head 
#This is not readable so we clean it
snomed_lymfom = c('M96003', 'M96803', 'M96903') #snomed codes for the three lymphoma subtypes
#pato_text_clean_cd20_intermediate1 = pato_text_clean %>% 
#mutate(text2 = clean_abbreviations(string = text))  # ETA: forever. Read by the datamanager in the Database 

##dataset with filtered sentence made by the datamanager cmf
load('/ngc/projects2/dalyca_r/sangho_r/data/intermediate1_cmf.RData')
pato_text_clean = SDS_t_mikro_ny_subset %>% clean_SDS_t_mikro() 

#Only rituximab at 1 line included
cohort_rituximab <- pato_text_clean %>%
  select(k_inst, k_rekvnr)

rm(pato_text_clean)
head(cohort_rituximab)

#fwrite(cohort_rituximab, "rekv_cohort_rituximab.tsv", sep = "\t")#We need to filter for only biopsi from the time of second line therapy

test=pato_text_clean_cd20_intermediate2 %>% # filters sentences (i.e. text between to punctuations) containing the pattern
  left_join(SDS_pato_subset %>% select(patientid, k_inst, k_rekvnr, date_pato = d_rekvdato, c_snomedkode) %>% distinct()) %>% 
  select(patientid, everything()) %>% 
  filter(!str_detect(c_snomedkode, "8$"),
         !str_detect(c_snomedkode, "X$")) %>% # exclude remission and obs pro
  left_join(RR_CD20_cohort %>% select(patientid, subtype, date_diagnosis, date_treatment_1st_line, date_treatment_2nd_line)) 

#joining test=pato_text_clean_cd20_intermediate
#Only 30 days after and 60 days before 2L treatment
pato_text_clean_cd20 = test %>% 
  mutate(time_2L_pato=diff_days(date_treatment_2nd_line,date_pato)) %>% 
  filter(time_2L_pato<=30,
         time_2L_pato>-60) %>% 
  select(-c_snomedkode) %>% 
  distinct()

nrow(pato_text_clean_cd20)

#write.table(pato_text_clean_cd20, file = "pato_text_clean.tsv", sep = "\t")
#write.table(pato_text_clean_cd20, file = "pato_text_clean.csv", sep = ";")


#Filter so text is where CD20 is mentiond
pato_text_clean_cd20_filtered<- pato_text_clean_cd20 %>% 
  filter_sentence(string = text, pattern = 'CD20')
head(pato_text_clean_cd20_filtered)

n_distinct(pato_text_clean_cd20_filtered$patientid)

#Unique patients in the filtered dataframe n=2076

#write.table(pato_text_clean_cd20_filtered, file = "pato_text_clean.csv", sep = ";",row.names = FALSE)

#for smooth text_mining i remove the 6 colomns i dont need to look at: date_diagnosis, date_treatment_1st_line, date_treatment_2st_line, time_2L_pato

cleaned_pato_text_data <- pato_text_clean_cd20_filtered %>% 
  select(patientid,k_inst, k_rekvnr,text, text2)
n_distinct(cleaned_pato_text_data$patientid)

#write.table(cleaned_pato_text_data, file = "pato_text_cd20_clean.tsv", sep = "\t")
#write.table(cleaned_pato_text_data, file = "pato_text_clean.csv", sep = ";",row.names = FALSE)

#write.csv(cleaned_pato_text_data[1:500,c("patientid","k_inst","k_rekvnr","text","text2")],"pato_text_clean")