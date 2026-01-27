# title: "CD20-cleaning cohort"
# output: html_document
# date: "2025-12-17"
# Author: Sanaz M. Gholy
#---

library(data.table)
library(tidyverse)
library(survival)
library(dplyr)
library(survminer)

source('/ngc/projects2/dalyca_r/clean_r/load_dalycare_package.R') 

load_dataset(c('t_dalycare_diagnoses', 'patient', 'RKKP_LYFO','clean_RKKP_LYFO', 'SDS_pato'))
RKKP_LYFO <- readRDS('/ngc/projects2/dalyca_r/andkat_r/RT/Data/RKKP_LYFO_oct25.rds')
LYFO_clean = RKKP_LYFO %>% 
  clean_RKKP_LYFO()
names(LYFO_clean)

setwd("/ngc/projects2/dalyca_r/sangho_r/data")
cohort <- fread("studycohort_inkl_cd20_JBD_20251127.txt") %>% #Jojo has fixed the problem with missing date_pato in this script
  select(patientid, k_rekvnr, date_pato, CD20)
#CD20 status here is the final 

#All pato from the patients in the cohort
pato <- fread("/ngc/projects2/dalyca_r/sangho_r/data/brug_disse_patobeskrivelser_full_20250915.txt") 

#All the 2664 pathology report that has been manually reviewed
pato1 <- pato %>%
  group_by(patientid) %>%
  slice(1) %>% #one row each patient
  ungroup() %>% 
  select(patientid)

cohort1 <- cohort %>%
  left_join(pato1, select("patientid", "k_rekvnr", "date_pato"),
            by= "patientid") 

any(is.na(cohort1))# k_rekvnr, date_pato and CD20 status on all patients in the cohort

table(cohort1$CD20)

#I want colomns for CD20 as three categories (pos, red, neg), and one which is the origian (E (expressed), L (low), N (negative))
cohort2 <- cohort1 %>% 
 dplyr::rename(CD20_raw = CD20) %>% 
  mutate(CD20= recode(CD20_raw, 
                      "E" = "Positive",
                      "N" = "Negative",
                      "L" = "Negative"),
         CD20_3cat = recode(CD20_raw,
                            "E" = "Positive",
                            "N" = "Negative",
                            "L" = "Reduced"))

#Sanity check
n_distinct(cohort2$patientid)
cohort2 %>% count(CD20) #this is the cohort before excluding PBL, and the extra 6 PCNSL patients

#ASCT patient making trouble, reads from another file
setwd("/ngc/projects2/dalyca_r/sangho_r/data/CD20_final_data_SG")
ASCT_id <- fread("ASCT_pt.txt") %>% 
  select(c('patientid', 'k_rekvnr', 'date_pato', 'CD20', 'CD20_raw', 'CD20_3cat'))

cohort2<- cohort2 %>% mutate(patientid= as.character(patientid))
cohort2<- cohort2 %>% mutate(date_pato= as.Date(date_pato))

ASCT_id<- ASCT_id %>% mutate(patientid= as.character(patientid))
ASCT_id<- ASCT_id %>% mutate(k_rekvnr= as.character(k_rekvnr))
ASCT_id<- ASCT_id %>% mutate(date_pato= as.Date(date_pato))
cohort3 <- cohort2%>% 
  bind_rows(ASCT_id)
#Excluding  PBL 
getwd()
setwd("/ngc/projects2/dalyca_r/sangho_r/data/")
PBL_id <- fread("PBL_id.txt") #identified and saved in another file

cohort4 <- cohort3 %>% 
  filter(!patientid %in% PBL_id$patientid)

table(cohort4$CD20)
table(cohort4$CD20_raw)
n_distinct(cohort4)

#Primary CNS lymphoma
#This is from a script by PB and CB with PCNSL definition. This includes aditional 6 patients being excluded due to PCNSL 
LYFO_clean= RKKP_LYFO %>% 
  clean_RKKP_LYFO() %>% 
  mutate(CNS_inclusion = factor(ifelse(CNS_diagnosis == 'Yes'
                                       | CNS_involvement_diagnosis == 'Yes' 
                                       | leptomeninges_diagnosis == 'Yes', 
                                       'Yes', 'No'), c('Yes', 'No'))) %>% 
  mutate(nodal_disease_diagnosis = ifelse(cervical_diagnosis == 'No'
                                          & supraclavicular_diagnosis == 'No'
                                          & infraclavicular_diagnosis == 'No' 
                                          & axillary_diagnosis == 'No'
                                          & mediastinum_diagnosis == 'No'
                                          & hilar_diagnosis == 'No' 
                                          & retroperitoneum_diagnosis == 'No'
                                          & mesenteric_diagnosis == 'No'
                                          & pelvic_diagnosis == 'No'
                                          & inguinal_diagnosis == 'No'
                                          & bone_marrow_diagnosis == 'No'
                                          & spleen_diagnosis == 'No', 
                                          'No', 'Yes')) %>% 
  mutate(PCNSL = ifelse(CNS_inclusion == 'Yes' & nodal_disease_diagnosis == 'No' & n_extranodal_regions_diagnosis ==1, 
                        'Yes', 'No'))
LYFO_clean <- LYFO_clean %>% 
  mutate(patientid = as.character(patientid),
    PCNSL = as.character(PCNSL))

cohort5 <- cohort4 %>%
  left_join(LYFO_clean %>% 
              select(patientid,PCNSL), by='patientid') %>% 
  mutate(PCNSL = as.character(PCNSL))

cohort5$PCNSL %>%  table

#filtering an additional pcnsl patients out since this patient receive CNS directed therapy at first line- revised with CB
setwd("/ngc/projects2/dalyca_r/sangho_r/data/CD20_final_data_SG/")
PCNSL_id<- fread("PCNSL.txt")
cohort6 <- cohort5 %>% 
  filter(PCNSL != "Yes"|is.na(PCNSL))
cohort7 <- cohort6 %>% 
  filter(!patientid %in% PCNSL_id)

n_distinct(cohort7)

#fwrite(cohort7, "CD20_cohort_20251203.txt", sep = "\t") #This is the final cohort



