# title: "CD20-data_Baseline clinical data"
# output: html_document
# date: "2025-12-17"
# Author: Sanaz M. Gholy
#---

library(data.table)
library(tidyverse)

setwd("/ngc/projects2/dalyca_r/sangho_r/data/CD20_final_data_SG")
study_population <- fread("CD20_cohort_20251203.txt")

n_distinct(study_population) #1396
sum(is.na(study_population)) #0 So no missing values for k_rekvnr, date_pato, CD20 status 

source('/ngc/projects2/dalyca_r/clean_r/load_dalycare_package.R') 
load_dataset(c("t_dalycare_diagnoses","patient", "RKKP_LYFO", "IPI"))

LYFO_clean = RKKP_LYFO %>% 
  clean_RKKP_LYFO()

LYFO_clean_studypop <- LYFO_clean %>%
  filter(patientid %in% study_population$patientid)
#Saving all LYFO data on the cohort
#fwrite(LYFO_clean_studypop, "CD20_LYFO_data_20251203.txt", sep = "\t")  #03.121.2025 all LYFO data is saved on the cohort


#1L_Chemo####
setwd("/ngc/projects2/dalyca_r/sangho_r/data/CD20_final_data_SG")
chemo <- fread("chemo_1stline.txt")
sum(is.na(chemo$treatment_type_1st_line)) #So 1L therapy registered for all patients
n_distinct(chemo) 

#Baseline####
names(LYFO_clean_studypop)

#Gattering data:
study_population1 <- study_population %>%
  left_join(patient %>% select(patientid, sex, date_birth), by = "patientid") %>%
  left_join(LYFO_clean_studypop %>% 
              select(patientid, subtype, date_diagnosis, age_diagnosis,date_response_1st_line, date_ASCT_2nd_line,LDH_elevated_diagnosis, 
                     n_extranodal_regions_diagnosis, AA_stage_diagnosis, PS_diagnosis, PS_2nd_line,response_1st_line, response_2nd_line, date_death_fu,dead, 
                     date_treatment_2nd_line, date_treatment_1st_line), 
            by = "patientid") %>%
  left_join(IPI %>% select(patientid, IPI, subtype = Disease), by = c('patientid', 'subtype')) %>% 
  left_join(chemo %>% select(patientid, treatment_type_1st_line), by ='patientid')%>% 
  mutate(
    date_birth = as_date(date_birth),
    date_death_fu = as_date(date_death_fu),
    age_at_2L = round(diff_years(date_birth,coalesce(date_treatment_2nd_line, date_ASCT_2nd_line)), 1))

colSums(is.na(study_population1))


#### PolyRx ####
# polyrx.sg.20250624 script made by CB:
setwd("/ngc/projects2/dalyca_r/sangho_r/data/")
DX.FIRST.POLY <- readRDS("DX_FIRST_POLY.rds")
study_population3 <- left_join(
  study_population1,
  DX.FIRST.POLY %>% select(patientid,n.ATC) %>% 
    mutate (ATC_category = case_when(n.ATC <5 ~ "<5",
                                     n.ATC >= 5 ~ "≥5")),
  by = "patientid"
)

#Sanity check
table(study_population3$filter,useNA = "always" )
n_distinct(study_population3$patientid) 

#respons####
sum(is.na(study_population3$response_1st_line))#only 3 NA
#grouping and recoding
#CMR (complete metabolik response) and CMR with residual mass recoded as CR (complete response)
#Unknow = NA AND "MORS" = NA
#CR:Complete Response(incl. CRi, CRu, CMR, osv.),PR:Partial Response (incl. PRi, PRu),SD:Stable Disease(incl. SDi, SDu),PD:Progressive Disease(incl.PDi,PDu,PD)
study_population4 <- study_population3 %>% 
  mutate(response = recode_factor(response_1st_line, CRi = "CR",
                                  CRu = "CR", Cru ="CR",
                                  "CMR with residual mass" = "CR","CMR and CR" = "CR",
                                  PRi = "PR",
                                  PRu = "PR",
                                  SDi = "SD",
                                  SDu = "SD",
                                  PDi = "PD",
                                  PDu = "PD",
                                  UNK = NA_character_,
                                  mors = NA_character_,))

#Making AA numeric so i can categorize
study_population4$AA_stage_diagnosis<- factor(study_population4$AA_stage_diagnosis,
                                              levels = 1:4,
                                              labels = c(" I", "II", "III", "IV"))
#Missing date_treatment####
# study_population4[is.na(study_population4$date_response_1st_line), ]

patientids_with_missing_first_line_treatment <- study_population4 %>% filter(
  is.na(date_treatment_1st_line)
) %>% pull(patientid)
new_dates_to_merge <- c(as.Date("2005-03-22"),
                        as.Date("2013-05-14"),
                        as.Date("2022-04-28"))
data_frame_to_merge <- data_frame(patientids_with_missing_first_line_treatment, 
                                  new_dates_to_merge) %>%
  select(
    patientid = patientids_with_missing_first_line_treatment,
    date_treatment_1st_line_new = new_dates_to_merge
  )

# merge new dates and then drop new date column - there's probably a cleaner way of doing this

study_population4 <- study_population4 %>%
  left_join(data_frame_to_merge) %>%
  mutate(date_treatment_1st_line = ifelse(is.na(date_treatment_1st_line), date_treatment_1st_line_new, date_treatment_1st_line)) %>%
  select(-date_treatment_1st_line_new)

#Time to relapse####
study_population5 <- study_population4 %>% 
  mutate(date_treatment_1st_line = as.Date(date_treatment_1st_line),
         date_treatment_2nd_line = if_else(is.na(date_treatment_2nd_line),
                                           as.Date(date_ASCT_2nd_line),
                                           as.Date(date_treatment_2nd_line)),
         TT2L = as.numeric(difftime(date_treatment_2nd_line,date_treatment_1st_line, units = "days"))/30.44)

summary(study_population5$TT2L)

#For DLBCL: within 12 month early relapse, MCL and FL 24
study_population6 <- study_population5 %>% 
  mutate(TT2L_group = case_when(
    subtype == "DLBCL" & TT2L <=12 ~ "POD12",
    subtype == "DLBCL" & TT2L > 12 ~ "Late",
    subtype == "FL" & TT2L <= 24 ~ "POD24",
    subtype == "FL" & TT2L > 24 ~ "Late",
    subtype == "MCL" & TT2L <=12 ~ "POD12",
    subtype == "MCL" & TT2L >12 ~ "Late",
    TRUE~NA_character_))

study_population6 [is.na(study_population6 $TT2L_group), ]


#Survival calculations####
study_population6$time_to_death <- as.numeric(difftime(as.Date(study_population6$date_death_fu), 
                                            as.Date(study_population6$date_diagnosis), unit = "days")) 
study_population6 %>% count(CD20)

colSums(is.na(study_population6))

setwd("/ngc/projects2/dalyca_r/sangho_r/data/CD20_final_data_SG")
#fwrite(study_population6, "CD20_df.txt", sep = "\t")
