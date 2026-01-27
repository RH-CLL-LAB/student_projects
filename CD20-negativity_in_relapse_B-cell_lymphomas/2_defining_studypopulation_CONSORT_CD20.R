# title: "CD20-data_Studypopulation"
# output: html_document
# date: "2025-09-15"
# Author: Sanaz M. Gholy, Jojo Biel-Nielsen Dietz 
# ---
#   
library(data.table)
library(tidyverse)
library(dplyr)
source('/ngc/projects2/dalyca_r/clean_r/load_dalycare_package.R') 
load_dataset(c('t_dalycare_diagnoses', 'patient', 'RKKP_LYFO')) 

LYFO_clean = RKKP_LYFO %>% 
  clean_RKKP_LYFO() 

#diagnoser i LYFO + antal
LYFO_clean$subtype %>% table

#Filtrer på sybtyper af interesse
LYFO_subtypes <- LYFO_clean %>% 
  filter(subtype %in% c('DLBCL', 'FL', 'MCL')) 
table(LYFO_subtypes$subtype, useNA = "always")
nrow(LYFO_subtypes)

#Removing CNS-DLBCL based on their 1L treatment
CNS_pt <- fread("/ngc/projects2/dalyca_r/sangho_r/data/chemo_CNS.txt")  %>% 
  filter(treatment_type == "CNS")
head(CNS_pt)
nrow(CNS_pt)

LYFO_subtype_uCNS <- LYFO_subtypes %>%
  filter(!patientid %in% CNS_pt$patientid)

table(LYFO_subtype_uCNS$subtype, useNA = "always")

#Only those recieving rituximab at first line
LYFO_subtypes_rituximab <- LYFO_subtype_uCNS %>%
  filter(immunotherapy_type_1st_line =='rituximab')
table(LYFO_subtypes_rituximab$subtype, useNA = "always")
nrow(LYFO_subtypes_rituximab)

#All relapsed patients who received secondline 
LYFO_subtypes_rituximab_relaps <- LYFO_subtypes_rituximab %>%
  filter(!is.na(as.character(date_treatment_2nd_line)) | 
        !is.na(as.character(date_ASCT_2nd_line)))

table(LYFO_subtypes_rituximab_relaps$subtype, useNA = "always")

nrow(LYFO_subtypes_rituximab_relaps)
n_distinct(LYFO_subtypes_rituximab_relaps$patientid)

#------------------------------------------------

#Patologi reports on this population => SDS_pato_subset
load_dataset('SDS_pato', LYFO_subtypes_rituximab_relaps$patientid)
nrow(SDS_pato_subset)
head(SDS_pato_subset)
#n = 247265 (OBS more rows each report)

pato_patientid <- SDS_pato_subset %>%
  select(k_inst, k_rekvnr, patientid) %>%
  distinct()
nrow(pato_patientid)
n_distinct(pato_patientid$patientid)

#more rows each report, fixing that issue here
load_dataset('SDS_t_mikro_ny', value = SDS_pato_subset$k_rekvnr, column = 'k_rekvnr') # load pato-text, specifically for k_rekvnr in CD20_cohort 
pato_text_clean <- SDS_t_mikro_ny_subset %>% 
  clean_SDS_t_mikro() 
pato_svar <- pato_patientid %>%
  left_join(pato_text_clean, by = c("k_inst", "k_rekvnr")) %>%
  select(!text)
  
head(pato_svar)
n_distinct(pato_svar$patientid)

population_subtype <- LYFO_subtypes_rituximab_relaps %>%
  select(patientid, subtype, date_treatment_2nd_line, date_ASCT_2nd_line)

pato_clean_subtype <- pato_svar %>%
  left_join(population_subtype, by = c("patientid"))

nrow(pato_clean_subtype)
head(pato_clean_subtype)
n_distinct(pato_clean_subtype$patientid)


table(pato_clean_subtype$subtype)
#DLBCL = 19402, FL = 8472, MCL = 6123

#join again, missing some values
manglende_info <-SDS_pato_subset %>%
  group_by(k_inst, k_rekvnr, patientid) %>%
  slice(1) %>%
  ungroup() %>%
  select(k_inst, k_rekvnr, patientid, d_rekvdato, c_snomedkode)
head(manglende_info)

pato_clean_subtype_2 <- pato_clean_subtype %>%
  left_join(manglende_info, by = c("k_inst", "k_rekvnr", "patientid"))
head(pato_clean_subtype_2)
nrow(pato_clean_subtype_2)

#Excluding remission or obs pro
pato_clean_subtype_3 <- pato_clean_subtype_2 %>% 
  dplyr::rename(date_pato = d_rekvdato)  %>%
  filter(!str_detect(c_snomedkode, "8$"),
         !str_detect(c_snomedkode, "X$")) # exclude remission and obs pro
table(pato_clean_subtype_3$subtype) 
nrow(pato_clean_subtype_3)
n_distinct(pato_clean_subtype_3$patientid)


names(pato_clean_subtype_3)
#fwrite(pato_clean_subtype_3, "pato_datoer_20251124.txt", sep="\t")
getwd()
#problem with those without 2.line treatment date, but ASCT.
#2.line is priorotized if they have both 
test <- pato_clean_subtype_3 %>%
  filter(is.na(date_treatment_2nd_line))

test2 <- pato_clean_subtype_3 %>%
  mutate(date_treatment_2nd_line = 
           if_else(is.na(date_treatment_2nd_line), date_ASCT_2nd_line, date_treatment_2nd_line))

test3 <- test2 %>%
  filter(is.na(date_treatment_2nd_line))

#only patologi report 60 dys before and 30 days after 2. line treatment date:
pato_clean_subtype_4 <- test2 %>%
  mutate(time_2L_pato=diff_days(date_treatment_2nd_line,date_pato)) %>% 
  filter(time_2L_pato<=30,
         time_2L_pato>-60) %>%
  mutate(kohorte = "yes")

#Are ASCT patients stille there after 30/60 days filtring?
setwd('/ngc/projects2/dalyca_r/sangho_r/data/split_text_output/')
chunk10 <- read_csv2("pato_3ASCT_20250903.csv")

normalize_id <- function(x){
  x %>% as.character() %>% str_trim() %>% str_replace("^0+", "")}

asct <- chunk10 %>%
  dplyr::rename(patientid_1 = k_inst,
         k_rekvnr_1 = date_pato) %>%
  select(patientid_1, k_rekvnr_1) %>%
  dplyr:: rename(patientid = patientid_1,
         k_rekvnr = k_rekvnr_1) %>%
  mutate(k_rekvnr = normalize_id(k_rekvnr))
head(asct)

asct_pato <- pato_clean_subtype_4 %>%
  mutate(k_rekvnr = normalize_id(k_rekvnr))

asct3 <- asct_pato %>%
  semi_join(asct, by = c("patientid", "k_rekvnr"))

ukendt_date <- pato_clean_subtype_4 %>% 
  filter(is.na(date_treatment_2nd_line))
  
table(pato_clean_subtype_4$subtype) 
nrow(pato_clean_subtype_4)
head(pato_clean_subtype_4)
n_distinct(pato_clean_subtype_4$patientid)

subtype <- pato_clean_subtype_4 %>%
  group_by(patientid) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(patientid = as.numeric(patientid))
table(subtype$subtype)
nrow(subtype)

setwd("/ngc/projects2/dalyca_r/sangho_r/data/")
getwd()
#fil med cd20-textminede patologi raporter
cd20 <- fread("pato_text_cd20_clean.tsv") %>% 
  select(!c("V1")) %>%
  mutate(CD20 = "CD20")
nrow(cd20)

summary(cd20)
n_distinct(cd20$patientid)

#filtring based on textminede patology-reports
pato_clean_subtype_5 <- pato_clean_subtype_4 %>%
  left_join(cd20, by = c("patientid", "k_inst", "k_rekvnr")) %>%
  filter(CD20 == "CD20")
nrow(pato_clean_subtype_5)
table(pato_clean_subtype_5$subtype)

n_distinct(pato_clean_subtype_5$patientid)
names(pato_clean_subtype_5)

er_asct_der_stadig <- pato_clean_subtype_5 %>%
  semi_join(asct3, by = c("patientid", "k_rekvnr"))

# 2664 nots in total (2011 unique)
#Sanaz went through 1343 to many reports. We wil filter these out when making file for Torgerod

pato_clean_subtype_6 <- pato_clean_subtype_5 %>%
  select(patientid, k_rekvnr, k_inst, CD20)
getwd()
#fwrite(pato_clean_subtype_5,"brug_disse_patobeskrivelser_full_20250915.txt", sep = "\t")
#fwrite(pato_clean_subtype_6,"brug_disse_patobeskrivelser_jd_20250915.txt", sep = "\t")

subtype2 <- pato_clean_subtype_5 %>%
  group_by(patientid) %>%
  slice(1) %>%
  ungroup()
nrow(subtype2)

table(subtype2$subtype)

#Study populationen
pato_beskrivelser <- fread("brug_disse_patobeskrivelser_20250915.txt")
names(pato_beskrivelser)
head(pato_beskrivelser)

pato_beskrivelser2 <- pato_beskrivelser %>%
  select(!CD20)
asct4 <- asct3 %>%
  select(patientid, k_rekvnr, k_inst)

pato_beskrivelser3 <- rbind(pato_beskrivelser2, asct4)

#Patobeskrivelser + ASCT
#fwrite(pato_beskrivelser3, "brug_disse_patobeskrivelser_20250916.txt")


n_distinct(pato_beskrivelser$patientid)
nrow(pato_beskrivelser)
head(pato_beskrivelser)

study_kohorte <- pato_beskrivelser %>%
  group_by(patientid) %>%
  slice(1) %>%
  ungroup() %>%
  select(patientid)

#fwrite(study_kohorte, "studie_population_id_20250915.txt", sep = "\t")

