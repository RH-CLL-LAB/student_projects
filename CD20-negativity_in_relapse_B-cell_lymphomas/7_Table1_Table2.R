# title: "CD20-data_Table 1 and 2"
# output: html_document
# last revised date: "2026-01-19"
# Author: Sanaz M. Gholy
#---

source('/ngc/projects2/dalyca_r/clean_r/load_dalycare_package.R')
library(data.table)
library(tidyverse)
library(Publish)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(scales)
library(table1)
library(tableone)


#Tabel1####

setwd("/ngc/projects2/dalyca_r/sangho_r/data/CD20_final_data_SG")

df <- fread("CD20_df.txt")
n_distinct(df)
names(df)


df <- df %>% 
  mutate(extranodal_group = case_when(
    n_extranodal_regions_diagnosis== 0 ~ "0",
    n_extranodal_regions_diagnosis== 1 ~ "1",
    n_extranodal_regions_diagnosis >1 ~ ">1", 
    TRUE ~NA_character_))

df <- df %>% 
  mutate(PS_diagnosis = case_when(
    PS_diagnosis== 0 ~ "0",
    PS_diagnosis== 1 ~ "1",
    PS_diagnosis >1 ~ ">1", 
    TRUE ~NA_character_))

df <- df %>% 
  mutate(IPI = ifelse(is.na(IPI), "Missing", IPI),
         AA_stage_diagnosis = ifelse(AA_stage_diagnosis== "", "Missing", AA_stage_diagnosis),
         AA_stage_diagnosis = ifelse(AA_stage_diagnosis == "Stage 1", "I", AA_stage_diagnosis),
         AA_stage_diagnosis = ifelse(AA_stage_diagnosis == "Stage 2", "II", AA_stage_diagnosis),
         AA_stage_diagnosis = ifelse(AA_stage_diagnosis == "Stage 3", "III", AA_stage_diagnosis),
         AA_stage_diagnosis = ifelse(AA_stage_diagnosis == "Stage 4", "IV", AA_stage_diagnosis))


df$treatment_type_1st_line<- recode(df$treatment_type_1st_line,
                                    "CHOP + HD-MTX" = "CHOP-like",
                                    "Bio-CHIC" = "CHOP-like",
                                    "EPOCH-R" = "CHOP-like",
                                    "CHOP" = "CHOP-like",
                                    "MCL" = "R-maxi-chop/TRIANGLE",
                                    "maxi-CHOP" = "R-maxi-chop/TRIANGLE")
unique(df$treatment_type_1st_line)

df$subtype <- factor(df$subtype, levels = c("DLBCL", "FL", "MCL")) 
df$sex <- factor(df$sex, levels = c("F", "M"), labels = c("Female", "Male")) 
df$AA_stage_diagnosis <- factor(df$AA_stage_diagnosis, levels = c("I", "II", "III", "IV"))
df$IPI <- factor(df$IPI, levels = c("Low", "Intermediate", "High"))   
df$response <- factor(df$response, levels = c("CR", "PR", "SD", "PD"))
df$treatment_type_1st_line <- factor(df$treatment_type_1st_line, levels = c("CHOP-like", "Bendamustine/BAC", "R-monotherapy", "R-maxi-chop/TRIANGLE","CVP", "Other"))
df$ATC_category <- factor(df$ATC_category, levels = c("<5", "≥5"))
df$CD20 <- factor(df$CD20, levels = c("Positive", "Negative"))
df$LDH_elevated_diagnosis <- factor(df$LDH_elevated_diagnosis, levels = c("Yes", "No"))
df$extranodal_group <- factor(df$extranodal_group, levels = c("0", "1", ">1"))
df$PS_diagnosis <- factor(df$PS_diagnosis, levels = c("0", "1", ">1"))

label(df$AA_stage_diagnosis) <- "Ann Arbor stage at diagnosis"
label(df$age_diagnosis) <- "Age at diagnosis"
label(df$age_at_2L) <- "Age at second-line treatment"
label(df$IPI) <- "IPI score at diagnosis"
label(df$treatment_type_1st_line)<- "First-line treatment regimen"
label(df$response)<- "Response to first-line treatment"
label(df$sex)<- "Sex"
label(df$ATC_category)<- "Polypharmacy at diagnosis"
label(df$LDH_elevated_diagnosis) <- "High LDH level at diagnosis"
label(df$extranodal_group) <- "Number of extranodal sites at diagnosis"
label(df$PS_diagnosis) <- "ECOG at diagnosis"

n_distinct(df$patientid)
#Tabel 1####
table1 <- table1(~ age_diagnosis +age_at_2L+ sex + AA_stage_diagnosis +LDH_elevated_diagnosis+extranodal_group+ PS_diagnosis+ IPI +
                   treatment_type_1st_line + ATC_category 
                 | subtype * CD20, data = df, render.missing = NULL, 
                 render.continuous = c(.="Median [Q1;Q3]"), render.categorical="FREQ (PCTnoNA%)")
table1

#Tabel2####
table2total<- CreateTableOne(vars <- c("response", "TT2L"),
                         strata = "CD20",
                         data = df,
                         test= TRUE, 
                         factorVars =c("sex", "AA_stage_diagnosis", "IPI", 
                                       "response", "ATC_category") )



print(table2total, showAllLevels = TRUE, test = TRUE, exact = "factorVars", quote = FALSE, noSpaces = TRUE)

#Tabel 2, resposne og TT2L

Response_table <- df %>% 
  count(subtype, CD20, response) %>% 
  group_by(subtype, CD20) %>% 
  mutate(precent= round (100*n/sum(n),1),
         n_precent= paste0(n, "(", precent,"%)")) %>% 
  select(subtype, CD20, response, n_precent) %>% 
  pivot_wider(names_from= response,
              values_from = n_precent) %>% 
  ungroup()


TT2L_table <- df %>%
  group_by(subtype, CD20) %>% 
  summarise(median=median(TT2L, na.rm = TRUE),
                    q1= quantile(TT2L, 0.25, na.rm = TRUE),
                    q3= quantile(TT2L, 0.75, na.rm = TRUE),
                    median_IQR=paste0(median, "(", q1, "-", q3, ")"),
            .groups = "drop")

p_val_resp <- function(df){
  tab <- table(df$CD20, df$response)
  if(any(tab<5)){fisher.test(tab)$p.value} else{chisq.test(tab)$p.value}}

p_val_TT2L <- function(df){
  wilcox.test(TT2L ~ CD20, data= df)$p.value}

pval_table <- df %>% 
  group_by(subtype) %>% 
  summarise(p_reponse =p_val_resp(cur_data()),
            p_TT2L= p_val_TT2L (cur_data()),
            .groups= "drop")

Table2 <-  Response_table %>% 
  left_join(TT2L_table, by= c("subtype", "CD20")) %>% 
  left_join(pval_table, by= "subtype")

Table2


#TT2L sensi####
#CR og PD______________________________________________________
CR_df <- df[df$response == "CR", c("patientid", "TT2L", "CD20", "subtype", "response")]
PD_df <- df[df$response == "PD", c("patientid", "TT2L", "CD20", "subtype", "response")]

label(CR_df$TT2L) <- "Time to 2nd-line therapy"
units(CR_df$TT2L) <- "months"

tableCR<- CreateTableOne(vars <- c( "subtype",  
                                    "response", "TT2L"),
                         strata = "CD20",
                         data = CR_df,
                         test= TRUE, 
                         factorVars =c( "response", "subtype") )
print(tableCR, showAllLevels = TRUE, test = TRUE, exact = "factorVars", quote = FALSE, noSpaces = TRUE)

PD_df$response <- factor(PD_df$response, levels = c("PD"))
label(PD_df$TT2L) <- "Time to 2nd-line therapy"
units(PD_df$TT2L) <- "months"


tablePD<- CreateTableOne(vars <- c( "subtype",  
                                    "response", "TT2L"),
                         strata = "CD20",
                         data = PD_df,
                         test= TRUE, 
                         factorVars =c( "response", "subtype") )
print(tablePD, showAllLevels = TRUE, test = TRUE, exact = "factorVars", quote = FALSE, noSpaces = TRUE)
