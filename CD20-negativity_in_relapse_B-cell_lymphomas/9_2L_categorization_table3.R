# title: "CD20-data_2nd line therapy"
# output: html_document
# date: "2025-12-18"
# Author: Sanaz M. Gholy, Christian Brieghel 
#---


library(data.table)
library(tidyverse)
library(dplyr)
library(tidyr)
library(arrow)
library(lubridate)

setwd("/ngc/projects2/dalyca_r/sangho_r/data/CD20_final_data_SG")
cohort2 <- fread("CD20_df.txt")%>% 
  select(patientid, subtype, date_treatment_2nd_line, CD20, date_death_fu)

source('/ngc/projects2/dalyca_r/clean_r/load_dalycare_package.R') 
load_dataset(c('t_dalycare_diagnoses', 'patient', 'RKKP_LYFO' ))

RKKP_LYFO <- readRDS('/ngc/projects2/dalyca_r/andkat_r/RT/Data/RKKP_LYFO_oct25.rds') #RKKP is being updated
LYFO_clean = RKKP_LYFO %>% 
  clean_RKKP_LYFO()
names(LYFO_clean)
LYFO_clean %>% select(matches('2nd')) %>% names

#2L####
cohort3 <- cohort2%>%
  left_join(LYFO_clean %>% 
  select(patientid, chemo_regime_1_type_2nd_line, 
         chemo_regime_2_type_2nd_line, 
         chemo_regime_3_type_2nd_line,
         immunotherapy_type_2nd_line,
         chop_like_2nd_line,
         high_dosis_treatment_2nd_line,
         maintenance_treatment_initiated_2nd_line,
         date_RT_2nd_line, #to check if they recieved Radiation therapy
         RT_dosis_Gy_2nd_line,
         RT_type_2nd_line,
         ASCT_2nd_line,
         steroid_monotherapy_2nd_line,
         RT_dosis_mCkg_2nd_line), by='patientid')

cohort3$immunotherapy_type_2nd_line %>% unique %>% sort() 

cohort3$chemo_regime_1_type_2nd_line%>% unique %>% sort() 
#55 

cohort3$chemo_regime_2_type_2nd_line%>% unique %>% sort()
#43 

cohort3$chemo_regime_3_type_2nd_line%>% unique %>% sort()
#25
cohort3$chop_like_2nd_line

summary(cohort3$maintenance_treatment_initiated_2nd_line)


#How many NA?
cohort3 %>% 
  summarise(NA_immunotherapy= sum(is.na(immunotherapy_type_2nd_line)),
            NA_chemo1= sum(is.na(chemo_regime_1_type_2nd_line)),
            NA_chemo2= sum(is.na(chemo_regime_2_type_2nd_line)),
            NA_chemo3= sum(is.na(chemo_regime_3_type_2nd_line)))


#Categorising_chemoregime####

text = cohort3$treatment_type %>% unique()
text = gsub(' NA','',  text)
text %>% n_distinct()

# %>% paste0(collapse = "` = '', `")

#OBS- minibeam --> salvage not HDT

cohort3$high_dosis_treatment_2nd_line %>% table
cohort4 <- cohort3 %>% 
  mutate(antiCD20_2L = ifelse(str_detect(immunotherapy_type_2nd_line, str_flatten(c('ritux', 'humax', 'ofat', 'obin'), '|')),
                                         'Yes', 'No')) %>% 
  mutate(treatment_type_2 =paste(chemo_regime_1_type_2nd_line,
                               chemo_regime_2_type_2nd_line,
                               chemo_regime_3_type_2nd_line,
                               immunotherapy_type_2nd_line),
         treatment_type_2nd_line = recode(treatment_type_2,
                                 `bendamustin ice dhap rituximab` = 'Salvage', 
                                 `ccvp NA NA rituximab` = 'Palliative',
                                 `bendamustin NA NA rituximab` = 'BR', 
                                 `ice NA NA NA` = 'Salvage',
                                 `ice minibeam hypercvad rituximab` = 'Salvage',
                                 `chop NA NA rituximab` = 'R-CHOP',
                                 `chop cope NA rituximab` = 'R-CHOP',
                                 `dhap NA NA rituximab` = 'Salvage',
                                 `dhap ice NA rituximab` = 'Salvage',
                                 `ice beam NA NA` = 'Salvage',
                                 `nordiskcns temozolomid NA NA` = 'CNS',
                                 `abvd_copp NA NA rituximab` = 'Other',
                                 `ice NA NA rituximab` = 'Salvage',
                                 `chop NA NA NA` = 'CHOP',
                                 `ice dhap NA rituximab` = 'Salvage', 
                                 `dhap vim NA rituximab` = 'Salvage',
                                 `bendamustin cope gemox rituximab` = 'Palliative', 
                                 `hdmtx NA NA NA` = 'HD-MTX', 
                                 `bendamustin NA NA humax` = 'BR',
                                 `dhap NA NA humax` = 'Salvage',
                                 `cvp NA NA rituximab` = 'R-CVP',
                                 `NA NA NA rituximab` = 'anti_CD20', 
                                 `oncovin hdmtx NA rituximab` = 'HD-MTX', 
                                 `cns_matrix ice NA rituximab` = 'CNS', 
                                 `chop ceop NA rituximab` = 'R-CHOP', 
                                 `ice peben gemcitabin rituximab` = 'Salvage',
                                 `lenalidomid NA NA rituximab` = 'Targeted-anti_CD20', #R2
                                 `ice beam NA rituximab` = 'Salvage', 
                                 `NA NA NA NA` = 'NA', 
                                 `temozolomid NA NA rituximab` = 'CNS', 
                                 `chop hdmtx NA rituximab` = 'R-CHOP',
                                 `chop cvp NA rituximab` = 'R-CHOP', 
                                 `peben hdmtx NA rituximab` = 'Other', 
                                 `cvbp NA NA rituximab` = 'Palliative', 
                                 `choep NA NA rituximab` = 'R-CHOP', 
                                 `ccvp NA NA NA` = 'Palliative', 
                                 `cnsbonn NA NA rituximab` = 'CNS',
                                 `chop dhap minibeam rituximab` = 'Salvage', 
                                 `gemox NA NA rituximab` = 'Palliative', 
                                 `chop dhap NA rituximab` = 'Salvage', 
                                 `ibrutinib chop NA rituximab` = 'R-CHOP', 
                                 `cope chop NA rituximab` = 'R-CHOP', 
                                 `mvbpcns NA NA rituximab` = 'CNS', 
                                 `dhap NA NA ofatumumab` = 'Salvage',
                                 `epoch beam NA rituximab` = 'Salvage', 
                                 `ice itbehandling NA rituximab` = 'Salvage',
                                 `ice dhap minibeam rituximab` = 'Salvage', 
                                 `bendamustin cvp NA rituximab` = 'R-CVP',
                                 `cope itbehandling NA rituximab` = 'R-CHOP', 
                                 `cop NA NA rituximab` = 'R-CHOP', 
                                 `dhap beam NA rituximab` = 'Salvage', 
                                 `fcd NA NA rituximab` = 'Other', 
                                 `ice hdmtx NA NA` = 'Salvage', 
                                 `peben NA NA rituximab` = 'Other',
                                 `ice NA NA obinutuzumab` = 'Salvage', 
                                 `cns_mvp NA NA rituximab` = 'CNS',
                                 `dhap hdmtx ice rituximab` = 'Salvage',
                                 `bendamustin chop NA ofatumumab` = 'R-CHOP',
                                 `dhap NA NA NA` = 'Salvage', 
                                 `maxichop hdarac NA rituximab` = 'MCL', #Skal det hedde MCL?
                                 `chop cope hdmtx rituximab` = 'R-CHOP',
                                 `venetoclax NA NA obinutuzumab` = 'Targeted-anti_CD20', 
                                 `chop choep NA rituximab` = 'R-CHOP', 
                                 `cvp bendamustin NA rituximab` = 'R-CVP',
                                 `mime NA NA rituximab` = 'Other',
                                 `ceop NA NA rituximab` = 'R-CHOP', 
                                 `cope NA NA rituximab` = 'R-CHOP', 
                                 `ice bendamustin idelalisib NA` = 'Salvage',
                                 `dhap gemcitabin NA rituximab` = 'Salvage',
                                 `choep hdmtx NA NA` = 'R-CHOP',
                                 `gdp hdmtx NA rituximab` = 'Salvage', 
                                 `gemox gemcitabin NA rituximab` = 'Palliative',
                                 `NA NA NA glofitamab` = 'BsAbs', 
                                 `gemcitabin NA NA NA` = 'Other', 
                                 `bfm NA NA rituximab` = 'Other',
                                 `hdmtx ifos NA rituximab` = 'CNS',
                                 `ice itbehandling cns_matrix rituximab` = 'CNS', 
                                 `velcade NA NA rituximab` = 'Targeted-anti_CD20', 
                                 `ibrutinib NA NA NA` = 'Targeted', 
                                 `hdarac bendamustin NA rituximab` = 'R-BAC', 
                                 `chop hdmtx NA obinutuzumab` = 'R-CHOP',
                                 `dhap ibrutinib NA rituximab` = 'Salvage',
                                 `venetoclax lenalidomid NA rituximab` = 'Targeted-anti_CD20', # VALERIA
                                 `bendamustin velcade NA rituximab` = 'BR', # + velcade
                                 `bendamustin ice NA rituximab` = 'Salvage', 
                                 `peroral NA NA NA` = 'Palliative',
                                 `ibrutinib itbehandling NA NA` = 'Targeted',
                                 `ibrutinib NA NA rituximab` = 'Targeted-anti_CD20', 
                                 `pix-etoposid-benda hdmtx NA rituximab` = 'Other', 
                                 `bac NA NA rituximab` = 'R-BAC',
                                 `venetoclax NA NA ofatumumab` = 'Targeted-anti_CD20', 
                                 `ibrutinib lenalidomid NA rituximab` = 'Targeted-anti_CD20',
                                 `chlorambucil NA NA rituximab` = 'Palliative',
                                 `mantle3 NA NA rituximab` = 'MCL',
                                 `ice hdarac NA rituximab` = 'Salvage',
                                 `fcd mitoxantrone NA rituximab` = 'Other',
                                 `lenalidomid ibrutinib NA NA` = 'Targeted',
                                 `NA NA NA humax` = 'anti_CD20', 
                                 `bendamustin ibrutinib NA rituximab` = 'BR',
                                 `chop hdarac NA rituximab` = 'R-CHOP', 
                                 `ibrutinib bendamustin NA rituximab` = 'BR', 
                                 `bendamustin NA NA NA` = 'Other', 
                                 `thalidomid NA NA rituximab` = 'Targeted-anti_CD20', 
                                 `epoch NA NA rituximab` = 'R-CHOP', 
                                 `gdp NA NA rituximab` = 'Salvage', 
                                 `hdmtx NA NA rituximab` = 'HD-MTX', 
                                 `cns_matrix NA NA NA` = 'CNS', 
                                 `bfm ice minibeam rituximab` = 'Salvage', 
                                 `chlorambucil NA NA NA` = 'Palliative', 
                                 `dhap NA NA glofitamab` = 'BsAbs', 
                                 `ccvp ice NA rituximab` = 'Palliative',
                                 `chop itbehandling NA rituximab` = 'R-CHOP',
                                 `itbehandling temozolomid NA NA` = 'CNS',
                                 `dhap ice hdarac rituximab` = 'Salvage', 
                                 `bendamustin ice hdctx rituximab` = 'BR',
                                 `eshap NA NA NA` = 'Other',
                                 `bendamustin bendamustin NA rituximab` = 'BR',
                                 `gemox NA NA glofitamab` = 'BsAbs',
                                 `choep chop NA rituximab` = 'R-CHOP',
                                 `itbehandling NA NA NA` = 'CNS', 
                                 `ice dhap NA NA` = 'Salvage', 
                                 `fcd NA NA NA` = 'Other', 
                                 `dhap cope NA rituximab` = 'Salvage',
                                 `dhap beam NA obinutuzumab` = 'Salvage',
                                 `chop NA NA glofitamab` = 'BsAbs', 
                                 `peben NA NA NA` = 'Other',
                                 `gdp chop NA rituximab` = 'Salvage',
                                 `cns_matrix NA NA rituximab` = 'CNS', 
                                 `dhap minibeam NA rituximab` = 'Salvage', 
                                 `cns_matrix chop NA rituximab` = 'CNS', 
                                 `bendamustin NA NA ofatumumab` = 'BR', 
                                 `hdmtx ice NA rituximab` = 'HD-MTX', 
                                 `bendamustin NA NA obinutuzumab` = 'BR',
                                 `dhap ice NA NA` = 'Salvage', 
                                 `ice minibeam ccvp rituximab` = 'Salvage', 
                                 `oncovin NA NA rituximab` = 'Palliative',
                                 `dhap ice lenalidomid rituximab` = 'Salvage',
                                 `minichop NA NA rituximab` = 'R-CHOP', 
                                 `cyclofosfamid NA NA rituximab` = 'Palliative', 
                                 `cope NA NA NA` = 'CHOP', 
                                 `cns_matrix ice bcnu rituximab` = 'CNS', 
                                 `temozolomid NA NA NA` = 'CNS', 
                                 `cns_mvp temozolomid NA rituximab` = 'CNS', 
                                 `cns_matrix temozolomid NA rituximab` = 'CNS', 
                                 `chop NA NA obinutuzumab` = 'R-CHOP', 
                                 `peben hdmtx NA NA` = 'Other', 
                                 `dhap venetoclax NA ofatumumab` = 'Salvage', 
                                 `choep NA NA ofatumumab` = 'R-CHOP', 
                                 `NA NA NA pdl1` = 'Targeted', 
                                 `dhap minibeam ice rituximab` = 'Salvage',
                                 `oncovin ice NA rituximab` = 'Salvage', 
                                 `cyclofosfamid doxorubicin hdarac rituximab` = 'R-CHOP', 
                                 `bendamustin hdctx NA rituximab` = 'BR', 
                                 `bendamustin gemcitabin NA rituximab` = 'BR', 
                                 `lenalidomid NA NA glofitamab` = 'BsAbs', 
                                 `gemcitabin NA NA rituximab` = 'Other', 
                                 `cop NA NA NA` = 'CHOP',
                                 `hdctx itbehandling hdarac rituximab` = 'CNS',
                                 `hdctx fcd NA rituximab` = 'Other', 
                                 `interferon velcade NA humax` = 'Targeted-anti_CD20', 
                                 `mvbpcns itbehandling NA NA` = 'CNS',
                                 `gdp peben NA obinutuzumab` = 'Salvage', 
                                 `pixantrone gemcitabin NA rituximab` = 'Other',
                                 `dhap pixantrone NA rituximab` = 'Salvage', 
                                 `dhap ice minibeam humax` = 'Salvage',
                                 `maximime NA NA rituximab` = 'Other', 
                                 `oncovin bendamustin NA rituximab` = 'BR',
                                 `gdp NA NA NA` = 'Salvage', 
                                 `cnsbonn NA NA NA` = 'CNS', 
                                 `ice hdmtx NA rituximab` = 'Salvage',
                                 `bendamustin cyclofosfamid NA rituximab` = 'BR',
                                 `hdmtx hdarac NA NA` = 'HD-MTX', 
                                 `cns_mvp NA NA NA` = 'CNS', 
                                 `choep chop hdarac obinutuzumab` = 'R-CHOP', 
                                 `cns_matrix ice lenalidomid rituximab` = 'CNS', 
                                 `hdmtx ifos temozolomid NA` = 'CNS', 
                                 `dhap peben NA rituximab` = 'Salvage', 
                                 `mime NA NA NA` = 'Other', 
                                 `cvbp NA NA NA` = 'Palliative', 
                                 `abvd ccvp NA rituximab` = 'Other',
                                 `nordiskcns NA NA NA` = 'CNS', 
                                 `chop cyclofosfamid NA rituximab` = 'R-CHOP', 
                                 `ice hdarac dhap rituximab` = 'Salvage', 
                                 `chop chop NA rituximab` = 'R-CHOP',
                                 `ice minibeam codoxmivac rituximab` = 'Salvage', 
                                 `ice minibeam NA rituximab` = 'Salvage', 
                                 `abvd NA NA brentuximab` = 'Other', 
                                 `venetoclax NA NA rituximab` = 'Targeted-anti_CD20',
                                 `dhap vim NA NA` = 'Salvage',
                                 `cns_mvp hdarac NA rituximab` = 'CNS', 
                                 `abvd NA NA NA` = 'Other', 
                                 `itbehandling hdmtx NA NA` = 'CNS', 
                                 `gdp ice NA rituximab` = 'Salvage', 
                                 `pixantrone NA NA rituximab` = 'Other', 
                                 `cns_matrix NA NA humax` = 'CNS',
                                 `cyclofosfamid NA NA NA` = 'Palliative',
                                 `chop minichop NA rituximab` = 'R-CHOP', 
                                 `lenalidomid doxorubicin cop ofatumumab` = 'R-CHOP',
                                 `bendamustin choep NA rituximab` = 'R-CHOP',
                                 `mopp itbehandling hdctx rituximab` = 'Other',
                                 `dhap minibeam NA NA` = 'Salvage', 
                                 `dhap itbehandling NA rituximab` = 'Salvage', 
                                 `ice ccvp NA rituximab` = 'Salvage',
                                 `dhap velbe NA rituximab` = 'Salvage', 
                                 `choep dhap bendamustin rituximab` = 'Salvage', 
                                 `gdp cope NA rituximab` = 'Salvage', 
                                 `dhap hdarac NA rituximab` = 'Salvage', 
                                 `hdmtx temozolomid NA rituximab` = 'CNS',
                                 `dhap gemox gemcitabin rituximab` = 'Salvage', 
                                 `gdp ice hdmtx rituximab` = 'Salvage', 
                                 `ice bendamustin NA rituximab` = 'Salvage', 
                                 `hdmtx cns_mvp NA rituximab` = 'CNS', 
                                 `ice dhap choep rituximab` = 'Salvage',
                                 `chop ibrutinib NA NA` = 'CHOP',
                                 `dhap vim hdmtx rituximab` = 'Salvage', 
                                 `chop beam NA rituximab` = 'Salvage', 
                                 `chop itbehandling dhap rituximab` = 'Salvage',
                                 `ice ice NA rituximab` = 'Salvage', 
                                 `ice vim cyclofosfamid rituximab` = 'Salvage', 
                                 `chop gdp NA rituximab` = 'Salvage', 
                                 `chop NA NA brentuximab` = 'CHOP', 
                                 `chl_vpp NA NA rituximab` = 'Palliative', 
                                 `dhap hdarac beam rituximab` = 'Salvage',
                                 `hdmtx ibrutinib NA rituximab` = 'HD-MTX',
                                 `dhap gemox NA rituximab` = 'Salvage',
                                 `minichop NA NA obinutuzumab` = 'R-CHOP', 
                                 `ice chop NA rituximab` = 'Salvage', 
                                 `hdarac itbehandling NA NA` = 'CNS',
                                 `ice dhap NA brentuximab` = 'Salvage', 
                                 `ice NA NA humax` = 'Salvage', 
                                 `cvp chop NA rituximab` = 'R-CHOP', 
                                 `dhap cyclofosfamid NA rituximab` = 'Salvage', 
                                 `ice hdctx beam rituximab` = 'Salvage', 
                                 `mbvdcns itbehandling hdmtx NA` = 'CNS',
                                 `NA NA NA polatuzumab` = 'Targeted', 
                                 `gemox ibrutinib lenalidomid rituximab` = 'Palliative',
                                 `gemox bendamustin NA rituximab` = 'Palliative',
                                 `gemox lenalidomid idelalisib rituximab` = 'Palliative',
                                 `cns_mvp ice itbehandling rituximab` = 'CNS', 
                                 `NA NA NA obinutuzumab` = 'anti_CD20', 
                                 `lenalidomid venetoclax NA rituximab` = 'Targeted-anti_CD20',
                                 `ibrutinib NA NA ofatumumab` = 'Targeted-anti_CD20',
                                 `lenalidomid ibrutinib NA obinutuzumab` = 'Targeted-anti_CD20',
                                 `adriamycin cyclofosfamid velcade rituximab` = 'Other', 
                                 `hdarac velcade hypercvad rituximab` = 'Other',
                                 `velcade oncovin NA rituximab` = 'Targeted-anti_CD20', 
                                 `velcade NA NA NA` = 'Targeted',
                                 `ibrutinib itbehandling NA rituximab` = 'Targeted-anti_CD20',
                                 `bac ibrutinib venetoclax rituximab` = 'R-BAC', 
                                 `chop hdarac ice rituximab` = 'Salvage', 
                                 `bac bendamustin NA rituximab` = 'R-BAC', 
                                 `peroral_cytostatika NA NA polatuzumab` = 'Targeted',
                                 `ibrutinib lenalidomid chop rituximab` = 'R-CHOP',
                                 `dhap ice peben rituximab` = 'Salvage',
                                 `pix-etoposid-benda NA NA rituximab` = 'Other',
                                 `gemox NA NA NA` = 'Palliative', 
                                 `chlorambucil NA NA ofatumumab` = 'Palliative',
                                 `ibrutinib NA NA glofitamab` = 'BsAbs',
                                 `idelalisib NA NA NA` = 'Targeted',
                                 `hdarac NA NA rituximab` = 'Other',
                                 `lenalidomid ibrutinib NA rituximab` = 'Targeted-anti_CD20',
                                 `bendamustin peroral_cytostatika NA rituximab` = 'BR',
                                 `choep NA NA NA` = 'CHOP',
                                 `cyclofosfamid NA NA ofatumumab` = 'Palliative',
                                 `hdarac velcade NA NA` = 'Targeted',
                                 `chop ibrutinib hdarac rituximab` = 'R-CHOP', 
                                 `hdarac hdmtx ifos rituximab` = 'CNS', 
                                 `fcd velcade bendamustin rituximab` = 'Other',
                                 `gdp NA NA glofitamab` = 'BsAbs', 
                                 `choep gdp beam rituximab` = 'Salvage', 
                                 `lenalidomid NA NA NA` = 'Targeted',
                                 `choep hdmtx NA rituximab` = 'R-CHOP',
                                 `maximime acvdl hdarac rituximab` = 'Other',
                                 `bendamustin NA NA brentuximab` = 'Other', 
                                 `bendamustin ceop NA rituximab` = 'BR', 
                                 `ice lenalidomid NA glofitamab` = 'BsAbs', 
                                 `dhap lenalidomid NA rituximab` = 'Salvage', 
                                 `gdp NA NA obinutuzumab` = 'Salvage', 
                                 `venetoclax NA NA NA` = 'Targeted', 
                                 `dhap ice beam rituximab` = 'Salvage', 
                                 `ccvp chlorambucil NA rituximab` = 'Palliative',
                                 `chop hdarac ibrutinib rituximab` = 'R-CHOP', 
                                 `hdmtx cyclofosfamid NA NA` = 'HD-MTX', 
                                 `ice cns_matrix itbehandling rituximab` = 'CNS', 
                                 `lenalidomid NA NA ofatumumab` = 'Targeted-anti_CD20',
                                 `dhap cisplatin itbehandling rituximab` = 'Salvage', 
                                 `cyclofosfamid doxorubicin NA rituximab` = 'R-CHOP',
                                 `chop NA NA polatuzumab` = 'CHOP', 
                                 `ibrutinib NA NA pdl1` = 'Targeted',
                                 `ice peben NA rituximab` = 'Salvage',
                                 `bendamustin hdmtx NA rituximab` = 'BR', 
                                 `itbehandling ibrutinib NA NA` = 'Targeted', 
                                 `dhap itbehandling ice rituximab` = 'Salvage', 
                                 `ice cyclofosfamid NA rituximab` = 'Salvage',
                                 `lvpp NA NA rituximab` = 'Other',
                                 `itbehandling cns_matrix NA NA` = 'CNS', 
                                 `cns_matrix bcnu NA rituximab` = 'CNS', 
                                 `chop ice NA rituximab` = 'Salvage',
                                 `bendamustin cnsbonn NA rituximab` = 'CNS', 
                                 `fludarabin hdarac NA NA` = 'Other',
                                 `dhap beam NA NA` = 'Salvage',
                                 `bendamustin chop beam rituximab` = 'Salvage', 
                                 `chop beam NA NA` = 'Salvage',
                                 `ice hdctx minibeam rituximab` = 'Salvage',
                                 `mvbpcns NA NA NA` = 'CNS', 
                                 `chop cop NA rituximab` = 'R-CHOP',
                                 `ice hdctx NA rituximab` = 'Salvage', 
                                 `mopp NA NA rituximab` = 'Other', 
                                 `peroral NA NA rituximab` = 'Palliative', 
                                 `cns_matrix itbehandling ice rituximab`  = 'CNS')) %>% 
  mutate(treatment_type_2nd_line = ifelse(treatment_type_2nd_line == 'Salvage' & high_dosis_treatment_2nd_line == 1,  
                                          'Salvage with HDT', treatment_type_2nd_line),
         treatment_type_2nd_line = ifelse(treatment_type_2nd_line == 'Salvage' & high_dosis_treatment_2nd_line == 0,  
                                          'Salvage without HDT', treatment_type_2nd_line),
         treatment_type_2nd_line = ifelse(treatment_type_2nd_line == 'Salvage' & is.na(high_dosis_treatment_2nd_line),  
                                          'Salvage without HDT', treatment_type_2nd_line)) %>% 
  mutate(treatment_type_2nd_line = ifelse(treatment_type_2nd_line == 'NA', NA, treatment_type_2nd_line),
         treatment_type_2nd_line = recode(treatment_type_2nd_line,
                                         `R-BAC` = 'BR',
                                          `R-CHOP` = 'CHOP',
                                         `HD-MTX` = 'Other',
                                          Palliative = 'Palliative_other'), 
         treatment_type_2nd_line_factor = factor(treatment_type_2nd_line, 
                                                 c('Salvage with HDT', 'Salvage without HDT',
                                                   'CHOP', 'R-CVP','BR', 'Palliative_other', 'CNS', 'Other')))

#High dose therapy
cohort4 <- cohort4 %>% 
  inner_join(RKKP_LYFO %>% select(patientid, Rec_Hoejdosisbehandling),
             by= "patientid") %>% 
  distinct(patientid, .keep_all = TRUE)

n_distinct(cohort4$patientid)
sum(is.na(cohort4$antiCD20_2L))

NA_antiCD20<- cohort4 %>% 
   filter(is.na(antiCD20_2L)) #240 is either "Yes" or "NO"- checked and they are all "No"

cohort4<- cohort4 %>% 
   mutate(antiCD20_2L=ifelse(is.na(antiCD20_2L), "No", antiCD20_2L))

colSums(is.na(cohort4)) #S0 79 patients missing treatment_type_2nd_line

NA_patients <- cohort4 %>% 
  filter(is.na(treatment_type_2nd_line))

#Ankat_beh tabel####
NA_patients1 <-NA_patients  %>%  pull(patientid)
secondeL <- open_dataset('/ngc/projects2/dalyca_r/andkat_r/medicine data/data/output tables/table_treatment_clean_2025_10_16.parquet') %>% 
  filter(patientid %in% NA_patients1) %>% 
  collect()
n_distinct(secondeL$patientid) #75

#+-90 dage from 2L beh. date
NA_pt <- NA_patients %>% 
  left_join(secondeL %>% select(patientid, tx_text, date_start, source, code_icd10), by = "patientid") %>% 
  mutate(date_treatment_2nd_line = as.Date(date_treatment_2nd_line),
         date_start = as.Date(date_start)) %>% 
  filter(date_start > date_treatment_2nd_line - 90,
         date_start < date_treatment_2nd_line + 90) %>% 
  group_by(patientid, tx_text) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(patientid) %>% 
  summarise(source=first(source),
            ny_text= paste(tx_text, collapse =  ":")) %>% 
  ungroup()

NA_pt1 <- NA_pt %>% 
  left_join(NA_patients %>% select(patientid, subtype,CD20), by = "patientid") 

#Categorising so it matched with the rest
NA_pt1_clean <- NA_pt1 %>% 
  mutate(treatment_type_2nd_line_ankat= case_when(
    str_detect(ny_text, "^cd20 antibody")~ 'Targeted-anti_CD20',
    str_detect(ny_text, "^cd20 antibody:chop")~ 'R-CHOP',
    str_detect(ny_text, "^cd20 antibody:rituximab")~ 'Targeted-anti_CD20',
    str_detect(ny_text, "^cd20 antibody:chop:gemcitabine:oxaliplatin")~ 'R-CHOP', 
    str_detect(ny_text, "^asct:asct (conditioning):cd20 antibody:complex cytostatic treatment")~ 'ASCT',
    str_detect(ny_text, "^basic cytostatic treatment")~ 'Other',
    str_detect(ny_text, "^basic cytostatic treatment:bortezomib:complex cytostatic treatment
               :cop:cyclophosphamide:doxorubicin:rituximab")~ 'Other', 
    str_detect(ny_text, "^basic cytostatic treatment:cd20 antibody:cisplatin|gemcitabine")~ 'Salvage', 
    str_detect(ny_text, "^basic cytostatic treatment:cd20 antibody:complex cytostatic treatment:gdp:rituximabe")~ 'Other', 
    str_detect(ny_text, "^cd20 antibody:chop:complex cytostatic treatment:ice:rituximab")~ 'R-CHOP',
    str_detect(ny_text, "^cd20 antibody:chop:gemcitabine:oxaliplatin")~ 'R-CHOP',
    str_detect(ny_text, "^cd20 antibody:complex cytostatic treatment")~ 'Other', 
    str_detect(ny_text, "^cd20 antibody:complex cytostatic treatment:dhap:pixantrone:rituximab")~ 'Other',
    str_detect(ny_text, "^cd20 antibody:complex cytostatic treatment:ice:lenalidomide:rituximab:vinblastine")~ 'R-CHOP', 
    str_detect(ny_text, "^cd20 antibody:methotrexate:mpv|rituximab:procarbazine:rituximab:vincristine")~ 'Other', 
    str_detect(ny_text, "^cd20 antibody:venetoclax")~ 'Targeted-anti_CD20',
    TRUE ~ ny_text)) 

cohort5<- cohort4 %>% 
  left_join(NA_pt1_clean %>% select(patientid, treatment_type_2nd_line_ankat),
            by="patientid") %>% 
  select(patientid, subtype, date_treatment_2nd_line, 
         CD20, high_dosis_treatment_2nd_line, 
         date_RT_2nd_line,antiCD20_2L, treatment_type_2nd_line_factor,Rec_Hoejdosisbehandling,
         treatment_type_2nd_line, treatment_type_2nd_line_ankat)%>% 
  mutate(treatment_type_2nd_line_ankat = ifelse(treatment_type_2nd_line_ankat == 'Salvage' & high_dosis_treatment_2nd_line == 1,  
                                          'Salvage with HDT', treatment_type_2nd_line_ankat),
         treatment_type_2nd_line_ankat = ifelse(treatment_type_2nd_line_ankat == 'Salvage' & high_dosis_treatment_2nd_line == 0,  
                                          'Salvage without HDT', treatment_type_2nd_line_ankat),
         treatment_type_2nd_line_ankat = ifelse(treatment_type_2nd_line_ankat == 'Salvage' & is.na(high_dosis_treatment_2nd_line),  
                                          'Salvage without HDT', treatment_type_2nd_line_ankat))

colSums(is.na(cohort5))

#Those who only recived RT
RT_only <- cohort5 %>% 
  filter(is.na(treatment_type_2nd_line),
         is.na(treatment_type_2nd_line_ankat),
         !is.na(date_RT_2nd_line)) %>% 
  select(patientid, subtype, date_treatment_2nd_line, 
         CD20, high_dosis_treatment_2nd_line, 
         date_RT_2nd_line,antiCD20_2L) %>% 
  inner_join(LYFO_clean %>% select(patientid, hospital_id),
             by= "patientid")
#Saving for maybe later 
#fwrite(RT_only, file = "Only_RT_treatment.txt", sep="\t") 

#Making a clean dataframe
secondL_treatment_df<- cohort5 %>% 
  mutate(all_treatment_2L= case_when(
    !is.na(treatment_type_2nd_line_ankat)~treatment_type_2nd_line_ankat,
    !is.na(treatment_type_2nd_line)~treatment_type_2nd_line,
    !is.na(treatment_type_2nd_line_factor)~treatment_type_2nd_line_factor,
    !is.na(date_RT_2nd_line)~"Radiotherapy 2L", 
    TRUE ~ NA_character_)) %>% 
  select(patientid, all_treatment_2L, subtype, CD20, 
         date_treatment_2nd_line, date_RT_2nd_line, 
         high_dosis_treatment_2nd_line,antiCD20_2L)


n_distinct(secondL_treatment_df$patientid)#1396
table(secondL_treatment_df$all_treatment_2L, useNA = "ifany") #56 only RT, 

#p-value ####
p_from_table <- function(tab, B= 1e6){
  chi <- suppressWarnings(chisq.test(tab))
  if (any(chi$expected <5)){
    if (prod(dim(tab))>4){
      f <- fisher.test(tab, simulate.p.value = TRUE, B= B)
      method <- "Fishers exact test (MC)"
      p <- f$p.value
    }  else {
      f <- fisher.test(tab)
      p<- f$p.value
      method <- "Fishers exact test"}}else{
      method <- "Chi-sqaure test"
      p <- chi$p.value}
    p_label <- if (p<0.0001)
      "<0.001"else
        sprintf("%.4f", p)
  list(
    p_raw = p,
    p_formatted = p_label,
    method = method)}

#Tabel3####
library(Publish)
df_table3<- secondL_treatment_df %>% 
  filter(subtype %in% c("DLBCL", "FL", "MCL"),
         CD20 %in% c("Positive", "Negative"),
         !is.na(all_treatment_2L),
         !all_treatment_2L %in% c("asct:asct (conditioning):cd20 antibody:complex cytostatic treatment", "Radiotherapy 2L")) %>% 
  mutate(subtype_cd20= factor(paste(subtype, CD20, sep="_"),
                              levels= c("DLBCL_Positive", "DLBCL_Negative",
                                        "FL_Positive", "FL_Negative",
                                        "MCL_Positive", "MCL_Negative"),
                              labels = c("DLBCL CD20 Positive", "DLBCL CD20 Negative",
                                         "FL CD20 Positive", "FL CD20 Negative",
                                         "MCL CD20 Positive", "MCL CD20 Negative"))) %>% 
           mutate(all_treatment_2L = factor(all_treatment_2L,
                                            levels= c("Salvage with HDT",
                                                      "Salvage without HDT",
                                                      "CHOP",
                                                      "Palliative_other",
                                                      "R-CVP",
                                                      "BR",
                                                      "MCL",
                                                      "BsAbs",
                                                      "CNS",
                                                      "anti_CD20",
                                                      "Targeted-anti_CD20",
                                                      "Targeted",
                                                      "Other")))
         
table3 <- utable(subtype_cd20~all_treatment_2L, data=df_table3)

colSums(is.na(df_table3))

#write_utable(table3, table_n= "3") 
table3

table3_CD20total <- utable(CD20~all_treatment_2L, data=df_table3)
table3_CD20total


#Ekstra supp_tabel 4####
df_tableS4<-df_table3 %>% 
  mutate(CD20= factor(CD20, levels = c("Positive", "Negative"),
                      labels = c("CD20 Positive", "CD20 Negative")),
         antiCD20_2L=factor(antiCD20_2L, levels = c("Yes", "No"),
                            labels = c("Anti-CD20 therapy Yes", "Anti-CD20 therapy No")),
         subtype_CD20_group_1=interaction(CD20, antiCD20_2L, sep="|"))
Supp_table4 <- utable( subtype_CD20_group_1~all_treatment_2L, df_tableS4 %>% filter(!is.na(antiCD20_2L)) )

Supp_table4
getwd()
#write_utable(Supp_table4, table_n= "5")
n_distinct(secondL_treatment_df$patientid)


#p value for anti CD20 YES/NO
anti_CD20_tab <-table(df_tableS4$CD20,df_tableS4$antiCD20_2L)
anti_CD20_tab

anti_cd20_p_val<- p_from_table(anti_CD20_tab)
anti_cd20_p_val$p_raw 
anti_cd20_p_val$p_formatted 
anti_cd20_p_val$method 



