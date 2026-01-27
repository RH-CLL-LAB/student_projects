# title: "CD20-Master script"
#Project: CD20 expression in relpased B-cell lymphomas
# output: html_document
# Start date: "2025-04-02"
# Author: Sanaz M. Gholy
#---
#This master script runds the full CD20 projects pipeline.
#It sources all scripts in the correct order to:
#Clean and prepare data
#Define study population
#Generate baseline data
#Perform survival analyses
#Produce tables and figures


#_________________Textmining and cohort definition_____________________
source('/ngc/projects2/dalyca_r/sangho_r/scripts/Working_with_finaldata/github/1_CD20_textmining_SDSpato.R')
source('/ngc/projects2/dalyca_r/sangho_r/scripts/Working_with_finaldata/github/2_defining_studypopulation_CONSORT_CD20.R')
source('/ngc/projects2/dalyca_r/sangho_r/scripts/Working_with_finaldata/github/3_Cleaning_cohort.R')
#_________________Baseline charachteristics______________________________
source('/ngc/projects2/dalyca_r/sangho_r/scripts/Working_with_finaldata/github/4_1st_line_categorization.R')
source('/ngc/projects2/dalyca_r/sangho_r/scripts/Working_with_finaldata/github/5_ASCT_polyrx.R')
source('/ngc/projects2/dalyca_r/sangho_r/scripts/Working_with_finaldata/github/6_CD20_cohort_baseline_charac.R')
#_______________Tables, figures and survival analyses____________________
source('/ngc/projects2/dalyca_r/sangho_r/scripts/Working_with_finaldata/github/7_Table1_Table2.R')
source('/ngc/projects2/dalyca_r/sangho_r/scripts/Working_with_finaldata/github/8_all_survivalanalysis_and_sensitivity.R')
source('/ngc/projects2/dalyca_r/sangho_r/scripts/Working_with_finaldata/github/9_2L_categorization_table3.R')
