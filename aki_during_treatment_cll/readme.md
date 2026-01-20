# author: Emma Siig
#

uploaded to GitHub January 20th, 2026

Risk of Acute Kidney Injury and Cytopenias During Treatment of Chronic Lymphocytic Leukemia � data and scripts

1cohort.r
This script identifies patients with CLL or SLL and adds some characteristics
CLL_clean and SLL_clean (and LAB_IGHVIMGT) are used to get data on: patientid, date_diagnosis, hospital_id, date_treatment_1st_line, treatment_type, age, binetstage, B2M, IGHV, TP53_ab
Cleaning of treatment types are done. 

2Baseline_characteristics.R
This script adds (and cleans) baseline characteristics to the CLL and SLL cohort
Smoking and drinking history: CLL_TREAT and SP_Social_Hx
BMI: SP_VitaleVaerdier
Comorbidities (CCI score): diagnoses_all
FISH: CLL_clean
IGHV: CLL_TREAT

3AdverseEvents.R
This script finds biochemical adverse events and adds them to the CLL and SLL cohort. We also define end of follow up, start and more. 
Uses dataset: laboratorymeasurements
Anemia: AE_anemia 
Neutropenia: AE_neutrophil_count_decreased 
Lymphopenia: AE_lymphocyte_count_decreased 
Creatinine increased (proxy of acute kidney injure): AE_creatinine_increased 
Thrombocytopenia: AE_platelets_descreased 
Leukopenia: AE_WBC_decreased

4DescriptiveStatistics.R
This script makes descriptive statistics on the CLL and SLL cohort
Looks into background information, adverse event by treatment type, age and sex. And more. 

5Statistics.R
This script looks at statistics for the CLL and SLL cohort. 
-	Kaplan-Meier estimated cumulative risk of adverse events
-	Cumulative incidence curves (Aalen-Johansen) 
-	Cox proportional hazard models (univariable) 
-	Cox proportional hazard models (multivariable) 
-	Forest plots 

