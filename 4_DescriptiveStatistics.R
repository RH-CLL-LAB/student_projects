# This script looks at characteristics for the CLL and SLL cohort
# Date: 18.03.2025
# Author: Emma Siig 
# Comment: Run 1cohort.R first, then 2Baseline_characteristics.R, then 3.AdverseEvents.R. Or downloade table: cohort_riskfactors_AE1

cohort_riskfactors_AE1 <- readRDS("/ngc/projects2/dalyca_r/emmsii_r/sandbox/tables/3cohort_riskfactors_AE1.rds")

cohort_done <- cohort_riskfactors_AE1

#### Basic descriptive statistics ####
hist(CLL_clean$age)
str(CLL_clean)
summary(CLL_clean)
summary(cohort)
summary(RKKP_CLL)
summary(cohort_2)
summary(cohort3)
str(cohort)
nrow(cohort)
hist(cohort_done$BMI)

utable(~treatment_type, cohort_done)

utable(~age+sex+date_treatment_1st_line+type+B2M+IGHV+TP53_ab, cohort)
utable(~Q(age), cohort)
utable(~age+sex+date_treatment_1st_line+type+stage+B2M+IGHV+TP53_ab+BMI+del17p_at_treatment+FISH_status_all, cohort_done)

#### Visuals ####
plot(cohort_test$age, cohort_test$time_to_treatment_days)
ggplot(cohort_test, aes(x = age, y = time_to_treatment_days)) +
  geom_point() +
  labs(title = "Age x time to treatment", 
       x = "Age",
       y = "Time to treatment")
summary(cohort_test)

ggplot(cohort_test, aes(x = date_diagnosis, y = date_treatment_1st_line)) +
  geom_point() +
  labs(title = "Diagnosis x first treatment (dates)",
       x = "Diagnosis",
       y = "First treatment")
#OBS: Some had their first treatment before their diagnosis

pairs(cohort_test %>% select_if(is.numeric)) #Scattermatrix

#### Prognostic factors ####
ggplot(cohort_test, aes(x = stage)) + 
  geom_bar(fill = "skyblue") + 
  labs(title = "Antal patienter per stadie", x = "Stadie", y = "Antal patienter") + 
  theme_minimal()

table(cohort$stage)

table(cohort$type)  #Number of patients (check)
missing_SLL <- cohort_test %>%
  filter(type == "SLL" & is.na(stage))  #Looks into missing data
print(missing_SLL)

table(cohort_test$IGHV)
table(cohort$TP53_ab)
table(cohort$B2M)

utable(~age+B2M+IGHV+TP53_ab+BMI+smoking_history+del17p_at_treatment, cohort_riskfactors_AE1) #Pronogstic factors

hist(cohort_riskfactors_AE1$age)
hist(cohort_riskfactors_AE1$status)


#### visuals of background information ####
ggplot(cohort_riskfactors_AE1, aes(x = year(date_birth))) +           #Birth date
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  labs(x = "Birth date", y = "Patients") +
  theme_minimal()


ggplot(cohort_riskfactors_AE1, aes(x = year(date_diagnosis))) +       #Date of diagnosis
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  labs(x = "Date of diagnosis", y = "Patients") +
  theme_minimal()


ggplot(cohort_done, aes(x = year(date_treatment_1st_line))) + #Date of first treatment
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  labs(x = "Date of first treatment", y = "Patients") +
  theme_minimal()



ggplot(cohort_done, aes(x = year(date_treatment_1st_line))) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  stat_bin(
    binwidth = 1,
    geom = "text",
    aes(label = ..count..),
    vjust = -0.5,        # Justerer placering over søjler
    color = "black",
    size = 3.5
  ) +
  scale_x_continuous(breaks = seq(2005, 2025, by = 2)) +  # tilpas årstal efter dine data
  labs(
    x = "Date of first treatment",
    y = "Patients"
  ) +
  theme_minimal()



#### Looking into missing data (adverse events) ####
cohort_riskfactors_AE1_missing <- cohort_riskfactors_AE1 %>%
  left_join(
    CLL_clean %>% select(patientid, hospital_id),
    by = "patientid"
  )

cohort_riskfactors_AE1_missing %>%
  filter(is.na(ae_anemia)) %>%
  mutate(treatment_year = year(date_treatment_1st_line)) %>%
  count(treatment_year) %>%
  ggplot(aes(x = treatment_year, y = n)) +
  geom_col(fill = "steelblue") +
  labs(
    title = "Antal patienter uden Hgb-måling pr. behandlingsår",
    x = "Behandlingsår",
    y = "Antal patienter uden Hgb-måling"
  ) +
  theme_minimal()

#### ggplot of prognostic factors ####
ggplot(cohort_riskfactors_AE1, aes(x = fct_rev(fct_infreq(treatment_type)))) + #treatment type
    geom_bar(fill = "dark green") +
    labs(x = "Treatment type", y = "Patients") +
    theme_minimal()

utable(~smoking_history, cohort_riskfactors_AE1)


ggplot(cohort_riskfactors_AE1, aes(x = stage, fill = type)) + #stage
  geom_bar(position = "dodge") +
  labs(x = "Stage", y = "Patients", fill = "Diagnosis") +
  theme_minimal()

ggplot(cohort_riskfactors_AE1, aes(x = B2M)) + #Beta 2 microglobuline level at diagnosis
  geom_bar() +
  labs(x = "B2M", y = "Patients") +
  theme_minimal()

ggplot(data = filter(cohort_riskfactors_AE1, type == "CLL"), aes(x = IGHV)) + #IGHV mutational status (CLL)
  geom_bar(fill = "steelblue") +
  labs(x = "IGHV mutational status", y = "Patients") +
  theme_minimal()


ggplot(data = filter(cohort_riskfactors_AE1, type == "CLL"), aes(x = TP53_ab)) + #TP53 mutation at treatment (CLL)
  geom_bar(fill = "steelblue") +
  labs(x = "TP53 mutation at treatment", y = "Patients") +
  theme_minimal()

ggplot(cohort_riskfactors_AE1, aes(x = smoking_history)) + #Smoking history with NA's
  geom_bar() +
  labs(x = "Smoking history", y = "Patients") +
  theme_minimal()

cohort_riskfactors_AE1 %>%              #Smoking history without NA's
  filter(!is.na(smoking_history)) %>%
  ggplot(aes(x = smoking_history)) +
  geom_bar(fill = "steelblue") +
  labs(x = "Smoking history", y = "Patients", title = "Histogram of smoking history") +
  theme_minimal()

cohort_riskfactors_AE1 %>%              #Drinking history without NA's
  filter(!is.na(drinking_history)) %>%
  ggplot(aes(x = drinking_history)) +
  geom_bar(fill = "steelblue") +
  labs(x = "Drinking history", y = "Patients", title = "Histogram of drinking history") +
  theme_minimal()

ggplot(cohort_riskfactors_AE1, aes(x = del17p_at_treatment)) + #del17p_at_treatment
  geom_bar() +
  labs(x = "del17p_at_treatment", y = "Patients") +
  theme_minimal()

ggplot(cohort_done, aes(x = FISH_status_all)) + #FISH status
  geom_bar() +
  labs(x = "FISH mutational status", y = "Patients") +
  theme_minimal()


cohort_riskfactors_AE1 %>%              #FISH mutational status without NA's
  filter(!is.na(FISH_status)) %>%
  ggplot(aes(x = FISH_status)) +
  geom_bar(fill = "steelblue") +
  labs(x = "FISH mutational status", y = "Patients", title = "Histogram of FISH mutational status") +
  theme_minimal()


#### ggplots of adverse events (grades) ####
utable(~age+ae_anemia+ae_neutrophil_count_decreased+ae_lymphocyte_count_decreased+ae_creatinine_increased+ae_platelets_decreased+ae_WBC_decreased, cohort_done) #Adverse events 
table(cohort_done$anemia_binary)
table(cohort_done$leukopenia_binary)
628  /(1426  +628   )*100


ggplot(cohort_done, aes(x = ae_anemia)) + #Anemia
  geom_bar(fill="darkblue") +
  labs(x = "Anemia Grade", y = "Patients") +
  theme_minimal()+
  theme(
    axis.title = element_text(size=18),
    axis.text = element_text(size = 20)
  )

ggplot(cohort_done, aes(x = ae_neutrophil_count_decreased)) + #Neutropenia
  geom_bar(fill="darkblue") +
  labs(x = "Neutropenia Grade", y = "Patients") +
  theme_minimal()+
  theme(
    axis.title = element_text(size=18),
    axis.text = element_text(size = 20)
  )

ggplot(cohort_done, aes(x = ae_lymphocyte_count_decreased)) + #Lymphopenia
  geom_bar(fill="darkblue") +
  labs(x = "Lymphopenia Grade", y = "Patients") +
  theme_minimal()+
  theme(
    axis.title = element_text(size=18),
    axis.text = element_text(size = 20)
  )

ggplot(cohort_done, aes(x = ae_platelets_decreased)) + #Thrombocytopenia
  geom_bar(fill="darkblue") +
  labs(x = "Thrombocytopenia Grade", y = "Patients") +
  theme_minimal()+
  theme(
    axis.title = element_text(size=18),
    axis.text = element_text(size = 20)
  )

ggplot(cohort_done, aes(x = ae_WBC_decreased)) + #Leukopenia
  geom_bar(fill="darkblue") +
  labs(x = "Leukopenia Grade", y = "Patients") +
  theme_minimal()+
  theme(
    axis.title = element_text(size=18),
    axis.text = element_text(size = 20)
  )

ggplot(cohort_done, aes(x = ae_creatinine_increased)) + #Creatinine
  geom_bar(fill="darkblue") +
  labs(x = "Creatinine Increase Grade", y = "Patients") +
  theme_minimal()+
  theme(
    axis.title = element_text(size=18),
    axis.text = element_text(size = 20)
)

#### Adverse events by treatment type ####

#Creatinine (acute kidney injury)
creatinine_by_treatment_type <- cohort_riskfactors_AE1 %>%    #Making a dataset for my plot. Looks at percentage with AKI by each treatment type
  group_by(treatment_type) %>%
  summarise(
    n_total_creatinine = n(),
    n_creatinine = sum(creatinine_binary == 1, na.rm = TRUE),
    pct_creatinine = round(100 * n_creatinine / n_total_creatinine, 1)
  ) %>%
  arrange(desc(pct_creatinine))  #Sorts by highest %

creatinine_by_treatment_type

ggplot(creatinine_by_treatment_type, aes(x = reorder(treatment_type, -pct_creatinine), y = pct_creatinine)) + #Making a bar chart
  geom_col(fill = "steelblue") +
  geom_text(aes(label = paste0(n_creatinine, "/", n_total_creatinine)),                        #Adds the number of patients to the chart
            vjust = -0.5, size = 3.5) +                                                          #Adjusting the text size
  labs(
    title = "Percentage with creatinine increase (AKI) by treatment type",
    x = "Treatment Type",
    y = "Percentage with creatinine increase"
  ) +
  theme_minimal()

#Anemia
anemia_by_treatment_type <- cohort_riskfactors_AE1 %>%  #Making a dataset for my plot. Looks at percentage with anemia by each treatment type
  group_by(treatment_type) %>%
  summarise(
    n_total_anemia = n(),
    n_anemia = sum(anemia_binary == 1, na.rm = TRUE),
    pct_anemia = round(100 * n_anemia / n_total_anemia, 1)
  ) %>%
  arrange(desc(pct_anemia))                                      #Sorts by highest %

anemia_by_treatment_type                                         #Looks at the table

ggplot(anemia_by_treatment_type, aes(x = reorder(treatment_type, -pct_anemia), y = pct_anemia)) +   #Making a bar chart
  geom_col(fill = "steelblue") +
  geom_text(aes(label = paste0(n_anemia, "/", n_total_anemia)),                               #Adds the number of patients to the chart
            vjust = -0.5, size = 3.5) +                                                       #Adjusting the text size
  labs(
    title = "Percentage with anemia by treatment type",
    x = "Treatment Type",
    y = "Percentage with Anemia"
  ) +
  theme_minimal()

#Neutropenia
neutropenia_by_treatment_type <- cohort_riskfactors_AE1 %>%    #Making a dataset for my plot. Looks at percentage with Neutropenia by each treatment type
  group_by(treatment_type) %>%
  summarise(
    n_total_neutropenia = n(),
    n_neutropenia = sum(neutropenia_binary == 1, na.rm = TRUE),
    pct_neutropenia = round(100 * n_neutropenia / n_total_neutropenia, 1)
  ) %>%
  arrange(desc(pct_neutropenia))  #Sorts by highest %

neutropenia_by_treatment_type

ggplot(neutropenia_by_treatment_type, aes(x = reorder(treatment_type, -pct_neutropenia), y = pct_neutropenia)) + #Making a bar chart
  geom_col(fill = "steelblue") +
  geom_text(aes(label = paste0(n_neutropenia, "/", n_total_neutropenia)),                        #Adds the number of patients to the chart
            vjust = -0.5, size = 3.5) +                                                          #Adjusting the text size
  labs(
    title = "Percentage with neutropenia by treatment type",
    x = "Treatment Type",
    y = "Percentage with Neutropenia"
  ) +
  theme_minimal()

#Lymphopenia 
lymphopenia_by_treatment_type <- cohort_riskfactors_AE1 %>%    #Making a dataset for my plot. Looks at percentage with lymphopenia by each treatment type
  group_by(treatment_type) %>%
  summarise(
    n_total_lymphopenia = n(),
    n_lymphopenia = sum(lymphopenia_binary == 1, na.rm = TRUE),
    pct_lymphopenia = round(100 * n_lymphopenia / n_total_lymphopenia, 1)
  ) %>%
  arrange(desc(pct_lymphopenia))  #Sorts by highest %

lymphopenia_by_treatment_type

ggplot(lymphopenia_by_treatment_type, aes(x = reorder(treatment_type, -pct_lymphopenia), y = pct_lymphopenia)) + #Making a bar chart
  geom_col(fill = "steelblue") +
  geom_text(aes(label = paste0(n_lymphopenia, "/", n_total_lymphopenia)),                        #Adds the number of patients to the chart
            vjust = -0.5, size = 3.5) +                                                          #Adjusting the text size
  labs(
    title = "Percentage with lymphopenia by treatment type",
    x = "Treatment Type",
    y = "Percentage with lymphopenia"
  ) +
  theme_minimal()

#Thrombocytopenia 
thrombocytopenia_by_treatment_type <- cohort_riskfactors_AE1 %>%    #Making a dataset for my plot. Looks at percentage with thrombocytopenia by each treatment type
  group_by(treatment_type) %>%
  summarise(
    n_total_thrombocytopenia = n(),
    n_thrombocytopenia = sum(thrombocytopenia_binary == 1, na.rm = TRUE),
    pct_thrombocytopenia = round(100 * n_thrombocytopenia / n_total_thrombocytopenia, 1)
  ) %>%
  arrange(desc(pct_thrombocytopenia))  #Sorts by highest %

thrombocytopenia_by_treatment_type

ggplot(thrombocytopenia_by_treatment_type, aes(x = reorder(treatment_type, -pct_thrombocytopenia), y = pct_thrombocytopenia)) + #Making a bar chart
  geom_col(fill = "steelblue") +
  geom_text(aes(label = paste0(n_thrombocytopenia, "/", n_total_thrombocytopenia)),                        #Adds the number of patients to the chart
            vjust = -0.5, size = 3.5) +                                                          #Adjusting the text size
  labs(
    title = "Percentage with thrombocytopenia by treatment type",
    x = "Treatment Type",
    y = "Percentage with thrombocytopenia"
  ) +
  theme_minimal()

#Leukopenia 
leukopenia_by_treatment_type <- cohort_riskfactors_AE1 %>%    #Making a dataset for my plot. Looks at percentage with leukopenia by each treatment type
  group_by(treatment_type) %>%
  summarise(
    n_total_leukopenia = n(),
    n_leukopenia = sum(leukopenia_binary == 1, na.rm = TRUE),
    pct_leukopenia = round(100 * n_leukopenia / n_total_leukopenia, 1)
  ) %>%
  arrange(desc(pct_leukopenia))  #Sorts by highest %

leukopenia_by_treatment_type

ggplot(leukopenia_by_treatment_type, aes(x = reorder(treatment_type, -pct_leukopenia), y = pct_leukopenia)) + #Making a bar chart
  geom_col(fill = "steelblue") +
  geom_text(aes(label = paste0(n_leukopenia, "/", n_total_leukopenia)),                        #Adds the number of patients to the chart
            vjust = -0.5, size = 3.5) +                                                          #Adjusting the text size
  labs(
    title = "Percentage with leukopenia by treatment type",
    x = "Treatment Type",
    y = "Percentage with leukopenia"
  ) +
  theme_minimal()

#### Adverse events by age at treatment and sex ####

cohort_done <- cohort_done %>%        #Making age groups using age at treatment
  mutate(age_group = case_when(
    age_at_treatment < 60 ~ "<60",
    age_at_treatment >= 60 & age_at_treatment < 70 ~ "60-69",
    age_at_treatment >= 70 & age_at_treatment < 80 ~ "70-79",
    age_at_treatment >= 80 ~ "80+"
  ))


utable(~age_group, cohort_done)

#Anemia
table_age_sex__anemia <- cohort_done %>%
        group_by(sex, age_group) %>%
        summarise(n = n(),
                  anemia_cases = sum(anemia_binary == 1, na.rm = TRUE),
                  anemia_prop = anemia_cases / n) %>%
        ungroup()
      
print(table_age_sex__anemia)

ggplot(table_age_sex__anemia, aes(x = age_group, y = anemia_prop, fill = sex)) +
  geom_col(position = "dodge") +
  labs(
    title = "Proportion of patients with anemia by age group and sex",
    x = "Age group", y = "Proportion with anemia"
  ) +
  scale_fill_manual(
    values = c("M" = "forestgreen", "F" = "pink")
  ) +
  theme_minimal()


#Neutropenia
table_age_sex__neutropenia <- cohort_done %>%
  group_by(sex, age_group) %>%
  summarise(n = n(),
            neutropenia_cases = sum(neutropenia_binary == 1, na.rm = TRUE),
            neutropenia_prop = neutropenia_cases / n) %>%
  ungroup()

print(table_age_sex__neutropenia)

ggplot(table_age_sex__neutropenia, aes(x = age_group, y = neutropenia_prop, fill = sex)) +
  geom_col(position = "dodge") +
  labs(
    title = "Proportion of patients with neutropenia by age-group and sex",
    x = "Age-group", y = "Proportion with neutropenia"
  ) +
  scale_fill_manual(
    values = c("M" = "forestgreen", "F" = "pink")
  ) +
  theme_minimal()

#Lymphopenia
table_age_sex__lymphopenia <- cohort_done %>%
  group_by(sex, age_group) %>%
  summarise(n = n(),
            lymphopenia_cases = sum(lymphopenia_binary == 1, na.rm = TRUE),
            lymphopenia_prop = lymphopenia_cases / n) %>%
  ungroup()

print(table_age_sex__lymphopenia)

ggplot(table_age_sex__lymphopenia, aes(x = age_group, y = lymphopenia_prop, fill = sex)) +
  geom_col(position = "dodge") +
  labs(
    title = "Proportion of patients with lymphopenia by age-group and sex",
    x = "Age-group", y = "Proportion with lymphopenia"
  ) +
  scale_fill_manual(
    values = c("M" = "forestgreen", "F" = "pink")
  ) +
  theme_minimal()

#Thrombocytopenia
table_age_sex__thrombocytopenia <- cohort_done %>%
  group_by(sex, age_group) %>%
  summarise(n = n(),
            thrombocytopenia_cases = sum(thrombocytopenia_binary == 1, na.rm = TRUE),
            thrombocytopenia_prop = thrombocytopenia_cases / n) %>%
  ungroup()

print(table_age_sex__thrombocytopenia)

ggplot(table_age_sex__thrombocytopenia, aes(x = age_group, y = thrombocytopenia_prop, fill = sex)) +
  geom_col(position = "dodge") +
  labs(
    title = "Proportion of patients with thrombocytopenia by age-group and sex",
    x = "Age-group", y = "Proportion with thrombocytopenia"
  ) +
  scale_fill_manual(
    values = c("M" = "forestgreen", "F" = "pink")
  ) +
  theme_minimal()

#Leukopenia
table_age_sex__leukopenia <- cohort_done %>%
  group_by(sex, age_group) %>%
  summarise(n = n(),
            leukopenia_cases = sum(leukopenia_binary == 1, na.rm = TRUE),
            leukopenia_prop = leukopenia_cases / n) %>%
  ungroup()

print(table_age_sex__leukopenia)

ggplot(table_age_sex__leukopenia, aes(x = age_group, y = leukopenia_prop, fill = sex)) +
  geom_col(position = "dodge") +
  labs(
    title = "Proportion of patients with leukopenia by age-group and sex",
    x = "Age-group", y = "Proportion with leukopenia"
  ) +
  scale_fill_manual(
    values = c("M" = "forestgreen", "F" = "pink")
  ) +
  theme_minimal()

#Creatinine increase (AKI)
table_age_sex__creatinine <- cohort_done %>%
  group_by(sex, age_group) %>%
  summarise(n = n(),
            creatinine_cases = sum(creatinine_binary == 1, na.rm = TRUE),
            creatinine_prop = creatinine_cases / n) %>%
  ungroup()

print(table_age_sex__creatinine)

ggplot(table_age_sex__creatinine, aes(x = age_group, y = creatinine_prop, fill = sex)) +
  geom_col(position = "dodge") +
  labs(
    title = "Proportion of patients with creatinine increase (AKI) by age-group and sex",
    x = "Age-group", y = "Proportion with creatinine increase (AKI)"
  ) +
  scale_fill_manual(
    values = c("M" = "forestgreen", "F" = "pink")
  ) +
  theme_minimal()

#### Looks into asymmetry regarding time to treatment ####
summary(cohort_test$time_to_treatment_days)  # Looks into time to treatment
hist(cohort_test$time_to_treatment_days)  # Visualization 

plot(cohort_test$age, cohort_test$time_to_treatment_days, 
     xlab = "Alder ved diagnose", ylab = "Tid til behandling (dage)", 
     main = "Alder vs. tid til behandling", col = "blue", pch = 16)

hist(cohort_test$time_to_treatment_days, breaks = 30, main = "Tid fra diagnose til behandling", xlab = "Dage", col = "lightblue")

boxplot(cohort_test$time_to_treatment_days, main = "Tid fra diagnose til behandling", ylab = "Dage", col = "lightblue")



#### Changing variables from character til numeric variables ####
class(cohort_test$sex) #datatypen for sex er character
cohort_test <- cohort_test %>%
  mutate(sex_numeric = ifelse(sex == "M", 1, 0))
table(cohort_test$sex, cohort_test$sex_numeric) #M = 1 og K = 0


cohort_test %>% group_by(stage) %>% summarise(mean_age = mean(age))
cohort_test %>% group_by(time_to_treatment) %>% summarise (mean_age = mean(age))
cohort_test %>% group_by(status) %>% summarise(mean_age = mean(age)) #mean age for dead patients
cohort_test %>% group_by(sex) %>% summarise (mean_age = mean(age))

cohort_test %>%
  count(status)

ggplot(cohort_test, aes(x = age, fill = sex)) + 
  geom_histogram(binwidth = 5, position = "dodge", alpha = 0.7) +  
  scale_fill_manual(values = c("M" = "blue", "K" = "red")) +  
  labs(title = "Fordeling af nye tilfælde efter alder og køn", 
       x = "Alder", y = "Antal nye tilfælde") +  
  theme_minimal()

#### Visuals of treatment types ####
ggplot(cohort_done, aes(x = treatment_type)) +
  geom_bar(fill = "steelblue", color = "black") +
  labs(title = "Fordeling af patienter på forskellige behandlingstyper",
       x = "Behandlingstype",
       y = "Antal patienter") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  #Rotating labels





#### Summary ####
summary(cohort_done)
table(cohort_done$sex)
utable(cohort_done$treatment_type)
unique(cohort_done$treatment_type) #83 different

#### All adverse events combined ####
cohort_done <- cohort_done %>%                 #Combining all adverse events
  mutate(any_ae = if_else(
    rowSums(across(ends_with("_binary")) == 1, na.rm = TRUE) > 0,
    1, 0
  ))

cohort_done %>%                               #Looking into the data
  count(any_ae) %>%
  mutate(percent = round(100 * n / sum(n), 1))

cohort_done <- cohort_done %>%                #Looking into number of adverse events
  mutate(n_ae = rowSums(across(ends_with("_binary")) == 1, na.rm = TRUE))

cohort_done %>%                               #Looking into number of adverse events
  count(n_ae) %>%
  mutate(percent = round(100 * n / sum(n), 1))

#By sex
cohort_done %>%
  group_by(sex) %>%
  count(any_ae) %>%
  mutate(percent = round(100 * n / sum(n), 1))

#By age
cohort_done %>%
  group_by(age_group) %>%
  count(any_ae) %>%
  mutate(percent = round(100 * n / sum(n), 1))

ae_by_age <- cohort_done %>%            #Beregn andel med AE per aldersgruppe
  group_by(age_group) %>%
  summarise(
    percent_with_ae = 100 * mean(any_ae, na.rm = TRUE)
  )

ggplot(ae_by_age, aes(x = age_group, y = percent_with_ae)) + 
  geom_col(fill = "dark green") +
  labs(
    title = "Proportion of Patients with ≥1 Adverse Event by Age Group",
    x = "Age Group",
    y = "Percent with Adverse Event (%)"
  ) +
  geom_text(aes(label = round(percent_with_ae, 1)), vjust = -0.5) +
  theme_minimal()

#By treatment type
cohort_done %>%
  group_by(treatment_type) %>%
  count(any_ae) %>%
  mutate(percent = round(100 * n / sum(n), 1))


ae_by_treatment <- cohort_done %>%            #Beregn andel med AE per treatment type
  group_by(treatment_type) %>%
  summarise(
    percent_with_ae_treatment = 100 * mean(any_ae, na.rm = TRUE)
  )

ae_by_treatment$treatment_type <- reorder(ae_by_treatment$treatment_type, -ae_by_treatment$percent_with_ae_treatment) #Sorting treatment_type by procent

ggplot(ae_by_treatment, aes(x = treatment_type, y = percent_with_ae_treatment)) + 
  geom_col(fill = "dark green") +
  labs(
    title = "Proportion of Patients with ≥1 Adverse Event by Treatment Type",
    x = "Treatment type",
    y = "Percent with Adverse Event (%)"
  ) +
  geom_text(aes(label = round(percent_with_ae_treatment, 1)), vjust = -0.5) +
  theme_minimal()

#Bendamustine
bendamustine_patients <- cohort_done %>%
  filter(treatment_type == "bendamustine")

bendamustine_patients %>% #Count
  summarise(
    Anemia = sum(anemia_binary == 1, na.rm = TRUE),
    Neutropenia = sum(neutropenia_binary == 1, na.rm = TRUE),
    Lymphopenia = sum(lymphopenia_binary == 1, na.rm = TRUE),
    Thrombocytopenia = sum(thrombocytopenia_binary == 1, na.rm = TRUE),
    Leukopenia = sum(leukopenia_binary == 1, na.rm = TRUE),
    Creatinine = sum(creatinine_binary == 1, na.rm = TRUE)
  )

bendamustine_patients %>% #%
  summarise(
    Anemia = mean(anemia_binary == 1, na.rm = TRUE) * 100,
    Neutropenia = mean(neutropenia_binary == 1, na.rm = TRUE) * 100,
    Lymphopenia = mean(lymphopenia_binary == 1, na.rm = TRUE) * 100,
    Thrombocytopenia = mean(thrombocytopenia_binary == 1, na.rm = TRUE) * 100,
    Leukopenia = mean(leukopenia_binary == 1, na.rm = TRUE) * 100,
    Creatinine = mean(creatinine_binary == 1, na.rm = TRUE) * 100
  )



#### Dead without adverse event, before follow-up ####
all_deaths_without_AE_in_followup <- bind_rows(dead_without_anemia_in_followup,
                         dead_without_neutropenia_in_followup,
                         dead_without_creatinine_in_followup,
                         dead_without_leukopenia_in_followup,
                         dead_without_lymphopenia_in_followup,
                         dead_without_thrombocytopenia_in_followup) 

nrow_npatients(all_deaths_without_AE_in_followup)



#### SAVE TABLES ####
saveRDS(cohort_done, "/ngc/projects2/dalyca_r/emmsii_r/sandbox/tables/4cohort_riskfactors_AE2.rds")
