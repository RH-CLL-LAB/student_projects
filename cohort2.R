## This script analyzes survival in 1L treated SLL/CLL
# christian brieghel/emma siig
# date: 2025-02-24

cohort <- read_csv2('/ngc/projects2/dalyca_r/emmsii_r/sandbox/tables/cohort.csv')
cohort3 = cohort %>% 
  mutate(age_cut = cut(age, c(0, 50, 75, Inf)),
         time_treatment_to_death = diff_days(date_treatment_1st_line, date_death_fu)) %>% 
  filter(time_treatment_to_death >= 0)

cohort3 %>% nrow_npatients()
cohort3 %>% names() 
cohort3$age_cut %>% table

plot_os = KM_plot(survfit(Surv(time_treatment_to_death/365.25, status) ~ age_cut, cohort3),
                  title = 'Age',
                  labs = c('<50', '50-75', '>75 years'),
                  xlab = 'Time (years from 1L treatment)',
                  xlim = c(0,10))
plot_os # view plot

ggsave('/ngc/projects2/dalyca_r/emmsii_r/sandbox/figures/figure_1.png',
       arrange_ggsurvplots(list(plot_os), nrow = 1, ncol = 1),
       height = 12,
       width = 10,
       units = 'cm',
       dpi = 300) # save ggsurv-plot

pairwise_survdiff(Surv(time_treatment_to_death/365.25, status) ~ age_cut, cohort3)
plot_pairwiseLR = pairwise_survdiff(Surv(time_treatment_to_death/365.25, status) ~ age_cut, cohort3) %>% 
  tile_pairwise_survdiff(position = 'LL',
                         palette = c(1:3)) # pairwise log-rank

ggsave('/ngc/projects2/dalyca_r/emmsii_r/sandbox/figures/figure_1_labs.png',
       plot_pairwiseLR,
       height = 5,
       width = 6,
       units = 'cm',
       dpi = 300) # save regular ggplot

ggsurvplot(survfit(Surv(time_treatment_to_death/365.25, status) ~ sex, cohort3)) # KM_plot same as ggsurvplot
coxph(Surv(time_treatment_to_death/365.25, status) ~ sex, cohort3) %>% publish() #cox regression
