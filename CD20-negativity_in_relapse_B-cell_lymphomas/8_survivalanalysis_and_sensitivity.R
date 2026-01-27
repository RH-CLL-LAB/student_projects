# title: "CD20-data_Survival analysis"
# output: html_document
# date: "2025-12-17"
# Author: Sanaz M. Gholy
#---

library(data.table)
library(tidyverse)
library(survival)
library(dplyr)
library(survminer)
library(patchwork)
library(Publish)
library(forcats)
library(ggsurvfit)

source('/ngc/projects2/dalyca_r/clean_r/load_dalycare_package.R')
load_dataset(c('t_dalycare_diagnoses', 'patient', 'RKKP_LYFO' ))

RKKP_LYFO <- readRDS('/ngc/projects2/dalyca_r/andkat_r/RT/Data/RKKP_LYFO_oct25.rds')

LYFO_clean = RKKP_LYFO %>% 
 clean_RKKP_LYFO() 
names(LYFO_clean)

setwd("/ngc/projects2/dalyca_r/sangho_r/data")
cox_df1 <- fread("CD20_final_data_SG/CD20_df.txt")
colSums(is.na(cox_df1))
n_distinct(cox_df1)
names(cox_df1)

#Survival time####
data1 <- cox_df1 %>%
  mutate(date_treatment = if_else(is.na(date_treatment_2nd_line), date_ASCT_2nd_line, date_treatment_2nd_line))
sum(is.na(data1$date_treatment))  

cox_df2 <- data1 %>% 
  mutate(across(c(date_treatment, date_death_fu, date_pato),as.Date),
         start_date = pmax(date_treatment, date_pato, na.rm = TRUE),
         exit_date = date_death_fu,
         surv_time_relapse = as.numeric(difftime(exit_date,start_date,units = "days")))%>% 
  filter(!is.na(surv_time_relapse)& surv_time_relapse>=0)

cox_df2$TTR_group_start <- ifelse(cox_df2$TTR_group_start %in% c("POD12", "POD24"),
                                  "Early",
                                  "Late")
#Follow-up time####
cox_df2 <- cox_df2 %>% 
  mutate(surv_time_years = surv_time_relapse /365.25)

quantile_followup <- quantile(cox_df2$surv_time_years, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

quantile_followup

#Making df ready####
cox_df <- cox_df2 %>% 
  mutate(
    CD20 = factor(CD20, levels = c("Positive", "Negative")),
    Subtype = factor(subtype, levels = c("DLBCL", "FL", "MCL")),
    IPI = factor(IPI, levels = c("Low", "Intermediate", "High")),
    Sex = factor(sex, levels = c("F", "M")),
    Response = factor(response, levels = c("CR", "PR", "SD", "PD")),
    Polypharmacy = factor(ATC_category, levels = c("<5","≥5" )),
    ECOG = factor(PS_diagnosis, levels = c("0", "1", "2", "3", "4")),
    Age = age_at_2L,
    TT2L=factor(cox_df2$TTR_group_start, levels = c("Late", "Early")))
colSums(is.na(cox_df))
#MAIN MODEL Multi Cox IPI####
#ALle subtyper IPI
cox_all_subtype_IPI <- (coxph(Surv(surv_time_relapse, dead)~CD20 + Age + Sex + IPI, data=cox_df)) 
publish (cox_all_subtype_IPI)


#DLBCL, FL, MCL
cox_DLBCL_IPI <- (coxph(Surv(surv_time_relapse, dead)~CD20 + Age + Sex + IPI, data=subset (cox_df, subtype == "DLBCL")))
publish(cox_DLBCL_IPI)

cox_FL_IPI <- (coxph(Surv(surv_time_relapse, dead)~CD20 + Age + Sex + IPI, data=subset (cox_df, subtype == "FL")))
publish(cox_FL_IPI)


cox_MCL_IPI <- (coxph(Surv(surv_time_relapse, dead)~CD20 + Age + Sex + IPI, data=subset (cox_df, subtype == "MCL")))
publish(cox_MCL_IPI)

#MANUALLY REMAKE GGFOREST FUNCTION####
library(broom)
library(grid)

.get_data <- function (fit, data = NULL, complain = TRUE) 
{
  if (is.null(data)) {
    if (complain) 
      warning("The `data` argument is not provided. Data will be extracted from model fit.")
    data <- eval(fit$call$data)
    if (is.null(data)) 
      stop("The `data` argument should be provided either to ggsurvfit or survfit.")
  }
  data
}

ggforestfixed <- function (model, data = NULL, main = "Hazard ratio", cpositions = c(0.02, 
                                                                    0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2, hr_xlim = NULL) 
{
  
  conf.high <- conf.low <- estimate <- NULL
  stopifnot(inherits(model, "coxph"))
  data <- .get_data(model, data = data)
  terms <- attr(model$terms, "dataClasses")[-1]
  coef <- as.data.frame(tidy(model, conf.int = TRUE))
  gmodel <- glance(model)
  allTerms <- lapply(seq_along(terms), function(i) {
    var <- names(terms)[i]
    if (terms[i] %in% c("factor", "character")) {
      adf <- as.data.frame(table(data[, var]))
      cbind(var = var, adf, pos = 1:nrow(adf))
    }
    else if (terms[i] == "numeric") {
      data.frame(var = var, Var1 = "", Freq = nrow(data), 
                 pos = 1)
    }
    else {
      vars = grep(paste0("^", var, "*."), coef$term, value = TRUE)
      data.frame(var = vars, Var1 = "", Freq = nrow(data), 
                 pos = seq_along(vars))
    }
  })
  allTermsDF <- do.call(rbind, allTerms)
  colnames(allTermsDF) <- c("var", "level", "N", "pos")
  inds <- apply(allTermsDF[, 1:2], 1, paste0, collapse = "")
  rownames(coef) <- gsub(coef$term, pattern = "`", replacement = "")
  toShow <- cbind(allTermsDF, coef[inds, ])[, c("var", "level", 
                                                "N", "p.value", "estimate", "conf.low", "conf.high", 
                                                "pos")]
  toShowExp <- toShow[, 5:7]
  toShowExp[is.na(toShowExp)] <- 0
  toShowExp <- format(exp(toShowExp), digits = noDigits)
  toShowExpClean <- data.frame(toShow, pvalue = signif(toShow[, 
                                                              4], noDigits + 1), toShowExp)
  toShowExpClean$stars <- paste0(round(toShowExpClean$p.value, 
                                       noDigits + 1), " ", ifelse(toShowExpClean$p.value < 
                                                                    0.05, "*", ""), ifelse(toShowExpClean$p.value < 0.01, 
                                                                                           "*", ""), ifelse(toShowExpClean$p.value < 0.001, "*", 
                                                                                                            ""))
  toShowExpClean$ci <- paste0("(", toShowExpClean[, "conf.low.1"], 
                              " - ", toShowExpClean[, "conf.high.1"], ")")
  toShowExpClean$estimate.1[is.na(toShowExpClean$estimate)] = refLabel
  toShowExpClean$stars[which(toShowExpClean$p.value < 0.001)] = "<0.001 ***"
  toShowExpClean$stars[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$ci[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$estimate[is.na(toShowExpClean$estimate)] = 0
  toShowExpClean$var = as.character(toShowExpClean$var)
  toShowExpClean$var[duplicated(toShowExpClean$var)] = ""
  toShowExpClean$N <- paste0("(N=", toShowExpClean$N, ")")
  toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1, 
  ]
  
  ### FIX TO BE ABLE TO SET XLIM
  if (is.null(hr_xlim)) {
    rangeb <- range(toShowExpClean$conf.low, toShowExpClean$conf.high, na.rm = T)
  } else {
    rangeb <- log(hr_xlim) # on log scale so that you can just set it and it works (probably)
  }
  breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
  rangeplot <- rangeb
  rangeplot[1] <- rangeplot[1] - diff(rangeb)
  rangeplot[2] <- rangeplot[2] + 0.15 * diff(rangeb)
  width <- diff(rangeplot)
  y_variable <- rangeplot[1] + cpositions[1] * width
  y_nlevel <- rangeplot[1] + cpositions[2] * width
  y_cistring <- rangeplot[1] + cpositions[3] * width
  y_stars <- rangeb[2]
  x_annotate <- seq_len(nrow(toShowExpClean))
  annot_size_mm <- fontsize * as.numeric(convertX(unit(theme_get()$text$size, 
                                                       "pt"), "mm"))
  p <- ggplot(toShowExpClean, aes(seq_along(var), exp(estimate))) + 
    geom_rect(aes(xmin = seq_along(var) - 0.5, xmax = seq_along(var) + 
                    0.5, ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]), 
                  fill = ordered(seq_along(var)%%2 + 1))) + scale_fill_manual(values = c("#FFFFFF33", 
                                                                                         "#00000033"), guide = "none") + geom_point(pch = 15, 
                                                                                                                                    size = 4) + geom_errorbar(aes(ymin = exp(conf.low), 
                                                                                                                                                                  ymax = exp(conf.high)), width = 0.15) + geom_hline(yintercept = 1, 
                                                                                                                                                                                                                     linetype = 3) + coord_flip(ylim = exp(rangeplot)) + 
    ggtitle(main) + scale_y_log10(name = "", labels = sprintf("%g", 
                                                              breaks), expand = c(0.02, 0.02), breaks = breaks) + 
    theme_light() + theme(panel.grid.minor.y = element_blank(), 
                          panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), 
                          legend.position = "none", panel.border = element_blank(), 
                          axis.title.y = element_blank(), axis.text.y = element_blank(), 
                          axis.ticks.y = element_blank(), plot.title = element_text(hjust = 0.5)) + 
    xlab("") + annotate(geom = "text", x = x_annotate, y = exp(y_variable), 
                        label = toShowExpClean$var, fontface = "bold", hjust = 0, 
                        size = annot_size_mm) + annotate(geom = "text", x = x_annotate, 
                                                         y = exp(y_nlevel), hjust = 0, label = toShowExpClean$level, 
                                                         vjust = -0.1, size = annot_size_mm) + annotate(geom = "text", 
                                                                                                        x = x_annotate, y = exp(y_nlevel), label = toShowExpClean$N, 
                                                                                                        fontface = "italic", hjust = 0, vjust = ifelse(toShowExpClean$level == 
                                                                                                                                                         "", 0.5, 1.1), size = annot_size_mm) + annotate(geom = "text", 
                                                                                                                                                                                                         x = x_annotate, y = exp(y_cistring), label = toShowExpClean$estimate.1, 
                                                                                                                                                                                                         size = annot_size_mm, vjust = ifelse(toShowExpClean$estimate.1 == 
                                                                                                                                                                                                                                                "reference", 0.5, -0.1)) + annotate(geom = "text", 
                                                                                                                                                                                                                                                                                    x = x_annotate, y = exp(y_cistring), label = toShowExpClean$ci, 
                                                                                                                                                                                                                                                                                    size = annot_size_mm, vjust = 1.1, fontface = "italic") + 
    annotate(geom = "text", x = x_annotate, y = exp(y_stars), 
             label = toShowExpClean$stars, size = annot_size_mm, 
             hjust = -0.2, fontface = "italic") + annotate(geom = "text", 
                                                           x = 0.5, y = exp(y_variable), label = paste0("# Events: ", 
                                                                                                        gmodel$nevent, "; Global p-value (Log-Rank): ", 
                                                                                                        format.pval(gmodel$p.value.log, eps = ".001"), " \nAIC: ", 
                                                                                                        round(gmodel$AIC, 2), "; Concordance Index: ", round(gmodel$concordance, 
                                                                                                                                                             2)), size = annot_size_mm, hjust = 0, vjust = 1.2, 
                                                           fontface = "italic")
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  ggpubr::as_ggplot(gt)
}

#Main model Forest####

hr_lims <- c(0.5, 7)

A <- ggforestfixed(cox_all_subtype_IPI, data= model.frame(cox_all_subtype_IPI), main = "A- All subtypes", noDigits = 3, hr_xlim = hr_lims)+
  theme(text = element_text(size = 8))

B <- ggforestfixed(cox_DLBCL_IPI, data= model.frame(cox_DLBCL_IPI), main = "B- DLBCL", noDigits = 3, hr_xlim = hr_lims)+
  theme(text = element_text(size = 8))

C<- ggforestfixed(cox_FL_IPI, data= model.frame(cox_FL_IPI),main = "C- FL", noDigits = 3, hr_xlim = hr_lims)+
  theme(text = element_text(size = 8))

D <- ggforestfixed(cox_MCL_IPI, data= model.frame(cox_MCL_IPI), main = "D- MCL", noDigits = 3, hr_xlim = hr_lims)+
  theme(text = element_text(size = 8))
combined_plot <-  (A| B)/(C|D)

ggsave('/ngc/projects2/dalyca_r/sangho_r/figures/forestplot_IPI.png',
       plot = combined_plot,
       dpi=300,
       height = 10,
       width = 12) 

ggsave('/ngc/projects2/dalyca_r/sangho_r/figures/forestplot_IPI.pdf',
       plot = combined_plot,
       dpi=300,
       height = 10,
       width = 12) 


#Cox Sensitivity Models
#Sensitivity 1####
cox_all_subtype_s1 <- (coxph(Surv(surv_time_relapse, dead)~CD20 + Age + Sex + IPI+ Response ,
                          data=cox_df, ties = "efron"))
publish (cox_all_subtype_s1)

#DLBCL, FL, MCL
cox_DLBCL_s1 <- (coxph(Surv(surv_time_relapse, dead)~CD20 + Age + Sex + IPI+
                      Response, data=subset (cox_df, subtype == "DLBCL"), ties = "efron"))
publish(cox_DLBCL_s1)

cox_FL_s1 <- (coxph(Surv(surv_time_relapse, dead)~CD20 + Age + Sex + IPI+
                   Response, data=subset (cox_df, subtype == "FL"), ties = "efron"))
publish(cox_FL_s1)

cox_MCL_s1 <- (coxph(Surv(surv_time_relapse, dead)~CD20 + Age + Sex + IPI+
                    Response, data=subset (cox_df, subtype == "MCL"), ties = "efron"))
publish(cox_MCL_s1)


#Forest S1

AS1 <- ggforestfixed(cox_all_subtype_s1, data= model.frame(cox_all_subtype_s1), main = "A- All subtypes", noDigits = 3, hr_xlim = hr_lims)+
  theme(text = element_text(size = 8))

BS1 <- ggforestfixed(cox_DLBCL_s1, data= model.frame(cox_DLBCL_s1), main = "B- DLBCL" , noDigits = 3, hr_xlim = hr_lims)+
  theme(text = element_text(size = 8))

CS1<- ggforestfixed(cox_FL_s1, data= model.frame(cox_FL_s1),main = "C- FL", noDigits = 3, hr_xlim = hr_lims)+
  theme(text = element_text(size = 8))

DS1<- ggforestfixed(cox_MCL_s1, data= model.frame(cox_MCL_s1), main = "D- MCL" , noDigits = 3, hr_xlim = hr_lims)+
  theme(text = element_text(size = 8))
combined_plot_S1 <- (AS1| BS1)/(CS1|DS1) 

ggsave('/ngc/projects2/dalyca_r/sangho_r/figures/forestplot_S1.png',
       plot = combined_plot_S1,
       dpi=300,
       height = 10,
       width = 12)

ggsave('/ngc/projects2/dalyca_r/sangho_r/figures/forestplot_S1.pdf',
       plot = combined_plot_S1,
       dpi=300,
       height = 10,
       width = 12) 

#Sensitivity 2####
#DLBCL, FL, MCL

cox_all_subtype_s2 <- (coxph(Surv(surv_time_relapse, dead)~CD20 + Age + Sex + IPI+ TT2L ,
                             data=cox_df)) 
publish (cox_all_subtype_s2)


cox_DLBCL_s2 <- (coxph(Surv(surv_time_relapse, dead)~CD20 + Age + Sex + IPI+
                         TT2L, data=subset (cox_df, subtype == "DLBCL")))
publish(cox_DLBCL_s2)

cox_FL_s2 <- (coxph(Surv(surv_time_relapse, dead)~CD20 + Age + Sex + IPI+
                      TT2L, data=subset (cox_df, subtype == "FL")))
publish(cox_FL_s2)

cox_MCL_s2 <- (coxph(Surv(surv_time_relapse, dead)~CD20 + Age + Sex +IPI+ 
                       TT2L, data=subset (cox_df, subtype == "MCL")))
publish(cox_MCL_s2)

#Forest S2

AS2 <- ggforestfixed(cox_all_subtype_s2, data= model.frame(cox_all_subtype_s2), main = "A- All subtypes", noDigits = 3, hr_xlim = hr_lims )+
  theme(text = element_text(size = 8))

BS2 <- ggforestfixed(cox_DLBCL_s2, data= model.frame(cox_DLBCL_s2), main = "B- DLBCL", noDigits = 3, hr_xlim = hr_lims )+
  theme(text = element_text(size = 8))

CS2<- ggforestfixed(cox_FL_s2, data= model.frame(cox_FL_s2),main = "C- FL", noDigits = 3, hr_xlim = hr_lims )+
  theme(text = element_text(size = 8))

DS2<- ggforestfixed(cox_MCL_s2, data= model.frame(cox_MCL_s2), main = "D- MCL", noDigits = 3, hr_xlim = hr_lims )+
  theme(text = element_text(size = 8))
combined_plot_S2 <- (AS2| BS2)/(CS2|DS2) 

ggsave('/ngc/projects2/dalyca_r/sangho_r/figures/forestplot_S2.png',
       plot = combined_plot_S2,
       dpi=300,
       height = 10,
       width = 12
)

ggsave('/ngc/projects2/dalyca_r/sangho_r/figures/forestplot_S2.pdf',
       plot = combined_plot_S2,
       dpi=300,
       height = 10,
       width = 12) 

#KM Main model####
surv_obj <- Surv(time = cox_df$surv_time_relapse, event = cox_df$dead)
cox_df$surv_time_relapse %>% summary
levels(cox_df$CD20)
n_distinct(cox_df)

#KM plot for de forskellige subtyper

#All subtypes
fig2A <- KM_plot(survfit(Surv(surv_time_relapse/365.25, dead) ~ CD20, 
                         data = cox_df),
                 breaks = 2,
                 xlim = c(0,10),
                 title = 'A- All subtypes',
                 labs = c('Positive', 'Negative'),
                 pval = TRUE,
                 xlab= "Time (years from second-line therapy)",
                 ylab= "Overall survival probability")

#For visual abstract
fig2Avisuelabstract <- KM_plot(survfit(Surv(surv_time_relapse/365.25, dead) ~ CD20, 
                                       data = cox_df),
                               breaks = 1,
                               xlim = c(0,5),
                               title = 'CD20',
                               labs = c('Positive', 'Negative'),
                               pval = F,
                               surv.median.line = "hv",
                               xlab= "Time (years from second-line therapy)",
                               ylab= "Overall survival probability")


#different subtypes
fig2B <- KM_plot(survfit(Surv(surv_time_relapse/365.25, dead) ~ CD20, 
                         data = cox_df %>% filter(subtype == "DLBCL")),
                 breaks = 2,
                 xlim = c(0,10),
                 title = 'B- DLBCL',
                 labs = c('Positive', 'Negative'),
                 pval = T,
                 xlab= "Time (years from second-line therapy)",
                 ylab= "Overall survival probability")


fig2C <- KM_plot(survfit(Surv(surv_time_relapse/365.25, dead) ~ CD20, 
                         data = cox_df %>% filter(subtype == "FL")),
                 breaks = 2,
                 xlim = c(0,10),
                 title = 'C- FL',
                 labs = c('Positive', 'Negative'),
                 pval = T,
                 xlab= "Time (years from second-line therapy)",
                 ylab= "Overall survival probability")

fig2D <- KM_plot(survfit(Surv(surv_time_relapse/365.25, dead) ~ CD20, 
                         data = cox_df %>% filter(subtype == "MCL")),
                 breaks = 2,
                 xlim = c(0,10),
                 title = 'D- MCL',
                 labs = c('Positive', 'Negative'),
                 pval = T,
                 xlab= "Time (years from second-line therapy)",
                 ylab= "Overall survival probability")

library("gridExtra")
grid.arrange(fig2A$plot, fig2B$plot, fig2C$plot, fig2D$plot, nrow = 2, ncol=2)

ggsave('/ngc/projects2/dalyca_r/sangho_r/figures/figure_1_os.png',
       arrange_ggsurvplots(list(fig2A, fig2C, fig2B, fig2D), nrow = 2, ncol = 4),
       dpi=300,
       height = 10,
       width = 14
) 

ggsave('/ngc/projects2/dalyca_r/sangho_r/figures/KM_figure_2.pdf',
       arrange_ggsurvplots(list(fig2A, fig2C, fig2B, fig2D), nrow = 2, ncol = 4),
       dpi=300,
       height = 10,
       width = 18
) 

ggsave('/ngc/projects2/dalyca_r/sftp/fromNGC/figure_2A_abstract.pdf',
       arrange_ggsurvplots(list(fig2Avisuelabstract), nrow = 1, ncol=2),
       dpi=300,
       height = 7.5,
       width = 18*0.75) 

fit <- survfit(Surv(surv_time_relapse/365.25, dead) ~ CD20, 
               data = cox_df)
summary(fit, quantile=c(0.25, 0.5, 0.75))$table



#Sanity check
summary(cox_df$surv_time_relapse/365.25)
summary((cox_df$surv_time_relapse/365.25)[cox_df$dead == 0])
summary((cox_df$surv_time_relapse/365.25)[cox_df$dead == 1])

#KM sensitivity Supplemental F1####
cox_df3 <- cox_df %>% 
  mutate(CD20_3cat= factor(CD20_3cat, levels = c("Positive", "Reduced", "Negative")))
#DLBCL:
A <- KM_plot(survfit(Surv(surv_time_relapse/365.25, dead) ~ CD20_3cat, 
                     data = cox_df3),
             breaks = 2,
             xlim = c(0,10),
             title = 'A- All subtypes',
             labs = c('Positive','Reduced','Negative'),
             pval = TRUE,
             xlab= "Time (years from second-line therapy)",
             ylab= "Overall survival probability",
             palette = c("#C41130", "#36B44B","#135284"))



B <- KM_plot(survfit(Surv(surv_time_relapse/365.25, dead) ~ CD20_3cat, 
                     data = cox_df3 %>% filter(subtype == "DLBCL")),
             breaks = 2,
             xlim = c(0,10),
             title = 'B- DLBCL',
             labs = c('Positive','Reduced','Negative'),
             pval = T,
             xlab= "Time (years from second-line therapy)",
             ylab= "Overall survival probability",
             palette = c("#C41130", "#36B44B","#135284"))


C <- KM_plot(survfit(Surv(surv_time_relapse/365.25, dead) ~ CD20_3cat, 
                     data = cox_df3 %>% filter(subtype == "FL")),
             breaks = 2,
             xlim = c(0,10),
             title = 'C- FL',
             labs = c('Positive','Reduced','Negative'),
             pval = T,
             xlab= "Time (years from second-line therapy)",
             ylab= "Overall survival probability",
             palette = c("#C41130", "#36B44B","#135284"))

D <-  KM_plot(survfit(Surv(surv_time_relapse/365.25, dead) ~ CD20_3cat, 
                     data = cox_df3 %>% filter(subtype == "MCL")),
             breaks = 2,
             xlim = c(0,10),
             title = 'D- MCL',
             labs = c('Positive','Reduced','Negative'),
             pval = T,
             xlab= "Time (years from second-line therapy)",
             ylab= "Overall survival probability",
             palette = c("#C41130", "#36B44B","#135284"))


ggsave('/ngc/projects2/dalyca_r/sangho_r/figures/figure_KM_LOW_os.png',
       arrange_ggsurvplots(list(A, C, B, D), nrow = 2, ncol = 4),
       dpi=300,
       height = 10,
       width = 14
) 
ggsave('/ngc/projects2/dalyca_r/sangho_r/figures/figure_KM_LOW_os.pdf',
       arrange_ggsurvplots(list(A, C, B, D), nrow = 2, ncol = 4),
       dpi=300,
       height = 10,
       width = 18
) 
#Pairwise log rank:
survdiff(Surv(surv_time_relapse/365.25, dead)~CD20_3cat, data = cox_df3)

pw_all = pairwise_survdiff(Surv(surv_time_relapse/365.25, dead)~CD20_3cat,
                  data = cox_df3,
                  p.adjust.method= "none") %>%
  tile_pairwise_survdiff(position = "LL", palette=c(1,3,2), labs= F)
#p-værdi for hver subtype
pw_DLBCL<- pairwise_survdiff(Surv(surv_time_relapse/365.25, dead)~CD20_3cat,
                          data = cox_df3 %>%  filter(subtype == "DLBCL"),
                          p.adjust.method= "none")%>%
  tile_pairwise_survdiff(position = "LL", palette=c(1,3,2), labs= F)

pw_FL<- pairwise_survdiff(Surv(surv_time_relapse/365.25, dead)~CD20_3cat,
                  data = cox_df3 %>%  filter(subtype == "FL"),
                  p.adjust.method= "none")%>%
  tile_pairwise_survdiff(position = "LL", palette=c(1,3,2), labs= F)


pw_MCL<- pairwise_survdiff(Surv(surv_time_relapse/365.25, dead)~CD20_3cat,
                          data = cox_df3 %>%  filter(subtype == "MCL"),
                          p.adjust.method= "none")%>%
  tile_pairwise_survdiff(position = "LL", palette=c(1,3,2), labs= F)


ggarrange(pw_all, pw_DLBCL, pw_FL, pw_MCL)

#KM- Sensitivitet without secondary CNS relpas####
n_CNS_relaps <- cox_df3 %>% 
  inner_join(RKKP_LYFO %>% select(patientid, Rec_HavdePatientenCNS),
             by= "patientid") %>% 
  filter(Rec_HavdePatientenCNS == "Y") %>% 
  distinct(patientid, .keep_all = TRUE)
table(n_CNS_relaps$CD20)#Så vi har 121 med CNS relap hvoraf 20 af dem er CD20 negative

cox_df3_uden_CNS_relaps <- cox_df3 %>% 
  anti_join(n_CNS_relaps %>% select(patientid),
            by="patientid")
#All subtypes
Supfig2A <- KM_plot(survfit(Surv(surv_time_relapse/365.25, dead) ~ CD20, 
                         data = cox_df3_uden_CNS_relaps),
                 breaks = 2,
                 xlim = c(0,10),
                 title = 'A- All subtypes',
                 labs = c('Positive', 'Negative'),
                 pval = TRUE,
                 xlab= "Time (years from second-line therapy)",
                 ylab= "Overall survival probability")


#different subtypes
Supfig2B <- KM_plot(survfit(Surv(surv_time_relapse/365.25, dead) ~ CD20, 
                         data = cox_df3_uden_CNS_relaps %>% filter(subtype == "DLBCL")),
                 breaks = 2,
                 xlim = c(0,10),
                 title = 'B- DLBCL',
                 labs = c('Positive', 'Negative'),
                 pval = T,
                 xlab= "Time (years from second-line therapy)",
                 ylab= "Overall survival probability")


Supfig2C <- KM_plot(survfit(Surv(surv_time_relapse/365.25, dead) ~ CD20, 
                         data = cox_df3_uden_CNS_relaps %>% filter(subtype == "FL")),
                 breaks = 2,
                 xlim = c(0,10),
                 title = 'C- FL',
                 labs = c('Positive', 'Negative'),
                 pval = T,
                 xlab= "Time (years from second-line therapy)",
                 ylab= "Overall survival probability")

Supfig2D <- KM_plot(survfit(Surv(surv_time_relapse/365.25, dead) ~ CD20, 
                         data = cox_df3_uden_CNS_relaps %>% filter(subtype == "MCL")),
                 breaks = 2,
                 xlim = c(0,10),
                 title = 'D- MCL',
                 labs = c('Positive', 'Negative'),
                 pval = T,
                 xlab= "Time (years from second-line therapy)",
                 ylab= "Overall survival probability")

library("gridExtra")
grid.arrange(Supfig2A$plot, Supfig2B$plot, Supfig2C$plot, Supfig2D$plot, nrow = 2, ncol=2)

ggsave('/ngc/projects2/dalyca_r/sangho_r/figures/Suppfigure_2_uden_CNSrelaps.png',
       arrange_ggsurvplots(list(Supfig2A, Supfig2C, Supfig2B, Supfig2D), nrow = 2, ncol = 4),
       dpi=300,
       height = 10,
       width = 14
) 
ggsave('/ngc/projects2/dalyca_r/sangho_r/figures/Suppfigure_2_uden_CNSrelaps.pdf',
       arrange_ggsurvplots(list(Supfig2A, Supfig2C, Supfig2B, Supfig2D), nrow = 2, ncol = 4),
       dpi=300,
       height = 10,
       width = 18
) 

