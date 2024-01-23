library(tidyverse)
library(reshape2)
library(ggplot2)
library(dplyr)
library(geepack)
library(qgcomp)
library(stats)
library(ranger)
library(mice)
library(ggpubr)


#Run models
exp_baseline_list_new<-c(
  "tertile_ma_pfda",                     
  "tertile_ma_pfhps",                     
  "tertile_ma_pfhxs",                      
  "tertile_ma_pfna",                      
  "tertile_ma_pfoa",                       
  "tertile_ma_pfos",                      
  "tertile_ma_pfpes"
)

formula_base <- "ma_sleep_wkday_durations ~ ma_race_eth + ma_sex + ma_age +ma_hei+ma_parental_edu_2cat+ma_cig_life+ma_alc+ma_self_pa"
# Run the logistic regression models using lapply
models <- lapply(exp_baseline_list_new, function(x_var) {
  formula <- paste(formula_base, x_var, sep = " + ")
  model <- summary(glm(as.formula(formula), data = chs_new_ma_baseline_tertiles, family = gaussian))$coefficients[, 1:4]
})

combined_models <- do.call(rbind, models)
exposure_variable_column <- rep(exp_baseline_list_new, each = nrow(models[[1]]))
combined_models_baseline_duration <- as.data.frame(cbind(Exposure_Variable = exposure_variable_column, combined_models))

rows_to_keep <- rep(FALSE, nrow(combined_models_baseline_duration))
rows_to_keep[seq(11, nrow(combined_models_baseline_duration), by = 11)] <- TRUE
combined_models_baseline_duration <- combined_models_baseline_duration[rows_to_keep, ]


#PFAS mixture: qgcomp
mixture_sleep_duration <- matrix(nrow = 1, ncol = 5)

exposures <- c("ma_pfda","ma_pfhps","ma_pfhxs","ma_pfna","ma_pfoa","ma_pfos","ma_pfpes")
outcome <- "ma_sleep_wkday_durations"
covariates <- c("ma_race_eth","ma_sex","ma_age","ma_hei","ma_parental_edu_2cat","ma_cig_life","ma_alc","ma_self_pa")

formula <- as.formula(paste(outcome, "~", paste(c(exposures, covariates), collapse = " + ")))

# Run the qgcomp model with bootstrapping
mixture_baseline <- qgcomp.glm.boot(formula,
                                    family = gaussian(),
                                    expnms = exposures,
                                    data = chs_new_clean,
                                    q = 3,
                                    B = 100)  
mixture_sleep_duration[1,]<-t(summary(mixture_baseline)$coefficients[2,])



#Create plots
p1<-ggplot(combined_models_baseline_linear_duration, aes(x = Exposure_Variable, y = RR, color = level)) +
  geom_point(position = position_dodge(width = 0.5),size=2.5) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, size=0.8, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
  ylab("Beta (95% CI)") +
  xlab("Baseline") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),text= element_text(size=14))+
  ylim(-2, 2)+
  scale_color_manual(values = c("#D55E00","#0072b2"))
p1




#Repeat for follow-up visit
exp_baseline_list_new<-c("tertile_mc_pfda",
                         "tertile_mc_pfhps",     
                        "tertile_mc_pfhxs",     
                        "tertile_mc_pfna",      
                        "tertile_mc_pfoa",     
                        "tertile_mc_pfos",     
                        "tertile_mc_pfpes")
#Run models

formula_base <- "mc_sleep_disturbance_tscore ~ ma_race_eth + ma_sex + mc_age +mc_hei+ma_parental_edu_2cat+mc_cig_life+mc_alc+mc_met_cat"

models <- lapply(exp_baseline_list_new, function(x_var) {
  formula <- paste(formula_base, x_var, sep = " + ")
  model <- summary(glm(as.formula(formula), data = chs_new_mc_follow_tertiles, family = gaussian))$coefficients[, 1:4]
})

combined_models <- do.call(rbind, models)
exposure_variable_column <- rep(exp_baseline_list_new, each = nrow(models[[1]]))
combined_models_followup_disturb <- as.data.frame(cbind(Exposure_Variable = exposure_variable_column, combined_models))

rows_to_keep <- rep(FALSE, nrow(combined_models_followup_disturb))
rows_to_keep[seq(12, nrow(combined_models_followup_disturb), by = 12)] <- TRUE
combined_models_followup_disturb <- combined_models_followup_disturb[rows_to_keep, ]

formula_base <- "mc_sleep_impair_tscore ~ ma_race_eth + ma_sex + mc_age +mc_hei+ma_parental_edu_2cat+mc_cig_life+mc_alc+mc_met_cat"

models <- lapply(exp_baseline_list_new, function(x_var) {
  formula <- paste(formula_base, x_var, sep = " + ")
  model <- summary(glm(as.formula(formula), data = chs_new_mc_follow_tertiles, family = gaussian))$coefficients[, 1:4]
})

combined_models <- do.call(rbind, models)
exposure_variable_column <- rep(exp_baseline_list_new, each = nrow(models[[1]]))
combined_models_followup_impair <- as.data.frame(cbind(Exposure_Variable = exposure_variable_column, combined_models))

rows_to_keep <- rep(FALSE, nrow(combined_models_followup_impair))
rows_to_keep[seq(12, nrow(combined_models_followup_impair), by = 12)] <- TRUE
combined_models_followup_impair <- combined_models_followup_impair[rows_to_keep, ]

formula_base <- "mc_nightly_sleep ~ ma_race_eth + ma_sex + mc_age +mc_hei+ma_parental_edu_2cat+mc_cig_life+mc_alc+mc_met_cat"

models <- lapply(exp_baseline_list_new, function(x_var) {
  formula <- paste(formula_base, x_var, sep = " + ")
  model <- summary(glm(as.formula(formula), data = chs_new_mc_follow_tertiles, family = gaussian))$coefficients[, 1:4]
})

combined_models <- do.call(rbind, models)
exposure_variable_column <- rep(exp_baseline_list_new, each = nrow(models[[1]]))
combined_models_followup_night <- as.data.frame(cbind(Exposure_Variable = exposure_variable_column, combined_models))

rows_to_keep <- rep(FALSE, nrow(combined_models_followup_night))
rows_to_keep[seq(12, nrow(combined_models_followup_night), by = 12)] <- TRUE
combined_models_followup_night <- combined_models_followup_night[rows_to_keep, ]





#Mixture
mixture_followup_sleep_disturb <- matrix(nrow = 1, ncol = 5)

exposures <- c("mc_pfda","mc_pfhps","mc_pfhxs","mc_pfna","mc_pfoa","mc_pfos","mc_pfpes")
outcome <- "mc_sleep_disturbance_tscore"
covariates <- c("ma_race_eth","ma_sex","mc_age","mc_hei","ma_parental_edu_2cat","mc_cig_life","mc_alc","mc_met_cat")

# Construct the formula for the outcome and covariates
formula <- as.formula(paste(outcome, "~", paste(c(exposures, covariates), collapse = " + ")))

# Run the qgcomp model with bootstrapping
mixture_baseline <- qgcomp.glm.boot(formula,
                                    family = gaussian(),
                                    expnms = exposures,
                                    data = chs_new_clean,
                                    q = 3,
                                    B = 100)  
mixture_followup_sleep_disturb[1,]<-t(summary(mixture_baseline)$coefficients[2,])


mixture_followup_sleep_impair <- matrix(nrow = 1, ncol = 5)

exposures <- c("mc_pfda","mc_pfhps","mc_pfhxs","mc_pfna","mc_pfoa","mc_pfos","mc_pfpes")
outcome <- "mc_sleep_impair_tscore"
covariates <- c("ma_race_eth","ma_sex","mc_age","mc_hei","ma_parental_edu_2cat","mc_cig_life","mc_alc","mc_met_cat")

formula <- as.formula(paste(outcome, "~", paste(c(exposures, covariates), collapse = " + ")))

# Run the qgcomp model with bootstrapping
mixture_baseline <- qgcomp.glm.boot(formula,
                                    family = gaussian(),
                                    expnms = exposures,
                                    data = chs_new_clean,
                                    q = 3,
                                    B = 100)  
mixture_followup_sleep_impair[1,]<-t(summary(mixture_baseline)$coefficients[2,])

mixture_followup_sleep_night <- matrix(nrow = 1, ncol = 5)

exposures <- c("mc_pfda","mc_pfhps","mc_pfhxs","mc_pfna","mc_pfoa","mc_pfos","mc_pfpes")
outcome <- "mc_nightly_sleep"
covariates <- c("ma_race_eth","ma_sex","mc_age","mc_hei","ma_parental_edu_2cat","mc_cig_life","mc_alc","mc_met_cat")

formula <- as.formula(paste(outcome, "~", paste(c(exposures, covariates), collapse = " + ")))
# Run the qgcomp model with bootstrapping
mixture_baseline <- qgcomp.glm.boot(formula,
                                    family = gaussian(),
                                    expnms = exposures,
                                    data = chs_new_clean,
                                    q = 3,
                                    B = 100)  
mixture_followup_sleep_night[1,]<-t(summary(mixture_baseline)$coefficients[2,])



#Plot
p2<-ggplot(combined_models_followup_linear_duration, aes(x = Exposure_Variable, y = RR, color = level)) +
  geom_point(position = position_dodge(width = 0.5),size=2.5) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, size=0.8, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.5) +  
  ylab("Beta (95% CI)") +
  xlab("Follow-up") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),text= element_text(size=14),legend.position = "top")+
  ylim(-2, 2)+
  scale_color_manual(values = c("#D55E00","#0072b2"))
p2


p3<-ggplot(combined_models_followup_linear_other, aes(x = Exposure_Variable, y = RR, color = level)) +
  geom_point(position = position_dodge(width = 0.5),size=2.5) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, size=0.8, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
  facet_grid(outcome ~ .) +
  ylab("Beta (95% CI)") +
  xlab("Follow-up") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),text= element_text(size=14),legend.position = "top")+
  scale_color_manual(values = c("#D55E00","#0072b2"))
p3



#binary outcome
exp_baseline_list_new<-c("tertile_ma_pfda",                       
                         "tertile_ma_pfhps",                     
                         "tertile_ma_pfhxs",                      
                         "tertile_ma_pfna",                      
                         "tertile_ma_pfoa",                       
                         "tertile_ma_pfos",                      
                         "tertile_ma_pfpes"  )

formula_base <- "ma_sleep_7 ~ ma_race_eth + ma_sex + ma_age +ma_hei+ma_parental_edu_2cat+ma_cig_life+ma_alc+ma_self_pa"

models <- lapply(exp_baseline_list_new, function(x_var) {
  formula <- paste(formula_base, x_var, sep = " + ")
  model <- summary(glm(as.formula(formula), data = chs_new_ma_baseline_tertiles, family = binomial))$coefficients[, 1:4]
})

combined_models <- do.call(rbind, models)
exposure_variable_column <- rep(exp_baseline_list_new, each = nrow(models[[1]]))
combined_models_baseline_7 <- as.data.frame(cbind(Exposure_Variable = exposure_variable_column, combined_models))

rows_to_keep <- rep(FALSE, nrow(combined_models_baseline_7))
rows_to_keep[seq(11, nrow(combined_models_baseline_7), by = 11)] <- TRUE
combined_models_baseline_7 <- combined_models_baseline_7[rows_to_keep, ]

formula_base <- "mc_sleep_7 ~ ma_race_eth + ma_sex + ma_age +mc_hei+ma_parental_edu_2cat+mc_cig_life+mc_alc+mc_met_cat"
# Run the logistic regression models using lapply
models <- lapply(exp_baseline_list_new, function(x_var) {
  formula <- paste(formula_base, x_var, sep = " + ")
  model <- summary(glm(as.formula(formula), data = chs_new_ma_follow_tertiles, family = binomial))$coefficients[, 1:4]
})
combined_models <- do.call(rbind, models)
exposure_variable_column <- rep(exp_baseline_list_new, each = nrow(models[[1]]))
combined_models_follow_7 <- as.data.frame(cbind(Exposure_Variable = exposure_variable_column, combined_models))

rows_to_keep <- rep(FALSE, nrow(combined_models_follow_7))
rows_to_keep[seq(12, nrow(combined_models_follow_7), by = 12)] <- TRUE
combined_models_follow_7 <- combined_models_follow_7[rows_to_keep, ]


# Mixture
mixture_sleep_7 <- matrix(nrow = 4, ncol = 6)

exposures <- c("ma_pfda","ma_pfhps","ma_pfhxs","ma_pfna","ma_pfoa","ma_pfos","ma_pfpes")
outcome <- "ma_sleep_7"
covariates <- c("ma_race_eth","ma_sex","ma_age","ma_hei","ma_parental_edu_2cat","ma_cig_life","ma_alc","ma_self_pa")

formula <- as.formula(paste(outcome, "~", paste(c(exposures, covariates), collapse = " + ")))

# Run the qgcomp model with bootstrapping
mixture_baseline <- qgcomp.glm.boot(formula,
                                    family = "binomial",
                                    expnms = exposures,
                                    data = chs_new_clean,
                                    q = 3,  
                                    B = 100,  
                                    rr = FALSE)  
mixture_sleep_7[1,]<-t(summary(mixture_baseline)$coefficients[2,])


mixture_sleep_follow_7 <- matrix(nrow = 4, ncol = 6)

exposures <- c("ma_pfda","ma_pfhps","ma_pfhxs","ma_pfna","ma_pfoa","ma_pfos","ma_pfpes")
outcome <- "mc_sleep_7"
covariates <- c("ma_race_eth","ma_sex","ma_age","mc_hei","ma_parental_edu_2cat","mc_cig_life","mc_alc","mc_met_cat")

formula <- as.formula(paste(outcome, "~", paste(c(exposures, covariates), collapse = " + ")))


# Run the qgcomp model with bootstrapping
mixture_baseline <- qgcomp.glm.boot(formula,
                                    family = binomial(),
                                    expnms = exposures,
                                    data = chs_new_clean,
                                    q = 3,
                                    B = 100)  
mixture_sleep_follow_7[1,]<-t(summary(mixture_baseline)$coefficients[2,])


#Plot
p4<-ggplot(combined_models_baseline_logistic, aes(x = Exposure_Variable, y = OR, color = level)) +
  geom_point(position = position_dodge(width = 0.5),size=2.5) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, size=0.8, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 0.5) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(10^-1.2,10^1.5)) +
  ylab("Odds Ratio (95% CI)") +
  xlab("Baseline") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),text= element_text(size=14))+
  scale_color_manual(values = c("#D55E00","#0072b2"))
p4


exp_baseline_list_new<-c(
  "tertile_mc_pfda"  ,
  "tertile_mc_pfhps" ,
  "tertile_mc_pfhxs" ,
  "tertile_mc_pfna"  ,
  "tertile_mc_pfoa"  ,
  "tertile_mc_pfos" ,
  "tertile_mc_pfpes"
)

#sleep 7
formula_base <- "mc_sleep_7 ~ ma_race_eth + ma_sex + mc_age +mc_hei+ma_parental_edu_2cat+mc_cig_life+mc_alc+mc_met_cat"

# Run the logistic regression models using lapply
models <- lapply(exp_baseline_list_new, function(x_var) {
  formula <- paste(formula_base, x_var, sep = " + ")
  model <- summary(glm(as.formula(formula), data = chs_new_mc_follow_tertiles, family = binomial))$coefficients[, 1:4]
})

combined_models <- do.call(rbind, models)
exposure_variable_column <- rep(exp_baseline_list_new, each = nrow(models[[1]]))
combined_models_followup_7 <- as.data.frame(cbind(Exposure_Variable = exposure_variable_column, combined_models))

rows_to_keep <- rep(FALSE, nrow(combined_models_followup_7))
rows_to_keep[seq(12, nrow(combined_models_followup_7), by = 12)] <- TRUE
combined_models_followup_7 <- combined_models_followup_7[rows_to_keep, ]


#Mixture
# chs_new_mc_sleep_disturb
mixture_followup_sleep_7 <- matrix(nrow = 1, ncol = 6)

exposures <- c("mc_pfda","mc_pfhps","mc_pfhxs","mc_pfna","mc_pfoa","mc_pfos","mc_pfpes")
outcome <- "mc_sleep_7"
covariates <- c("ma_race_eth","ma_sex","mc_age","mc_hei","ma_parental_edu_2cat","mc_cig_life","mc_alc","mc_met_cat")

# Construct the formula for the outcome and covariates
formula <- as.formula(paste(outcome, "~", paste(c(exposures, covariates), collapse = " + ")))

# Run the qgcomp model with bootstrapping
mixture_baseline <- qgcomp.glm.boot(formula,
                                    family = binomial,
                                    expnms = exposures,
                                    data = chs_new_clean,
                                    q = 3,
                                    B = 100)  
mixture_followup_sleep_7[1,]<-t(summary(mixture_baseline)$coefficients[2,])

#Plot
p5<-ggplot(combined_models_followup_linear, aes(x = Exposure_Variable, y = OR, color = level)) +
  geom_point(position = position_dodge(width = 0.5),size=2.5) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, size=0.8, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 0.5) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(10^-1.2,10^1.5)) +
  ylab("Odds Ratio (95% CI)") +
  labs(x="Follow-up") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),text= element_text(size=14))+
  scale_color_manual(values = c("#D55E00","#0072b2"))
p5


#Longitudinal
formula_base<-"sleep_duration_new ~ ma_race_eth+ma_sex+ma_parental_edu_2cat+follow_up+hei_new+cig_new+alc_new+mc_met_cat"

models <- lapply(exp_list, function(x_var) {
  formula <- paste(formula_base, x_var, sep = " + ")
  formula <- paste(formula, "follow_up", sep = " * ")
  model <- summary(geeglm(as.formula(formula),id=id,data=chs_long, corstr = "exchangeable"))$coefficients[, 1:4]
})

combined_models <- do.call(rbind, models)


