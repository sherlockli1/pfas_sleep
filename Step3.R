#Protein mediation
library(mediation)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(dplyr)
library(gtsummary)
library(lme4)
library(geepack)
library(qgcomp)
library(stats)
library(randomForest)
library(pROC)
library(caret)
library(grpreg)
library(ranger)
library(mice)
library(ggpubr)

#Load data and gene list
set.seed(100)
#Mediation analysis
mediator_list<-colnames(baseline_protein_pfda)
results_list <- list()
extract_mediation_summary <- function (x) { 
  
  clp <- 100 * x$conf.level
  isLinear.y <- ((class(x$model.y)[1] %in% c("lm", "rq")) || 
                   (inherits(x$model.y, "glm") && x$model.y$family$family == 
                      "gaussian" && x$model.y$family$link == "identity") || 
                   (inherits(x$model.y, "survreg") && x$model.y$dist == 
                      "gaussian"))
  
  printone <- !x$INT && isLinear.y
  
  if (printone) {
    
    smat <- c(x$d1, x$d1.ci, x$d1.p)
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
    
    rownames(smat) <- c("ACME", "ADE", "Total Effect", "Prop. Mediated")
    
  } else {
    smat <- c(x$d0, x$d0.ci, x$d0.p)
    smat <- rbind(smat, c(x$d1, x$d1.ci, x$d1.p))
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$z1, x$z1.ci, x$z1.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
    smat <- rbind(smat, c(x$n1, x$n1.ci, x$n1.p))
    smat <- rbind(smat, c(x$d.avg, x$d.avg.ci, x$d.avg.p))
    smat <- rbind(smat, c(x$z.avg, x$z.avg.ci, x$z.avg.p))
    smat <- rbind(smat, c(x$n.avg, x$n.avg.ci, x$n.avg.p))
    
    rownames(smat) <- c("ACME (control)", "ACME (treated)", 
                        "ADE (control)", "ADE (treated)", "Total Effect", 
                        "Prop. Mediated (control)", "Prop. Mediated (treated)", 
                        "ACME (average)", "ADE (average)", "Prop. Mediated (average)")
    
  }
  
  colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep = ""), 
                      paste(clp, "% CI Upper", sep = ""), "p-value")
  smat
  
}

for (mediator_var in mediator_list) {
  # Create the formula for the mediation models
  m_formula <- as.formula(paste(mediator_var, "~ tertile_ma_pfda + ma_race_eth + ma_sex + ma_age + ma_hei + ma_parental_edu_2cat + ma_cig_life + ma_alc + ma_self_pa"))
  y_formula <- as.formula(paste("ma_sleep_wkday_durations ~", mediator_var, "+ tertile_ma_pfda + ma_race_eth + ma_sex + ma_age + ma_hei + ma_parental_edu_2cat + ma_cig_life + ma_alc + ma_self_pa"))
  
  # Fit the mediation models
  m_model <- lm(m_formula, data = chs_new_ma_baseline_tertiles_pfda)
  y_model <- lm(y_formula, data = chs_new_ma_baseline_tertiles_pfda)
  
  # Perform mediation analysis
  contcont.cont <- mediate(m_model, y_model,sims = 50,
                           treat = "tertile_ma_pfda", mediator = mediator_var,
                           data = chs_new_ma_baseline_tertiles_pfda, control.value = 0, treat.value = 1)
  
  # Store the results in the list
  results_list[[mediator_var]] <- extract_mediation_summary(summary(contcont.cont))

}

combined_models <- do.call(rbind, results_list)
mediator_variable_column <- rep(mediator_list, each = nrow(results_list[[1]]))

combined_models_mediation_1 <- as.data.frame(cbind(Mediator = mediator_variable_column, combined_models))

#Repeat for other PFAS-sleep associations