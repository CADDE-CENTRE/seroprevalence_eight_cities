library('nimble')
library("rstudioapi")
library(tidyverse)
library(lubridate)
library(rstudioapi)
library(reshape2)
library(Hmisc)
setwd(dirname(getActiveDocumentContext()$path))

source("serorev_functions.R")

BSC_prev_val <- function(center, ppos, suffix = "", Nsamples_val = 1000, thr = 0.49, week = FALSE) 
#This function runs the seroreversion correction model
{
  print(paste0(center, "_", suffix))
  if(center == "HEMOAM")
    city_code <- 1302603
  else if(center == "FPS"){
    city_code <- 3550308
  }
  else if(center == "HEMOCE") {
    city_code <- 2304400
  }
  else if(center == "HEMOBA") {
    city_code <- 2927408
  }
  else if(center == "HEMOMINAS") {
    city_code <- 3106200
  }
  else if(center == "HEMOPE") {
    city_code <- 2611606
  }
  else if(center == "HEMEPAR") {
    city_code <- 4106902
  }
  else if(center == "HEMORIO") {
    city_code <- 3304557
  }
  else{
    stop("Invalid center. Choose FPS or HEMOAM")
  }
  all_centres <- readRDS("data/Bloodbank.rds")
  
  pop <- readRDS("data/population2020_eightcities.rds") %>% filter(code == city_code) %>% select(-code) %>% arrange(age_sex) 
  
  all_centres <- all_centres %>% filter(!is.na(idade), idade != "", !is.na(sex), sex != "", idade >= 16, idade <= 70)
  all_centres$age_group <- cut(all_centres$idade, c(0, 25, 35, 45, 55, 100), 
                               labels = c("15-24", "25-34", "35-44", "45-54", "55-69"))
  all_centres$age_sex <- paste0(all_centres$sex, "_", all_centres$age_group)
  
  all_centres <- all_centres %>% filter(blood_center == center, month <= 15)
  if(week) {
    dates <- seq(as.Date("2020-01-01"), as.Date("2021-04-08"), by = "week")
    all_centres$month <- cut(all_centres$donation_date, dates, label = 1:(length(dates)-1)) %>% as.numeric() #We keep the name month, 
                                                                                                             #but it now represents the number of weeks 
                                                                                                             #after January 1st 2020
    
  }
  
  prev <- all_centres %>% group_by(month, age_sex) %>% summarise(ntests = n(), npos = sum(result >= thr))
  prev <- expand.grid(month = 10:max(prev$month), age_sex = unique(prev$age_sex)) %>% left_join(prev, by = c("month", "age_sex"))
  prev$ntests[is.na(prev$ntests)] <- 0
  prev$npos[is.na(prev$npos)] <- 0
  prev <- prev %>% arrange(month, age_sex)
  
  res <- prevalence.blood_donors(prev, ppos, pop$population, cases_max = 2, nburnin = 50000, niter = 100000)
  months <- prev %>% arrange(month) %>% .$month %>% unique()
  
  colnames(res$prev_est) <- months
  colnames(res$cases_est) <- months
  
  prev_est <- res$prev_est %>% data.frame() %>% melt()
  
  write.csv(summarise.samples(res$prev_est) %>% select(-age_sex), paste0("outputs/prev_est_", suffix, ".csv"))
  write.csv(summarise.samples(res$cases_est) %>% select(-age_sex), paste0("outputs/cases_est_", suffix, ".csv"))
  write.csv(summarise.samples(res$drho), paste0("outputs/drho_est_", suffix, ".csv"))
  write.csv(summarise.samples(apply(res$drho, MARGIN = c(2,3), FUN = cumsum)), paste0("outputs/rho_est_", suffix, ".csv"))
  
  saveRDS(res$samples, paste0("outputs/samples_", suffix, ".rds"))
  
  # --- Validation --- #
  
  prev_validation <- simulate.epidemic(prev, res$drho[,,sample(1:dim(res$drho)[3], Nsamples_val, replace = FALSE)], ppos, 1)
  prev_validation <- prev_validation %>% mutate(prev_mes = npos/ntests)
  prev <- prev %>% mutate(prev_mes = npos/ntests) %>% group_by(month, age_sex) %>% 
    mutate(cil = binconf(x = npos, n = ntests)[2], ciu = binconf(x = npos, n = ntests)[3]) %>% ungroup()
  
  dfplot <- rbind(prev %>% select(month, age_sex, prev_mes, cil, ciu) %>% mutate(type = "Measured"), 
                  prev_validation %>% group_by(month, age_sex) %>% summarise(
                    cil = quantile(prev_mes, 0.025, na.rm = TRUE), 
                    ciu = quantile(prev_mes, 0.975, na.rm = TRUE),
                    prev_mes = mean(prev_mes, na.rm = TRUE)) %>% 
                    mutate(type = "Reconstructed"))
  
  dfplot$month = as.factor(dfplot$month)
  
  write.csv(dfplot, file = paste0("outputs/validation_", suffix, ".csv"))
  
  # --- Validation: All age-sex groups together --- #
  
  prev_all <- prev %>% group_by(month) %>% summarise(npos = sum(npos), ntests = sum(ntests)) %>% 
    group_by(month) %>% mutate(prev_mes = npos/ntests, cil = binconf(x = npos, n = ntests)[2], ciu = binconf(x = npos, n = ntests)[3]) %>%
    ungroup()
  
  prev_validation_all <- prev_validation %>% group_by(month, id_sample) %>% summarise(prev_mes = sum(npos)/sum(ntests))
  
  dfplot_all <- rbind(prev_all %>% select(month, prev_mes, cil, ciu) %>% mutate(type = "Measured"), 
                      prev_validation_all %>% group_by(month) %>% summarise(
                        cil = quantile(prev_mes, 0.025, na.rm = TRUE), 
                        ciu = quantile(prev_mes, 0.975, na.rm = TRUE),
                        prev_mes = mean(prev_mes, na.rm = TRUE)) %>% 
                        mutate(type = "Reconstructed"))
  
  write.csv(dfplot_all, file = paste0("outputs/validation_all_", suffix, ".csv"))
  
}

# --- Run the model --- #


ppos_df <- read.csv("data/ppos_week.csv")
Nsamples_val <- 1000 #Number of samples used in the validation analysis

for(center in c("HEMOAM", "FPS", "HEMOCE", "HEMEPAR", "HEMORIO", "HEMOBA", "HEMOPE", "HEMOMINAS")) {
  print(center)
  ppos <- ppos_df %>% filter(type == "Repeat Donors") %>% arrange(date) %>% .$ppos #p+[n] estimated with repeat donors (approach used in the paper)
  BSC_prev_val(center, ppos, paste0(center, "_repeat_week_agesex"), Nsamples_val, 0.49, TRUE)
  
  ppos <- ppos_df %>% filter(type == "Repeat Donors (deaths)") %>% arrange(date) %>% .$ppos #p+[n] estimated with repeat donors with incidence estimated using deaths 
  BSC_prev_val(center, ppos, paste0(center, "_deaths_week_agesex"), Nsamples_val, 0.49, TRUE)
  
  ppos <- ppos_df %>% filter(type == "Plasma Donors") %>% arrange(date) %>% .$ppos          #p+[n] estimated with plasma donors
  BSC_prev_val(center, ppos, paste0(center, "_plasma_week_agesex"), Nsamples_val, 0.49, TRUE)
}