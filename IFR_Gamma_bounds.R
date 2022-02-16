library(tidyverse)
library(Hmisc)
library(rstudioapi)
library(nimble)
setwd(dirname(getActiveDocumentContext()$path))

source("serorev_functions.R")
source("IFR_functions.R")
source("prev_functions.R")
source("Mamoeiro.R")

#This script produces estimates of the lower and upper bounds for the attack rate 
#and IFR of the Gamma wave in Manaus

cities <- c("HEMOAM", "FPS", "HEMOCE", "HEMOBA", "HEMOMINAS", "HEMOPE", "HEMEPAR", "HEMORIO")
codes <- c("HEMOAM" = 130260, "FPS" = 355030, "HEMOCE" = 230440, "HEMOBA" = 292740, 
           "HEMOMINAS" = 310620, "HEMOPE" = 261160, "HEMEPAR" = 410690, "HEMORIO" = 330455)

all_centres <- readRDS("data/Bloodbank.rds")
all_centres <- all_centres %>% filter(!is.na(idade), idade != "", !is.na(sex), sex != "", idade >= 16, idade <= 70)
all_centres$age_group <- cut(all_centres$idade, c(0, 25, 35, 45, 55, 100), 
                             labels = c("15-24", "25-34", "35-44", "45-54", "55-64"))
all_centres$age_sex <- paste0(all_centres$sex, "_", all_centres$age_group)

L <- list()
for(i in 1:8) {
  center <- names(codes)[i]
  city_code <- codes[i]
  
  thr <- 1.4
  df <- all_centres %>% filter(blood_center == center, month <= 15)
  prev <- df %>% group_by(month, age_sex) %>% summarise(ntests = n(), npos = sum(result >= thr))
  prev <- expand.grid(month = min(prev$month):max(prev$month), age_sex = unique(prev$age_sex)) %>% left_join(prev, by = c("month", "age_sex"))
  prev$ntests[is.na(prev$ntests)] <- 0
  prev$npos[is.na(prev$npos)] <- 0
  prev <- prev %>% arrange(month, age_sex)
  L[[i]] <- prev %>% mutate(center = center)
}
prev <- do.call(rbind, L)
prev <- prev %>% filter(month %in% c(12, 13, 14, 15), center == "HEMOAM")


pop <- readRDS("data/population2020.rds")
pop <- pop %>% mutate(code = floor(code/10)) %>% 
  filter(code %in% codes, age_bin %in% c("15-20", "20-25", "25-30", "30-35", "35-40", "40-45", "45-50", "50-55", "55-60", "60-65"))
agemap <- c("15-20" = "15-24", "20-25" = "15-24", "25-30" = "25-34", "30-35" = "25-34", 
            "35-40" = "35-44", "40-45" = "35-44", "45-50" = "45-54", "50-55" = "45-54", "55-60" = "55-64", "60-65" = "55-64")
pop <- pop %>% mutate(new_bin = agemap[age_bin])
pop <- pop %>% group_by(code, new_bin, CS_SEXO) %>% summarise(population = sum(population)) %>% 
  ungroup() %>% rename(age_bin = new_bin) %>% mutate(age_sex = paste0(CS_SEXO, "_", age_bin)) %>% arrange(age_sex)


srag <- read.csv("D:/Downloads/INFLUD20-14-02-2022.csv", stringsAsFactors = F, sep = ";")
srag21 <- read.csv("D:/Downloads/INFLUD21-14-02-2022.csv", stringsAsFactors = F, sep = ";")
srag21 <- srag21[, colnames(srag)]
srag <- rbind(srag, srag21)

srag <- srag %>% filter(CO_MUN_RES %in% codes, EVOLUCAO == 2) %>%
  mutate(DT_EVOLUCA = as.Date(DT_EVOLUCA, format = "%d/%m/%Y"),
         DT_SIN_PRI = as.Date(DT_SIN_PRI, format = "%d/%m/%Y"),
         age = as.numeric(ifelse(TP_IDADE == 3, NU_IDADE_N, 0)))
srag$CLASSI_FIN[is.na(srag$CLASSI_FIN)] <- 9
srag <- srag %>% filter(CLASSI_FIN %in% c(4,5,9), DT_SIN_PRI >= as.Date('2020-03-01'))

agebins <- c("15-24", "25-34", "35-44", "45-54", "55-64")
srag$age_group <- cut(srag$age, c(15, 25, 35, 45, 55, 65), labels = agebins)
srag <- srag %>% filter(!is.na(age_group), CS_SEXO %in% c("M","F")) %>% mutate(age_sex = paste0(CS_SEXO, "_", age_group))
deaths.mar <- srag %>% filter(DT_SIN_PRI <= as.Date("2021-03-15")) %>% group_by(age_sex, CO_MUN_RES) %>% tally()
deaths.dez <- srag %>% filter(DT_SIN_PRI < as.Date("2020-12-16")) %>% group_by(age_sex, CO_MUN_RES) %>% tally()
prev <- prev  %>% mutate(CO_MUN_RES = codes[center]) %>% 
  left_join(pop %>% rename(CO_MUN_RES = code) %>% select(-c(age_bin, CS_SEXO)), by = c("CO_MUN_RES", "age_sex"))

prev <- prev %>% arrange(month, age_sex)
prev.dez <- prev %>% filter(month == 12) %>% rename(npos.dez = npos, ntests.dez = ntests)
prev.jan <- prev %>% filter(month == 13) %>% select(npos, ntests) %>% rename(npos.jan = npos, ntests.jan = ntests) 
prev.fev <- prev %>% filter(month == 14) %>% select(npos, ntests) %>% rename(npos.fev = npos, ntests.fev = ntests)
prev.mar <- prev %>% filter(month == 15) %>% select(npos, ntests) %>% rename(npos.mar = npos, ntests.mar = ntests)
prev <- cbind(prev.dez, prev.jan, prev.fev, prev.mar)

prev <- prev %>% left_join(deaths.mar %>% rename(deaths.mar = n), by = c("CO_MUN_RES", "age_sex"))
prev <- prev %>% left_join(deaths.dez %>% rename(deaths.dez = n), by = c("CO_MUN_RES", "age_sex"))

prev$deaths <- prev$deaths.mar - prev$deaths.dez

#Merge sexes
prev$age_group <- substring(prev$age_sex, 3, nchar(prev$age_sex))
prev <- prev %>% group_by(age_group) %>% summarise(ntests.dez = sum(ntests.dez),
                                                   npos.dez = sum(npos.dez),
                                                   population = sum(population),
                                                   npos.jan = sum(npos.jan),
                                                   ntests.jan = sum(ntests.jan),
                                                   npos.fev = sum(npos.fev),
                                                   ntests.fev = sum(ntests.fev),
                                                   npos.mar = sum(npos.mar),
                                                   ntests.mar = sum(ntests.mar),
                                                   deaths = sum(deaths))

Nsamples <- 50E3
prev <- prev %>% mutate(q1.a = 0, q2.a = 0, q3.a = 0, q4.a = 0, q5.a = 0,
                        IFR.q1.a = 0, IFR.q2.a = 0, IFR.q3.a = 0, IFR.q4.a = 0, IFR.q5.a = 0,
                        q1.b = 0, q2.b = 0, q3.b = 0, q4.b = 0, q5.b = 0,
                        IFR.q1.b = 0, IFR.q2.b = 0, IFR.q3.b = 0, IFR.q4.b = 0, IFR.q5.b = 0,
                        IFR.q1.a.RR = 0, IFR.q2.a.RR = 0, IFR.q3.a.RR = 0, IFR.q4.a.RR = 0, IFR.q5.a.RR = 0,
                        IFR.q1.b.RR = 0, IFR.q2.b.RR = 0, IFR.q3.b.RR = 0, IFR.q4.b.RR = 0, IFR.q5.b.RR = 0)
TN <- 820
TP <- 164
FN <- 34
FP <- 1

se <- rbeta(Nsamples, TP, FN) #sensitivity

modelCode <- nimbleCode({
  un[1:3] ~ ddirch(alpha[1:3])
  umax ~ dunif(0,1)
  for(i in 1:3){
    u[i] <- un[i]*umax
  }
  rho.dez ~ dunif(0,1)
  rho.jan <- rho.dez + u[1]
  rho.fev <- rho.dez + u[1] + u[2]
  rho.mar <- rho.dez + u[1] + u[2] + u[3]
  
  npos.dez ~ dbinom(rho.dez, ntests.dez)
  npos.jan ~ dbinom(rho.jan, ntests.jan)
  npos.fev ~ dbinom(rho.fev, ntests.fev)
  npos.mar ~ dbinom(rho.mar, ntests.mar)
})

prev2020.samples <- readRDS("IFR2020_samples.rds") %>% filter(city == "HEMOAM")
a.samples <- list()
b.samples <- list()

IFR.samples.a <- list()
IFR.samples.b <- list()
for(i in 1:nrow(prev))
{
  samples <- nimbleMCMC(code = modelCode, 
                        constants = list(ntests.dez = prev$ntests.dez[i],
                                         ntests.jan = prev$ntests.jan[i],
                                         ntests.fev = prev$ntests.fev[i],
                                         ntests.mar = prev$ntests.mar[i],
                                         alpha = rep(1, 3)), 
                        data = list(npos.dez = prev$npos.dez[i],
                                    npos.jan = prev$npos.jan[i],
                                    npos.fev = prev$npos.fev[i],
                                    npos.mar = prev$npos.mar[i]), 
                        inits = list(umax = 0.5, un = c(0.4, 0.3, 0.3), rho.dez = prev$npos.dez[i]/prev$ntests.dez[i]),
                        nburnin = Nsamples/2, niter = 3*Nsamples/2, monitors = c("u", "rho.dez"))
  
  
  b.samples[[i]] <- (samples[, "u[1]"] + samples[, "u[2]"] + samples[, "u[3]"])/se
  a.samples[[i]] <- samples[, "rho.dez"] + b.samples[[i]]
  
  
  prev$q1.a[i] <- quantile(a.samples[[i]], 0.025)
  prev$q2.a[i] <- quantile(a.samples[[i]], 0.25)
  prev$q3.a[i] <- quantile(a.samples[[i]], 0.5)
  prev$q4.a[i] <- quantile(a.samples[[i]], 0.75)
  prev$q5.a[i] <- quantile(a.samples[[i]], 0.975)
  
  prev$q1.b[i] <- quantile(b.samples[[i]], 0.025)
  prev$q2.b[i] <- quantile(b.samples[[i]], 0.25)
  prev$q3.b[i] <- quantile(b.samples[[i]], 0.5)
  prev$q4.b[i] <- quantile(b.samples[[i]], 0.75)
  prev$q5.b[i] <- quantile(b.samples[[i]], 0.975)
  
  # -- IFR -- #
  IFR.samples.a[[i]] <- sapply(1:Nsamples, function(k) rbeta(1, prev$deaths[i]+1, 1+max(0, floor(prev$population[i]*b.samples[[i]][k]) - prev$deaths[i])) )
  IFR.samples.b[[i]] <- sapply(1:Nsamples, function(k) rbeta(1, prev$deaths[i]+1, 1+max(0, floor(prev$population[i]*a.samples[[i]][k]) - prev$deaths[i])) )
  
  prev$IFR.q1.a[i] <- quantile(IFR.samples.a[[i]], 0.025)
  prev$IFR.q2.a[i] <- quantile(IFR.samples.a[[i]], 0.25)
  prev$IFR.q3.a[i] <- quantile(IFR.samples.a[[i]], 0.5)
  prev$IFR.q4.a[i] <- quantile(IFR.samples.a[[i]], 0.75)
  prev$IFR.q5.a[i] <- quantile(IFR.samples.a[[i]], 0.975)
  
  prev$IFR.q1.b[i] <- quantile(IFR.samples.b[[i]], 0.025)
  prev$IFR.q2.b[i] <- quantile(IFR.samples.b[[i]], 0.25)
  prev$IFR.q3.b[i] <- quantile(IFR.samples.b[[i]], 0.5)
  prev$IFR.q4.b[i] <- quantile(IFR.samples.b[[i]], 0.75)
  prev$IFR.q5.b[i] <- quantile(IFR.samples.b[[i]], 0.975)
  
  
  dff <- prev2020.samples %>% filter(age_bin == prev$age_group[i])
  IFR.samples.a.RR <- IFR.samples.a[[i]]/dff$IFR
  IFR.samples.b.RR <- IFR.samples.b[[i]]/dff$IFR
  
  prev$IFR.q1.a.RR[i] <- quantile(IFR.samples.a.RR, 0.025)
  prev$IFR.q2.a.RR[i] <- quantile(IFR.samples.a.RR, 0.25)
  prev$IFR.q3.a.RR[i] <- quantile(IFR.samples.a.RR, 0.5)
  prev$IFR.q4.a.RR[i] <- quantile(IFR.samples.a.RR, 0.75)
  prev$IFR.q5.a.RR[i] <- quantile(IFR.samples.a.RR, 0.975)
  
  prev$IFR.q1.b.RR[i] <- quantile(IFR.samples.b.RR, 0.025)
  prev$IFR.q2.b.RR[i] <- quantile(IFR.samples.b.RR, 0.25)
  prev$IFR.q3.b.RR[i] <- quantile(IFR.samples.b.RR, 0.5)
  prev$IFR.q4.b.RR[i] <- quantile(IFR.samples.b.RR, 0.75)
  prev$IFR.q5.b.RR[i] <- quantile(IFR.samples.b.RR, 0.975)
  
}

#Global IFR/Attack Rate
a.samples.global <- Reduce("+", lapply(1:5, function(i) a.samples[[i]]*prev$population[i]))/sum(prev$population)
b.samples.global <- Reduce("+", lapply(1:5, function(i) b.samples[[i]]*prev$population[i]))/sum(prev$population)

IFR.samples.a.global <- Reduce("+", lapply(1:5, function(i) IFR.samples.a[[i]]*a.samples.global[[i]]*prev$population[i]))/
  Reduce("+", lapply(1:5, function(i) a.samples.global[[i]]*prev$population[i]))
IFR.samples.b.global <- Reduce("+", lapply(1:5, function(i) IFR.samples.b[[i]]*b.samples.global[[i]]*prev$population[i]))/
  Reduce("+", lapply(1:5, function(i) b.samples.global[[i]]*prev$population[i]))

prev$q1.a.global <- quantile(a.samples.global, 0.025)
prev$q2.a.global <- quantile(a.samples.global, 0.25)
prev$q3.a.global <- quantile(a.samples.global, 0.5)
prev$q4.a.global <- quantile(a.samples.global, 0.75)
prev$q5.a.global <- quantile(a.samples.global, 0.975)

prev$q1.b.global <- quantile(b.samples.global, 0.025)
prev$q2.b.global <- quantile(b.samples.global, 0.25)
prev$q3.b.global <- quantile(b.samples.global, 0.5)
prev$q4.b.global <- quantile(b.samples.global, 0.75)
prev$q5.b.global <- quantile(b.samples.global, 0.975)

prev$IFR.q1.a.global <- quantile(IFR.samples.a.global, 0.025)
prev$IFR.q2.a.global <- quantile(IFR.samples.a.global, 0.25)
prev$IFR.q3.a.global <- quantile(IFR.samples.a.global, 0.5)
prev$IFR.q4.a.global <- quantile(IFR.samples.a.global, 0.75)
prev$IFR.q5.a.global <- quantile(IFR.samples.a.global, 0.975)

prev$IFR.q1.b.global <- quantile(IFR.samples.b.global, 0.025)
prev$IFR.q2.b.global <- quantile(IFR.samples.b.global, 0.25)
prev$IFR.q3.b.global <- quantile(IFR.samples.b.global, 0.5)
prev$IFR.q4.b.global <- quantile(IFR.samples.b.global, 0.75)
prev$IFR.q5.b.global <- quantile(IFR.samples.b.global, 0.975)

saveRDS(prev, file = "prev_Gamma_IFR.rds")

samples.global.2020 <- readRDS("IFR2020_global_samples.rds")
samples.global.2020 <- samples.global.2020 %>% filter(city == "HEMOAM")
RR.IFR.global <- IFR.samples.b.global/samples.global.2020$IFR
#How many times Gamma's IFR was higher
print(quantile(RR.IFR.global, c(0.025, 0.25, 0.5, 0.75, 0.975))) 


library(egg)
prev.2020 <- read.csv("IFR2020.csv") %>% rename(age_group = age_bin)
prev.2020 <- prev.2020 %>% filter(city == "HEMOAM")


cities <- c("HEMOAM", "FPS", "HEMOCE", "HEMOBA", "HEMOMINAS", "HEMOPE", "HEMEPAR", "HEMORIO")

prev.2021 <- read.csv("IFR_Manaus_2021.csv")  %>% rename(age_group = age_bin)


dfplot <- rbind(prev %>% mutate(type = "Upper bound", q1 = IFR.q1.a, q2 = IFR.q2.a, q3 = IFR.q3.a, q4 = IFR.q4.a, q5 = IFR.q5.a) %>%
                  select(type, age_group, q1, q2, q3, q4, q5),
                prev %>% mutate(type = "Lower bound", q1 = IFR.q1.b, q2 = IFR.q2.b, q3 = IFR.q3.b, q4 = IFR.q4.b, q5 = IFR.q5.b) %>% 
                  select(type, age_group, q1, q2, q3, q4, q5),
                prev.2020 %>% select(age_group, IFR.q1, IFR.q2, IFR.q3, IFR.q4, IFR.q5) %>% 
                  rename(q1 = IFR.q1, q2 = IFR.q2, q3 = IFR.q3, q4 = IFR.q4, q5 = IFR.q5) %>% mutate(type = "Estimated by model (2020)"),
                prev.2021 %>% select(age_group, IFR.q1, IFR.q2, IFR.q3, IFR.q4, IFR.q5) %>% 
                  rename(q1 = IFR.q1, q2 = IFR.q2, q3 = IFR.q3, q4 = IFR.q4, q5 = IFR.q5) %>% mutate(type = "Estimated by model (2021)"))

ggplot(dfplot %>% filter(type != "Upper bound"), 
       aes(x = age_group, fill = type)) + geom_errorbar(aes(ymin = 100*q1, ymax = 100*q5), position = "dodge2", width = 0.65, lwd = 0.25) + 
  geom_crossbar(aes(ymin = 100*q2, ymax = 100*q4, y = 100*q3), position = "dodge2", width = 0.65, lwd = 0.25) +
  scale_y_log10(n.breaks = 10)  + theme_bw() + labs(x = "Age/Sex", y = "IFR (%)") + scale_fill_manual("", values = c(mypalette[5], mypalette[6], mypalette[1:2])) +
  theme(legend.position = "top")
ggsave(set_panel_size(last_plot(), width = unit(10, 'cm'), height = unit(10, 'cm')), file = "figs/Gamma_IFR_mes.png", width = 6, height = 6) 
ggsave(set_panel_size(last_plot(), width = unit(10, 'cm'), height = unit(10, 'cm')), file = "figs/Gamma_IFR_mes.pdf", width = 6, height = 6) 

write.csv(dfplot, file = "IFR_gamma_ub_lb.csv")

# --- Incidence --- #


prev.2020 <- read.csv("prevalence_2020.csv") %>% filter(city == "HEMOAM")
dfplot <- rbind(prev %>% mutate(type = "Upper bound", q1 = q1.a, q2 = q2.a, q3 = q3.a, q4 = q4.a, q5 = q5.a) %>%
                  select(type, age_group, q1, q2, q3, q4, q5),
                prev %>% mutate(type = "Lower bound", q1 = q1.b, q2 = q2.b, q3 = q3.b, q4 = q4.b, q5 = q5.b) %>% 
                  select(type, age_group, q1, q2, q3, q4, q5),
                prev.2021 %>% mutate(type = "Estimated by model (2021)", q1 = cases.q1, q2 = cases.q2, q3 = cases.q3, 
                                     q4 = cases.q4, q5 = cases.q5) %>% select(type, age_group, q1, q2, q3, q4, q5),
                prev.2020 %>% mutate(type = "Estimated by model (2020)", age_group = age_bin) %>% select(type, age_group, q1, q2, q3, q4, q5))

ggplot(dfplot %>% filter(type != "Lower bound"), aes(x = age_group, fill = type)) + 
  geom_errorbar(aes(ymin = 100*q1, ymax = 100*q5), position = "dodge2", width = 0.65, lwd = 0.25) + 
  geom_crossbar(aes(ymin = 100*q2, ymax = 100*q4, y = 100*q3), position = "dodge2", width = 0.65, lwd = 0.25)  + theme_bw() + 
  labs(x = "Age group", y = "Incidence (%)") + scale_fill_manual("", values = c( mypalette[5:6], mypalette[1])) +
  theme(legend.position = "top") + scale_y_continuous(n.breaks = 10)
ggsave(set_panel_size(last_plot(), width = unit(10, 'cm'), height = unit(10, 'cm')), file = "figs/incidence_gamma.png", width = 6, height = 6) 
ggsave(set_panel_size(last_plot(), width = unit(10, 'cm'), height = unit(10, 'cm')), file = "figs/incidence_gamma2.pdf", width = 6, height = 6) 

# --- Relative Risk --- #

prev.2021.RR <- read.csv("IFR_Manaus_2021_RR.csv")
dfplot <- rbind(prev %>% mutate(type = "Upper bound", q1 = IFR.q1.a.RR, q2 = IFR.q2.a.RR, q3 = IFR.q3.a.RR, q4 = IFR.q4.a.RR, q5 = IFR.q5.a.RR) %>%
                  select(type, age_group, q1, q2, q3, q4, q5),
                prev %>% mutate(type = "Lower bound", q1 = IFR.q1.b.RR, q2 = IFR.q2.b.RR, q3 = IFR.q3.b.RR, q4 = IFR.q4.b.RR, q5 = IFR.q5.b.RR) %>% 
                  select(type, age_group, q1, q2, q3, q4, q5),
                prev.2021.RR %>% select(age_bin, IFR.RR.q1, IFR.RR.q2, IFR.RR.q3, IFR.RR.q4, IFR.RR.q5) %>% 
                  rename(age_group = age_bin, q1 = IFR.RR.q1, q2 = IFR.RR.q2, q3 = IFR.RR.q3, q4 = IFR.RR.q4, q5 = IFR.RR.q5) %>%
                  mutate(type = "Estimated by model (2021)"))

ggplot(dfplot %>% filter(type != "Upper bound"), aes(x = age_group, fill = type)) + 
  geom_errorbar(aes(ymin = q1, ymax = q5), position = "dodge2", width = 0.65, lwd = 0.25) + 
  geom_crossbar(aes(ymin = q2, ymax = q4, y = q3), position = "dodge2", width = 0.65, lwd = 0.25) +
  scale_y_log10(n.breaks = 10)  + theme_bw() + labs(x = "Age/Sex", y = "Relative risk using the IFR\n obtained for 2020 as reference") + 
  scale_fill_manual("", values = c(mypalette[6], mypalette[1:2])) + geom_hline(yintercept = 1, linetype = 'dashed') +
  theme(legend.position = "top")
ggsave(set_panel_size(last_plot(), width = unit(10, 'cm'), height = unit(10, 'cm')), file = "figs/Gamma_IFR_mes_RR.png", width = 6, height = 6) 
ggsave(set_panel_size(last_plot(), width = unit(10, 'cm'), height = unit(10, 'cm')), file = "figs/Gamma_IFR_mes_RR.pdf", width = 6, height = 6) 

write.csv(dfplot, file = "IFR_gamma_ub_lb_rr.csv")
