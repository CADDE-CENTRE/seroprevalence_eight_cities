library(tidyverse)
library(Hmisc)
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

#This script validates the attack rates and IFRs obtained by our seroreversion correction model
#by estimating these quantities using a threshold of 0.1 and correcting for sensitivity and
#specificity

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
  
  thr <- 0.1
  df <- all_centres %>% filter(blood_center == center, month <= 12)
  prev <- df %>% group_by(month, age_sex) %>% summarise(ntests = n(), npos = sum(result >= thr))
  prev <- expand.grid(month = min(prev$month):max(prev$month), age_sex = unique(prev$age_sex)) %>% left_join(prev, by = c("month", "age_sex"))
  prev$ntests[is.na(prev$ntests)] <- 0
  prev$npos[is.na(prev$npos)] <- 0
  prev <- prev %>% arrange(month, age_sex)
  L[[i]] <- prev %>% mutate(center = center)
}
prev <- do.call(rbind, L)
prev <- prev %>% filter(month == 12)


pop <- readRDS("data/population2020.rds")
pop <- pop %>% mutate(code = floor(code/10)) %>% 
  filter(code %in% codes, age_bin %in% c("15-20", "20-25", "25-30", "30-35", "35-40", "40-45", "45-50", "50-55", "55-60", "60-65"))
agemap <- c("15-20" = "15-24", "20-25" = "15-24", "25-30" = "25-34", "30-35" = "25-34", 
            "35-40" = "35-44", "40-45" = "35-44", "45-50" = "45-54", "50-55" = "45-54", "55-60" = "55-64", "60-65" = "55-64")
pop <- pop %>% mutate(new_bin = agemap[age_bin])
pop <- pop %>% group_by(code, new_bin, CS_SEXO) %>% summarise(population = sum(population)) %>% 
  ungroup() %>% rename(age_bin = new_bin) %>% mutate(age_sex = paste0(CS_SEXO, "_", age_bin)) %>% arrange(age_sex)

#SARI data
srag <- read.csv("data/INFLUD20-14-02-2022.csv", stringsAsFactors = F, sep = ";")
srag21 <- read.csv("data/INFLUD21-14-02-2022.csv", stringsAsFactors = F, sep = ";")
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
deaths <- srag %>% filter(DT_SIN_PRI <= as.Date("2020-12-15")) %>% group_by(age_sex, CO_MUN_RES) %>% tally()

prev <- prev %>% mutate(CO_MUN_RES = codes[center]) %>% left_join(deaths %>% rename(deaths = n), by = c("CO_MUN_RES", "age_sex"))
prev <- prev %>% left_join(pop %>% rename(CO_MUN_RES = code) %>% select(-c(age_bin, CS_SEXO)), by = c("CO_MUN_RES", "age_sex"))

Nsamples <- 1E5
prev <- prev %>% mutate(q1 = 0, q2 = 0, q3 = 0, q4 = 0, q5 = 0,
                        IFR.q1 = 0, IFR.q2 = 0, IFR.q3 = 0, IFR.q4 = 0, IFR.q5 = 0,
                        q1.global = 0, q2.global = 0, q3.global = 0, q4.global = 0, q5.global = 0,
                        age = substring(age_sex, 3, nchar(age_sex)), 
                        sex =  substring(age_sex, 1, 1))
TN <- 709
TP <- 197
FN <- 11
FP <- 112
for(center in unique(prev$center)) {
  prev.samples <- list()
  for(age in unique(prev$age)) {
      for(sex in c("M", "F")) {
        agesex <- paste0(sex, "_", age)
        
        i <- which((prev$age == age) & (prev$sex == sex) & (prev$center == center))
        se <- rbeta(Nsamples, TP+1, FN+1)
        sp <- rbeta(Nsamples, TN+1, FP+1)
        Z <- sapply(1:Nsamples, function(k) rbeta(1, prev$npos[i]+1, prev$ntests[i] - prev$npos[i]+1))
        prev.samples[[agesex]] <- (Z + sp - 1)/(se+sp-1)
        prev.samples[[agesex]][prev.samples[[agesex]] < 0] <- 0
        
        prev$q1[i] <- quantile(prev.samples[[agesex]], 0.025)
        prev$q2[i] <- quantile(prev.samples[[agesex]], 0.25)
        prev$q3[i] <- quantile(prev.samples[[agesex]], 0.5)
        prev$q4[i] <- quantile(prev.samples[[agesex]], 0.75)
        prev$q5[i] <- quantile(prev.samples[[agesex]], 0.975)
      }
        
      iM <- which((prev$age == age) & (prev$sex == "M") & (prev$center == center))
      iF <- which((prev$age == age) & (prev$sex == "F") & (prev$center == center))
      
      # -- IFR -- #
      IFR.samples <- sapply(1:Nsamples, 
                            function(k) rbeta(1, prev$deaths[iM] + prev$deaths[iF] + 1, 
                            1+max(0, floor(prev$population[iM]*prev.samples[[paste0("M", "_",  age)]][k] + 
                                             prev$population[iF]*prev.samples[[paste0("F", "_", age)]][k] - prev$deaths[iM] - prev$deaths[iF])) ))
      
      prev$IFR.q1[iM] <- quantile(IFR.samples, 0.025)
      prev$IFR.q2[iM] <- quantile(IFR.samples, 0.25)
      prev$IFR.q3[iM] <- quantile(IFR.samples, 0.5)
      prev$IFR.q4[iM] <- quantile(IFR.samples, 0.75)
      prev$IFR.q5[iM] <- quantile(IFR.samples, 0.975)
      
      prev$IFR.q1[iF] <- prev$IFR.q1[iM]
      prev$IFR.q2[iF] <- prev$IFR.q2[iM]
      prev$IFR.q3[iF] <- prev$IFR.q3[iM]
      prev$IFR.q4[iF] <- prev$IFR.q4[iM]
      prev$IFR.q5[iF] <- prev$IFR.q5[iM]
  }
  
  prev.samples.global <- rep(0, Nsamples)
  w.aux <- prev[prev$center == center, ] %>% select(age_sex, population) %>% arrange(age_sex)
  w <- w.aux$population/sum(w.aux$population)
  names(w) <- w.aux$age_sex
  
  for(a in unique(prev$age_sex)) {
    prev.samples.global <- prev.samples.global + w[a]*prev.samples[[a]]
  }
  
  I <- which(prev$center == center)
  prev$q1.global[I] <- quantile(prev.samples.global, 0.025)
  prev$q2.global[I] <- quantile(prev.samples.global, 0.25)
  prev$q3.global[I] <- quantile(prev.samples.global, 0.5)
  prev$q4.global[I] <- quantile(prev.samples.global, 0.75)
  prev$q5.global[I] <- quantile(prev.samples.global, 0.975)
  
}

# ----- Plots ------ #

citynames <- c("HEMOAM" = "Manaus", "FPS" = "São Paulo", "HEMOCE" = "Fortaleza", "HEMOBA" = "Salvador", 
               "HEMOMINAS" = "Belo Horizonte", "HEMOPE" = "Recife", "HEMEPAR" = "Curitiba", "HEMORIO" = "Rio de Janeiro")
prev <- prev %>% mutate(city = citynames[center])

source("Mamoeiro.R")
ggplot(prev %>% filter(sex == "M"), aes(x = age, y = 100*IFR.q3, fill = city)) + geom_errorbar(aes(ymin = 100*IFR.q1, ymax =  100*IFR.q5), position = "dodge2") + 
  geom_crossbar(aes(ymin =  100*IFR.q2, ymax =  100*IFR.q4), position = "dodge2") + theme_bw() + labs(x = "Age", y = "IFR (%)") +
  scale_y_log10() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual("", values = mypalette)
ggsave(last_plot(), file = "figs/IFR_measured_01.png", width = 10, height = 5)
ggsave(last_plot(), file = "figs/IFR_measured_01.pdf", width = 10, height = 5)


ggplot(prev %>% filter(sex == "M"), aes(x = city, y = 100*q3, fill = age)) + geom_errorbar(aes(ymin = 100*q1, ymax =  100*q5, width = 0.5), position = "dodge2") + 
  geom_crossbar(aes(ymin =  100*q2, ymax =  100*q4), position = "dodge2", width = 0.5) + theme_bw() + labs(x = "", y = "Prevalence (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  scale_fill_manual("", values = mypalette)
ggsave(last_plot(), file = "figs/prev_measured_01.png", width = 10, height = 5)
ggsave(last_plot(), file = "figs/prev_measured_01.pdf", width = 10, height = 5)

write.csv(prev %>% filter(sex == "M"), "IFR2020_measured.csv")

