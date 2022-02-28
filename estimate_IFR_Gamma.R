library(tidyverse)
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

#This script generates bounds for the attack rate and IFR of the Gamma epidemic in Manaus
cities <- c("HEMOAM")
codes <- c("HEMOAM" = 130260, "FPS" = 355030, "HEMOCE" = 230440, "HEMOBA" = 292740, 
           "HEMOMINAS" = 310620, "HEMOPE" = 261160, "HEMEPAR" = 410690, "HEMORIO" = 330455)

fields <- paste0("drho[54, ", 1:10, "]") #drho[1] = week 10, hence drho[54] = week 63

fieldnames <- c("F_15-24", "F_25-34", "F_35-44", "F_45-54", "F_55-69",
                "M_15-24", "M_25-34", "M_35-44", "M_45-54", "M_55-69")

dfplot <- list()
i <- 1
for(city in cities) {
  dff_inc <- readRDS(paste0("outputs/samples_", city, "_repeat_week_agesex.rds"))
  dff <- 0*dff_inc
  for(i in 1:10) {
    dff[, paste0("drho[", 42, ", ", i, "]")] <- dff_inc[, paste0("drho[", 1, ", ", i, "]")] #drho[1] = week 10, hence drho[42] = week 51
    for(j in 43:54) {
      dff[, paste0("drho[", j, ", ", i, "]")] <- dff[, paste0("drho[", j - 1, ", ", i, "]")] + dff_inc[, paste0("drho[", j, ", ", i, "]")]
    }
  }
  
  dfplot[[city]] <- dff %>% as.data.frame() %>% select(fields) %>% mutate(city = city)
}

dfplot <- do.call(rbind, dfplot)
colnames(dfplot) <- c(fieldnames, "city")
dfplot <- dfplot %>% gather("age_sex", "value", -"city")

dfplot <- dfplot %>% mutate(CO_MUN_RES = codes[city])
dfplot$age_sex[dfplot$age_sex == "F_55-69"] <- "F_55-64"
dfplot$age_sex[dfplot$age_sex == "M_55-69"] <- "M_55-64"

dfplot <- dfplot %>% group_by(city, age_sex) %>% mutate(sampleid = 1:n())

# --- SRAG --- #

pop <- readRDS("data/population2020.rds")
pop <- pop %>% mutate(code = floor(code/10)) %>% 
  filter(code %in% codes, age_bin %in% c("15-20", "20-25", "25-30", "30-35", "35-40", "40-45", "45-50", "50-55", "55-60", "60-65"))
agemap <- c("15-20" = "15-24", "20-25" = "15-24", "25-30" = "25-34", "30-35" = "25-34", 
            "35-40" = "35-44", "40-45" = "35-44", "45-50" = "45-54", "50-55" = "45-54", "55-60" = "55-64", "60-65" = "55-64")
pop <- pop %>% mutate(new_bin = agemap[age_bin])
pop <- pop %>% group_by(code, new_bin, CS_SEXO) %>% summarise(population = sum(population)) %>% 
  ungroup() %>% rename(age_bin = new_bin) %>% mutate(age_sex = paste0(CS_SEXO, "_", age_bin)) %>% arrange(age_sex)
#pop$age_sex[pop$age_sex == "F_55-69"] <- "F_55-64"
#pop$age_sex[pop$age_sex == "M_55-69"] <- "M_55-64"


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
deaths <- srag %>% filter(DT_SIN_PRI >= as.Date("2020-12-16"), DT_SIN_PRI <= as.Date("2021-03-10")) %>% group_by(age_sex, CO_MUN_RES) %>% tally()

dfplot <- dfplot %>% left_join(deaths, by = c("age_sex", "CO_MUN_RES"))
dfplot <- dfplot %>% left_join(pop %>% rename(CO_MUN_RES = code), by = c("CO_MUN_RES", "age_sex"))


dfplot <- dfplot %>% group_by(city, age_bin, sampleid) %>% summarise(ninf = sum(value*population), n = sum(n), 
                                                                     u = ninf/sum(population))
dfplot$IFR <- sapply(1:nrow(dfplot), function(x) rbeta(1, dfplot$n[x] + 1, floor(dfplot$ninf[x]) - dfplot$n[x] + 1))

dfplot.samples <- dfplot #Save to compute the odds ratio

dfplot <- dfplot %>% group_by(age_bin) %>% dplyr::summarise(IFR.q1 = quantile(IFR, 0.025),
                                                            IFR.q2 = quantile(IFR, 0.25),
                                                            IFR.q3 = quantile(IFR, 0.5),
                                                            IFR.q4 = quantile(IFR, 0.75),
                                                            IFR.q5 = quantile(IFR, 0.975),
                                                            cases.q1 = quantile(u, 0.025),
                                                            cases.q2 = quantile(u, 0.25),
                                                            cases.q3 = quantile(u, 0.5),
                                                            cases.q4 = quantile(u, 0.75),
                                                            cases.q5 = quantile(u, 0.975))


# ------ Generate plots ----- #
dfplot <- dfplot %>% mutate(location = case_when(
  city == "FPS" ~ "São Paulo",
  city == "HEMEPAR" ~ "Curitiba",
  city == "HEMOAM" ~ "Manaus",
  city == "HEMOBA" ~ "Salvador",
  city == "HEMOCE" ~ "Fortaleza",
  city == "HEMOMINAS" ~ "Belo Horizonte",
  city == "HEMOPE" ~ "Recife",
  city == "HEMORIO" ~ "Rio de Janeiro"))

source("Mamoeiro.R")

ggplot(dfplot, aes(x = age_bin, y = 100*IFR.q3, fill = location)) + geom_errorbar(aes(ymin = 100*IFR.q1, ymax =  100*IFR.q5), position = "dodge2") + 
  geom_crossbar(aes(ymin =  100*IFR.q2, ymax =  100*IFR.q4), position = "dodge2") + theme_bw() + labs(x = "Age group", y = "IFR (%)") +
  scale_y_log10() + scale_fill_manual("", values = mypalette)
ggsave(last_plot(), file = "figs/IFR_Manaus_2021.png", width = 10, height = 5)

ggplot(dfplot, aes(x = age_bin, y = 100*cases.q3, fill = location)) + geom_errorbar(aes(ymin = 100*cases.q1, ymax =  100*cases.q5), position = "dodge2") + 
  geom_crossbar(aes(ymin =  100*cases.q2, ymax =  100*cases.q4), position = "dodge2") + theme_bw() + labs(x = "Age group", y = "Incidence (%)") +
  scale_y_log10() + scale_fill_manual("", values = mypalette)

write.csv(dfplot, "IFR_Manaus_2021.csv")


# --- Odds Ratio using 2020 as reference --- #

dfplot.samples.2020 <- readRDS("IFR2020_samples.rds") %>% filter(city == "HEMOAM") %>% ungroup()
dfplot.samples.2020 <- dfplot.samples.2020 %>% select(age_bin, sampleid, IFR) %>% rename(IFR.2020 = IFR)
dfplot.samples <- dfplot.samples %>% left_join(dfplot.samples.2020, by = c( "age_bin", "sampleid"))
dfplot.samples <- dfplot.samples %>% mutate(IFR.RR = IFR/IFR.2020)

dfplot <- dfplot.samples %>% group_by(age_bin) %>% summarise(IFR.RR.q1 = quantile(IFR.RR, 0.025),
                                                             IFR.RR.q2 = quantile(IFR.RR, 0.25),
                                                             IFR.RR.q3 = quantile(IFR.RR, 0.5),
                                                             IFR.RR.q4 = quantile(IFR.RR, 0.75),
                                                             IFR.RR.q5 = quantile(IFR.RR, 0.975))
write.csv(dfplot, "IFR_Manaus_2021_RR.csv")
