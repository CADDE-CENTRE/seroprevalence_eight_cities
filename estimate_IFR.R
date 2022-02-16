library(tidyverse)
library(Hmisc)
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
source("Mamoeiro.R")

cities <- c("HEMOAM", "FPS", "HEMOCE", "HEMOBA", "HEMOMINAS", "HEMOPE", "HEMEPAR", "HEMORIO")
codes <- c("HEMOAM" = 130260, "FPS" = 355030, "HEMOCE" = 230440, "HEMOBA" = 292740, 
           "HEMOMINAS" = 310620, "HEMOPE" = 261160, "HEMEPAR" = 410690, "HEMORIO" = 330455)

fields <- paste0("drho[42, ", 1:10, "]")
fieldnames <- c("F_15-24", "F_25-34", "F_35-44", "F_45-54", "F_55-69",
                "M_15-24", "M_25-34", "M_35-44", "M_45-54", "M_55-69")

dfplot <- list()
i <- 1
for(city in cities) {
  dff_inc <- readRDS(paste0("outputs/samples_", city, "_repeat_week_agesex.rds"))
  dff <- 0*dff_inc
  for(i in 1:10) {
    dff[, paste0("drho[", 1, ", ", i, "]")] <- dff_inc[, paste0("drho[", 1, ", ", i, "]")]
    for(j in 2:42) {
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
deaths <- srag %>% filter(DT_SIN_PRI <= as.Date("2020-12-16")) %>% group_by(age_sex, CO_MUN_RES) %>% tally()

dfplot <- dfplot %>% left_join(deaths, by = c("age_sex", "CO_MUN_RES"))
dfplot <- dfplot %>% left_join(pop %>% rename(CO_MUN_RES = code), by = c("CO_MUN_RES", "age_sex"))

#dfplot <- dfplot %>% group_by(city, age_bin) %>% 
#  summarise(ninf = (value[CS_SEXO == "F"]*population[CS_SEXO == "F"] + value[CS_SEXO == "M"]*population[CS_SEXO == "M"]),
#            n = n[CS_SEXO == "F"] + n[CS_SEXO == "M"])

dfplot.global <- dfplot %>% group_by(city, sampleid) %>% summarise(ninf = sum(value*population), n = sum(n))
dfplot.global.age <- dfplot %>% group_by(age_bin, sampleid) %>% summarise(ninf = sum(value*population), n = sum(n))

dfplot.out <- dfplot %>% group_by(city, age_bin, sampleid) %>% summarise(ninf = sum(value*population)/sum(population))
dfplot.out <- dfplot.out %>% group_by(city, age_bin) %>% summarise(q1 = quantile(ninf, 0.025),
                                                               q2 = quantile(ninf, 0.25),
                                                               q3 = quantile(ninf, 0.5),
                                                               q4 = quantile(ninf, 0.75),
                                                               q5 = quantile(ninf, 0.975))
write.csv(dfplot.out, file = "prevalence_2020.csv")

dfplot <- dfplot %>% group_by(city, age_bin, sampleid) %>% summarise(ninf = sum(value*population), n = sum(n))



dfplot$IFR <- sapply(1:nrow(dfplot), function(x) rbeta(1, dfplot$n[x] + 1, floor(dfplot$ninf[x]) - dfplot$n[x] + 1))


poptotal <- sum(pop$population)
pop <- pop %>% group_by(age_sex) %>% mutate(population_all = sum(population))

icodes <- names(codes)
names(icodes) <- codes

dfplot.global.all <- dfplot %>% left_join(pop %>% mutate(city = icodes[as.character(code)]))
dfplot.global.all <- dfplot.global.all %>% group_by(city, sampleid) %>% summarise(IFR = sum(IFR*population_all)/poptotal)
dfplot.global.all <- dfplot.global.all %>% group_by(city) %>% dplyr::summarise(IFR.q1 = quantile(IFR, 0.025),
                                                                               IFR.q2 = quantile(IFR, 0.25),
                                                                               IFR.q3 = quantile(IFR, 0.5),
                                                                               IFR.q4 = quantile(IFR, 0.75),
                                                                               IFR.q5 = quantile(IFR, 0.975))


# ninftotal <- dfplot %>% filter(sampleid == 1) %>% .$ninf %>% sum()
# dfplot.global.all2 <- dfplot %>% left_join(pop %>% mutate(city = icodes[as.character(code)]))
# dfplot.global.all2 <- dfplot.global.all2 %>% group_by(age_bin, sampleid) %>% mutate(ninfall = sum(ninf))
# dfplot.global.all2 <- dfplot.global.all2 %>% group_by(city, sampleid) %>% summarise(IFR = sum(IFR*ninfall)/ninftotal)
# dfplot.global.all2 <- dfplot.global.all2 %>% group_by(city) %>% dplyr::summarise(IFR.q1 = quantile(IFR, 0.025),
#                                                                                  IFR.q2 = quantile(IFR, 0.25),
#                                                                                  IFR.q3 = quantile(IFR, 0.5),
#                                                                                  IFR.q4 = quantile(IFR, 0.75),
#                                                                                  IFR.q5 = quantile(IFR, 0.975))


saveRDS(dfplot, "IFR2020_samples.rds")


dfplot <- dfplot %>% group_by(age_bin, city) %>% dplyr::summarise(IFR.q1 = quantile(IFR, 0.025),
                                                                 IFR.q2 = quantile(IFR, 0.25),
                                                                 IFR.q3 = quantile(IFR, 0.5),
                                                                 IFR.q4 = quantile(IFR, 0.75),
                                                                 IFR.q5 = quantile(IFR, 0.975))


dfplot.global$IFR <- sapply(1:nrow(dfplot.global), function(x) rbeta(1, dfplot.global$n[x] + 1, floor(dfplot.global$ninf[x]) - dfplot.global$n[x] + 1))

saveRDS(dfplot.global, "IFR2020_global_samples.rds")
dfplot.global <- dfplot.global %>% group_by(city) %>% dplyr::summarise(IFR.q1 = quantile(IFR, 0.025),
                                                                  IFR.q2 = quantile(IFR, 0.25),
                                                                  IFR.q3 = quantile(IFR, 0.5),
                                                                  IFR.q4 = quantile(IFR, 0.75),
                                                                  IFR.q5 = quantile(IFR, 0.975))

dfplot.global.age$IFR <- sapply(1:nrow(dfplot.global.age), function(x) rbeta(1, dfplot.global.age$n[x] + 1, 
                                                                             floor(dfplot.global.age$ninf[x]) - dfplot.global.age$n[x] + 1))

dfplot.global.age <- dfplot.global.age %>% group_by(age_bin) %>% dplyr::summarise(IFR.q1 = quantile(IFR, 0.025),
                                                                                 IFR.q2 = quantile(IFR, 0.25),
                                                                                 IFR.q3 = quantile(IFR, 0.5),
                                                                                 IFR.q4 = quantile(IFR, 0.75),
                                                                                 IFR.q5 = quantile(IFR, 0.975))



# ----------- #
dfplot <- dfplot %>% mutate(location = case_when(
  city == "FPS" ~ "São Paulo",
  city == "HEMEPAR" ~ "Curitiba",
  city == "HEMOAM" ~ "Manaus",
  city == "HEMOBA" ~ "Salvador",
  city == "HEMOCE" ~ "Fortaleza",
  city == "HEMOMINAS" ~ "Belo Horizonte",
  city == "HEMOPE" ~ "Recife",
  city == "HEMORIO" ~ "Rio de Janeiro"))


ggplot(dfplot, aes(x = age_bin, y = 100*IFR.q3, fill = location)) + geom_errorbar(aes(ymin = 100*IFR.q1, ymax = 100*IFR.q5), position = "dodge2", lwd = 0.25) + 
  geom_crossbar(aes(ymin =  100*IFR.q2, ymax =  100*IFR.q4), position = "dodge2", lwd = 0.25) + theme_bw() + labs(x = "Age group", y = "IFR (%)") +
  scale_y_log10() + scale_fill_manual("", values = mypalette)
ggsave(last_plot(), file = "figs/IFR_eightcities.png", width = 10, height = 5, dpi = 600)
ggsave(last_plot(), file = "figs/IFR_eightcities.pdf", width = 10, height = 5, dpi = 600)

write.csv(dfplot, "IFR2020.csv")

# --- Global IFRs --- #

dfplot.global <- dfplot.global %>% mutate(location = case_when(
  city == "FPS" ~ "São Paulo",
  city == "HEMEPAR" ~ "Curitiba",
  city == "HEMOAM" ~ "Manaus",
  city == "HEMOBA" ~ "Salvador",
  city == "HEMOCE" ~ "Fortaleza",
  city == "HEMOMINAS" ~ "Belo Horizonte",
  city == "HEMOPE" ~ "Recife",
  city == "HEMORIO" ~ "Rio de Janeiro"))

dfplot.global.all <- dfplot.global.all %>% mutate(location = case_when(
  city == "FPS" ~ "São Paulo",
  city == "HEMEPAR" ~ "Curitiba",
  city == "HEMOAM" ~ "Manaus",
  city == "HEMOBA" ~ "Salvador",
  city == "HEMOCE" ~ "Fortaleza",
  city == "HEMOMINAS" ~ "Belo Horizonte",
  city == "HEMOPE" ~ "Recife",
  city == "HEMORIO" ~ "Rio de Janeiro"))

dfplot.global <- rbind(dfplot.global %>% mutate(type = "(A) Crude global IFR"),
                       dfplot.global.all %>% mutate(type = "(B) Age-adjusted global IFR"))

#ggplot(dfplot.global, aes(x = location, y = 100*IFR.q3, fill = location)) + geom_errorbar(aes(ymin = 100*IFR.q1, ymax = 100*IFR.q5)) + 
#  geom_crossbar(aes(ymin = 100*IFR.q2, ymax = 100*IFR.q4)) + scale_fill_manual("", values = mypalette) + theme_bw() +
#  labs(x = "", y = "Global IFR (%)") + theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggplot(dfplot.global, aes(x = type, y = 100*IFR.q3, fill = location)) + geom_errorbar(aes(ymin = 100*IFR.q1, ymax = 100*IFR.q5), 
                                                                                      position = "dodge2", lwd = 0.25) + 
  geom_crossbar(aes(ymin = 100*IFR.q2, ymax = 100*IFR.q4), position = "dodge2", lwd = 0.25) + scale_fill_manual("", values = mypalette) + theme_bw() +
  labs(x = "", y = "Global IFR (%)") + theme(legend.position = "none")


ggsave(last_plot(), file = "figs/IFR_eightcities_global.png", width = 5, height = 5, dpi = 600)
ggsave(last_plot(), file = "figs/IFR_eightcities_global.pdf", width = 5, height = 5, dpi = 600)
write.csv(dfplot.global, "IFR2020global.csv")

# --- Compare with literature --- #

library(readxl)
lit <- read_excel("data/IFR.xlsx")
IFR.literature <- rbind(
  data.frame(age_bin = lit$age_bin, IFR.q1 = lit$Odriscoll_cil, IFR.q5 = lit$Odriscoll_ciu, IFR.q3 = lit$Odriscoll, type = "O'Driscoll"),
  data.frame(age_bin = lit$age_bin, IFR.q1 = lit$Brazeau_cil, IFR.q5 = lit$Brazeau_ciu, IFR.q3 = lit$Brazeau, type = "Brazeau")
)


dfplot.all <- dfplot.global.age %>% select(-c(IFR.q2, IFR.q4)) %>% mutate(type = "Our estimate")

pop.all <- readRDS("data/population2020.rds") %>% filter(floor(code/10) %in% codes)
pop.all <- pop.all %>% group_by(age_bin) %>% summarise(population = sum(population)) %>% filter(age_bin != "90+")
pop.all$age_bin <- sapply(1:nrow(pop.all), function(x) paste0(str_split(pop.all$age_bin[x], "-")[[1]][1], "-", 
                                                              as.numeric(str_split(pop.all$age_bin[x], "-")[[1]][2])-1)) 
IFR.literature <- IFR.literature %>% left_join(pop.all, by = "age_bin")

agemap <- c("15-19" = "15-24", "20-24" = "15-24", 
            "25-29" = "25-34", "30-34" = "25-34", "35-39" = "35-44", "40-44" = "35-44",
            "45-49" = "45-54", "50-54" = "45-54", "55-59" = "55-64", "60-64" = "55-64")
IFR.literature <- IFR.literature %>% mutate(age_bin = agemap[age_bin])
IFR.literature <- IFR.literature %>% group_by(age_bin, type) %>% summarise(IFR.q1 = sum(0.01*IFR.q1*population)/sum(population),
                                                                           IFR.q3 = sum(0.01*IFR.q3*population)/sum(population),
                                                                           IFR.q5 = sum(0.01*IFR.q5*population)/sum(population))



dfplot.all <- rbind(dfplot.all %>% mutate(type = "Our estimate"), IFR.literature)
dfplot.all$G <- 1:nrow(dfplot.all)

write.csv(dfplot.all, file = "IFR_global_literature.csv")

ggplot(dfplot.all, aes(x = age_bin, y = 100*IFR.q3, colour = type, group = G)) +
  geom_errorbar(aes(ymin = 100*IFR.q1, ymax =  100*IFR.q5), position = position_dodge(width = 0.75)) + 
 theme_bw() + labs(x = "Age group", y = "IFR (%)") + geom_point(position = position_dodge(width = 0.75)) +
  scale_y_log10() + scale_colour_manual("", values = mypalette) + theme(legend.position = "top")
ggsave(last_plot(), file = "figs/IFR_global_literature.png", width = 5, height = 5)
ggsave(last_plot(), file = "figs/IFR_global_literature.pdf", width = 5, height = 5)
