library(tidyverse)
library(zoo)
source("Mamoeiro.R")

dates <- as.Date(c("2020-04-15", "2020-05-15", "2020-06-15", "2020-07-15", "2020-08-15", "2020-09-15", 
                   "2020-10-15", "2020-11-15", "2020-12-15"))


pop <- readRDS("data/population2020.rds")
pop <- pop %>% filter(age_bin %in% c("15-20", "20-25", "25-30", "30-35", "35-40", "40-45", "45-50", "50-55", "55-60", "60-65", "65-70"))
agemap <- c("15-20" = "15-24", "20-25" = "15-24", "25-30" = "25-34", "30-35" = "25-34", 
            "35-40" = "35-44", "40-45" = "35-44", "45-50" = "45-54", "50-55" = "45-54", "55-60" = "55-70", "60-65" = "55-70", "65-70" = "55-70")
pop <- pop %>% mutate(new_bin = agemap[age_bin])
pop <- pop %>% group_by(code, new_bin, CS_SEXO) %>% summarise(population = sum(population)) %>% 
  ungroup() %>% rename(age_bin = new_bin) %>% mutate(age_sex = paste0(CS_SEXO, "_", age_bin)) %>% arrange(age_sex)

# --- Load SIVEP-Gripe --- #
srag <- read.csv("data/INFLUD20-14-02-2022.csv", stringsAsFactors = F, sep = ";")
srag21 <- read.csv("data/INFLUD21-14-02-2022.csv", stringsAsFactors = F, sep = ";")
srag21 <- srag21[, colnames(srag)]
srag <- rbind(srag, srag21)

srag <- srag %>% filter(CO_MUN_RES %in% c(355030, 410690, 230440, 310620, 261160, 330455, 130260, 292740), 
                        EVOLUCAO == 2) %>%
  mutate(DT_EVOLUCA = as.Date(DT_EVOLUCA, format = "%d/%m/%Y"),
         DT_SIN_PRI = as.Date(DT_SIN_PRI, format = "%d/%m/%Y"),
         age = as.numeric(ifelse(TP_IDADE == 3, NU_IDADE_N, 0)))
srag$CLASSI_FIN[is.na(srag$CLASSI_FIN)] <- 9
srag <- srag %>% filter(CLASSI_FIN %in% c(4,5,9), DT_SIN_PRI >= as.Date("2020-02-26"))

agebins <- c("15-24", "25-34", "35-44", "45-54", "55-70")
srag$age_group <- cut(srag$age, c(15, 25, 35, 45, 55, 70), labels = agebins)

srag <- srag %>% filter(!is.na(age_group), CS_SEXO %in% c("M","F")) %>% mutate(age_sex = paste0(CS_SEXO, "_", age_group))

srag <- srag %>% group_by(CO_MUN_RES, age_sex, DT_SIN_PRI) %>% summarise(deaths = n(), ID_MN_RESI = ID_MN_RESI[1])

pop <- pop %>% mutate(code = floor(code/10))

srag <- srag %>% left_join(pop %>% select(code, age_sex, population) %>% rename(CO_MUN_RES = code), by = c("CO_MUN_RES", "age_sex"))
srag <- srag %>% mutate(SMR = deaths/population)

popall <- pop %>% group_by(age_sex) %>% summarise(population_all = sum(population))
popcity <- pop %>% group_by(code) %>% summarise(popcity = sum(population))

srag <- srag %>% left_join(popall, by = "age_sex")
srag <- srag %>% left_join(popcity %>% rename(CO_MUN_RES = code), by = "CO_MUN_RES")

pop_total = popall$population_all %>% sum()

# --- Age-Standardised Number of Deaths --- #
sragn <- srag %>% group_by(CO_MUN_RES, DT_SIN_PRI) %>% summarise(SMR = sum(SMR*population_all/pop_total), ID_MN_RESI = ID_MN_RESI[1],
                                                                 MR = sum(deaths)/popcity[1])

sragn <- sragn %>% group_by(CO_MUN_RES) %>% arrange(DT_SIN_PRI) %>% mutate(SMR_mean = rollmean(SMR, 7, fill = 0, align = "center"), 
                                                                           SMR_cum = cumsum(SMR_mean),
                                                                           MR_mean = rollmean(MR, 7, fill = 0, align = "center"), 
                                                                           MR_cum = cumsum(MR_mean))


# --- Generate plots --- #
ggplot(sragn %>% filter(DT_SIN_PRI <= as.Date("2021-03-31")), aes(x = DT_SIN_PRI, y = 1E6*SMR, colour = ID_MN_RESI)) + facet_wrap(~ID_MN_RESI, nrow = 2, ncol = 4) + 
  geom_point(alpha = 0.25, size = 1) + geom_line(aes(y = 1E6*SMR_mean)) + theme_bw() + 
  labs(x = "Date of symptoms onset", y = "Age-sex adjusted daily mortality per million inhabitants") +
  theme(legend.position = "none") + scale_x_continuous(breaks = seq(as.Date("2020-01-01"), as.Date("2021-03-01"), by = "month"), 
                                                       labels = c("Jan", "Fev", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", 
                                                                  "Jan", "Fev", "Mar")) +
  ylim(c(0, 90)) + geom_line(aes(x = DT_SIN_PRI, y = 1E6*MR_mean), linetype = 2)
ggsave(last_plot(), file = "figs/age-standardised_MR.png", width= 15, height = 10)


citymap <- c("BELO HORIZONTE" = "Belo Horizonte", "CURITIBA" = "Curitiba", "FORTALEZA" = "Fortaleza", "MANAUS" = "Manaus", "RECIFE" = "Recife",
             "RIO DE JANEIRO" = "Rio de Janeiro", "SALVADOR" = "Salvador", "SAO PAULO" = "São Paulo")
sragn <- sragn %>% mutate(city = citymap[ID_MN_RESI])
ggplot(sragn %>% filter(DT_SIN_PRI <= as.Date("2021-03-31")), aes(x = DT_SIN_PRI, y = 1E6*MR, colour = city)) + facet_wrap(~city, nrow = 4, ncol = 2) + 
  geom_point(alpha = 0.25, size = 1) + geom_line(aes(x = DT_SIN_PRI, y = 1E6*MR_mean, group = city)) + theme_bw() + 
  labs(x = "Date of symptoms onset", y = "Mortality per million inhabitants") +
  theme(legend.position = "none" , strip.background =element_rect(fill="#f1e1ab")) + scale_x_continuous(breaks = seq(as.Date("2020-01-01"), as.Date("2021-03-01"), by = "month"), 
                                                       labels = c("Jan", "Fev", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", 
                                                                  "Jan", "Fev", "Mar")) +
  ylim(c(0, 75)) + scale_colour_manual("", values = mypalette)
ggsave(last_plot(), file = "figs/crude_MR.png", width= 10, height = 10)
