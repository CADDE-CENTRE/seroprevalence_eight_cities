library(tidyverse)
library(lubridate)
library(rstudioapi)
library(zoo)
setwd(dirname(getActiveDocumentContext()$path))
source("repeat_functions.R")

#This script estimates p+[n] from repeat blood donors and convalescent plasma donors data

df <- read.csv("data/repeat_blood_donors.csv")
df <- df %>% filter(covid_result != "", assay != "COV-2IgGII")
df$donation_date <- as.Date(paste0(df$donyr, "-", df$donmo, "-", df$donda))
df$result <- as.numeric(sapply(str_split(df$covid_result, pattern = " "), function(x) x[1]))

for(i in 1:nrow(df))
  df$result[i] <-  as.numeric(str_split(df$result[i], pattern = " ")[[1]][1])


thr <- 0.49
df <- compute.decay.rate(df)
df <- calculate.seroreversion.date(df, thr)

df <- df %>% group_by(donorid) %>% arrange(donation_date) %>% 
  mutate(kfirst = which(result >= thr)[1], 
         tinf.max = donation_date[kfirst],
         tinf.min = as.Date(ifelse(kfirst == 1, "2020-03-01", as.character(donation_date[kfirst - 1]))),
         ndon = n())

df$age_bin <- cut(df$age, c(-1, 25, 35, 45, 55, 70))


donors.all <- df %>% group_by(donorid) %>% 
  summarise(tinf.min = tinf.min[1], tinf.max = tinf.max[1])

df <- df %>% filter(kfirst < ndon, !is.na(serorev_date)) %>% group_by(donorid) %>% arrange(donation_date) %>%
  mutate(serorev_date_weibull = case_when(serorev_obs == FALSE ~ tail(donation_date, 1),
                                          serorev_obs == TRUE  ~ tail(donation_date[result >= thr], 1)),
         serorev_date_weibull_max = case_when(serorev_obs == FALSE ~ as.Date("2030-01-01"),
                                              serorev_obs == TRUE ~  donation_date[tail(which(result >= thr),1)+1]))

donors <- df %>% group_by(donorid) %>% 
  summarise(serorev_date = serorev_date[1], tinf.min = tinf.min[1], tinf.max = tinf.max[1],
            serorev_obs = serorev_obs[1], serorev_date_weibull = serorev_date_weibull[1],
            serorev_date_weibull_max = serorev_date_weibull_max[1])



dates <- seq(as.Date("2020-03-01"), as.Date("2021-03-31"), by = "day")
ptinf <- rep(0, length(dates))
for(i in 1:nrow(donors))
{
  ptinf <- ptinf + as.numeric(dates == donors.all$tinf.max[i])
}
ptinf <- ptinf/nrow(donors.all)
ptinf <- rollmean(ptinf, 30, align = "center")

srag <- read.csv("data/INFLUD20-14-02-2022.csv", sep = ";")
srag21 <- read.csv("data/INFLUD21-14-02-2022.csv", sep = ";")
srag <- rbind(srag, srag21[, colnames(srag)])
rm(srag21)

srag$CLASSI_FIN[is.na(srag$CLASSI_FIN)] <- 9
srag <- srag %>% filter(CLASSI_FIN %in% c(4,5,9), EVOLUCAO == 2, ID_MN_RESI == "MANAUS") %>% 
  mutate(DT_SIN_PRI = as.Date(DT_SIN_PRI, "%d/%m/%Y"))

ptinf_deaths <- rep(0, length(dates))
for(i in 1:length(dates))
{
  ptinf_deaths[i] <- sum(srag$DT_SIN_PRI == dates[i])
}
ptinf_deaths <- ptinf_deaths/sum(ptinf_deaths)
ptinf_deaths <- rollmean(ptinf_deaths, 7, align = "center")

Ns <- 1000
tts.exp <- matrix(data = 0, nrow = nrow(donors), ncol = Ns)
tts.deaths <- matrix(data = 0, nrow = nrow(donors), ncol = Ns)

for(i in 1:nrow(donors))
{
  datesi <- seq(1+donors$tinf.min[i], donors$tinf.max[i], by = "day") %>% as.character()
  i1 <- which(dates == 1+donors$tinf.min[i])[1]
  i2 <- which(dates == donors$tinf.max[i])[1]
  tinfi <- sample(datesi, size = Ns, replace = TRUE, prob = ptinf[i1:i2])
  tts.exp[i,] <- donors$serorev_date[i] - as.Date(tinfi)
  tinfi <- sample(datesi, size = Ns, replace = TRUE, prob = ptinf_deaths[i1:i2])
  tts.deaths[i,] <- donors$serorev_date[i] - as.Date(tinfi)
}

tts.exp <- as.vector(tts.exp)
tts.deaths <- as.vector(tts.deaths)

# --- Plasma donors --- #

dataset <- "AbbottRoche" #Same cohort used to estimate sensitivity
if(dataset == "AbbottRoche") {
plasma <- read.csv("data/convalescent_plasma_longitudinal_roche.csv")
hcw <- plasma %>% rename(id = donor_id, 
                         date_sample_collection = date_sample_collected, 
                         sc = abbott_sc, 
                         symp_start = date_symptom_onset,
                         hospit = hospitalized
                         )
} else{
  hcw <- read.csv("data/abbott_cohort_long.csv", stringsAsFactors = F) #Cohort used in Buss et al (not used in this paper)
}

hcw <- hcw %>% filter(hospit == 0) %>% mutate(donation_date = as.Date(date_sample_collection, format = "%m/%d/%y")) %>% rename(result = sc, donorid = id)
hcw$symp_start <- as.Date(hcw$symp_start, format = "%m/%d/%y")
hcw$date_sample_collection <- as.Date(hcw$date_sample_collection, format = "%m/%d/%y")
hcw <- hcw %>% mutate(days_from_symp = as.numeric(date_sample_collection - symp_start))

hcw <- hcw %>% group_by(donorid) %>% mutate(resmax = max(result), min_days_from_symp = min(days_from_symp), max_days_from_symp = min(days_from_symp))

hcw <- compute.decay.rate(hcw)
hcw <- calculate.seroreversion.date(hcw, thr)
hcw$tts_onset <- hcw$serorev_date - hcw$symp_start
tts.plasma <- hcw %>% group_by(donorid) %>% summarise(x = tts_onset[1]) %>% .$x
tts.plasma <- tts.plasma[!is.na(tts.plasma)]

tts.plasma <- tts.plasma - 8 #8 days is the mean time between onset and seroconversion
                             #ref: Orner et al - Comparison of SARS-CoV-2 IgM and IgG seroconversion profiles among 
                             #hospitalized patients in two US cities.

# --- Convert the observed times to seroreversion into a probability distribution --- #

ptts.exp <- rep(0, 106)
ptts.deaths <- rep(0, 106)
ptts.plasma <- rep(0, 106)
for(m in 1:106) {
  for(i in 1:7){
    for(j in (7*m + 1:7)) {
      if(j > i) {
        ptts.exp[m] <- ptts.exp[m] + sum(tts.exp == (j-i))
        ptts.deaths[m] <- ptts.deaths[m] + sum(tts.deaths == (j-i))
        ptts.plasma[m] <- ptts.plasma[m] + sum(tts.plasma == (j-i))
      }
    }
  }
}

ptts.exp <- ptts.exp/sum(ptts.exp)
ptts.deaths <- ptts.deaths/sum(ptts.deaths)
ptts.plasma <- ptts.plasma/sum(ptts.plasma)

ppos_rep <- c(1, 1 - cumsum(ptts.exp))
ppos_rep_deaths <- c(1, 1 - cumsum(ptts.deaths))
ppos_hcw_onset <- c(1, 1 - cumsum(ptts.plasma))

dfplot <- rbind(data.frame(date = seq(0, 720, 7), ppos = ppos_rep[1:103], type = "Repeat Donors"), 
                data.frame(date = seq(0, 720, 7), ppos = ppos_rep_deaths[1:103], type = "Repeat Donors (deaths)"),
                data.frame(date = seq(0, 720, 7), ppos = ppos_hcw_onset[1:103], type = "Plasma Donors"))

ggplot(dfplot, aes(x = date, y = ppos, color = type)) + geom_line(size = 1) + geom_point(size = 3) + theme_bw() + scale_color_discrete("") + 
  labs(x = "Days from first positive sample", y = "Probability of positive test") + 
  scale_x_continuous(breaks=seq(0, 720, 30)) + theme(legend.position = "top", legend.box = "horizontal", legend.direction = "vertical")
ggsave(last_plot(), file = "figs/ppos_week.png", width = 10, height = 8)
ggsave(last_plot(), file = "figs/ppos_week.pdf", width = 10, height = 8)
write.csv(dfplot, "data/ppos_week.csv")


dfplot <- rbind(data.frame(date = seq(0, 720, 7), pneg = c(0, -diff(ppos_rep[1:103])), type = "Repeat Donors"),
                data.frame(date = seq(0, 720, 7), pneg = c(0, -diff(ppos_rep_deaths[1:103])), type = "Repeat Donors (deaths)"), 
                data.frame(date = seq(0, 720, 7), pneg = c(0, -diff(ppos_hcw_onset[1:103])), type = "Plasma Donors"))

ggplot(dfplot, aes(x = date, y = pneg, color = type)) + geom_line(size = 1) + geom_point(size = 3) + theme_bw() + scale_color_discrete("") + 
  labs(x = "Days from reference date", y = "Probability of seroreversion given positivity at reference date") + 
  scale_x_continuous(breaks=seq(0, 720, 30))
ggsave(last_plot(), file = "figs/pneg_week.png", width = 10, height = 8)
ggsave(last_plot(), file = "figs/pneg_week.pdf", width = 10, height = 8)

write.csv(dfplot, "data/pneg_week.csv")

