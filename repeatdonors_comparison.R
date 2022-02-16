library(tidyverse)
library(cowplot)
library(rstudioapi)
source("Mamoeiro.R")

setwd(dirname(getActiveDocumentContext()$path))
#This script compares repeat blood donors and convalescent plasma donors
#It was used to generate Figure 1.

thr <- 0.49
df <- read.csv("data/repeat_blood_donors.csv")

df$donation_date <- as.Date(paste0(df$donyr, "-", df$donmo, "-", df$donda))
df <- df %>% filter(covid_result != "", assay != "COV-2IgGII")
df$result <- as.numeric(sapply(str_split(df$covid_result, pattern = " "), function(x) x[1]))
df <- df %>% group_by(donorid) %>% arrange(donation_date) %>%
  mutate(next_result = lead(result), next_date = lead(donation_date), halflife = -log(2)*as.numeric(next_date - donation_date)/(log(next_result) - log(result)))

resmax_repeat <- df %>% filter(result >= thr, donation_date <= as.Date("2020-05-31")) %>% group_by(donorid) %>% summarise(result_max = max(result)) %>% .$result_max
HL_repeat <- df %>% filter(result >= thr, !is.na(halflife), halflife > 0) %>% .$halflife

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
  hcw <- read.csv("data/abbott_cohort_long.csv", stringsAsFactors = F)
}

hcw <- hcw %>% filter(hospit == 0)
hcw <- hcw %>% mutate(donation_date = as.Date(date_sample_collection, format = "%m/%d/%y")) %>% rename(result = sc, donorid = id)

hcw <- hcw %>% group_by(donorid) %>% arrange(donation_date) %>%
  mutate(next_result = lead(result), next_date = lead(donation_date), halflife = -log(2)*as.numeric(next_date - donation_date)/(log(next_result) - log(result)))
resmax_plasma <- hcw %>% filter(result >= thr) %>% group_by(donorid) %>% summarise(result_max = max(result)) %>% .$result_max
HL_plasma <- hcw %>% filter(result >= thr, !is.na(halflife), halflife > 0) %>% .$halflife

# -- Plots -- #

mn <- readRDS("data/Bloodbank.rds") %>% filter(blood_center == "HEMOAM", month == 5)
  
dfplot <- rbind(data.frame(var = "Peak S/C", val = resmax_repeat, type = "Repeat Donors"),
                data.frame(var = "Peak S/C", val = resmax_plasma, type = "Plasma Donors"),
                data.frame(var = "Peak S/C", val = mn %>% filter(result >= 0.49) %>% .$result, type = "Manaus in May 2020"),
                data.frame(var = "Half-life (days)", val = HL_repeat, type = "Repeat Donors"),
                data.frame(var = "Half-life (days)", val = HL_plasma, type = "Plasma Donors"))

g1 <- ggplot(dfplot %>% filter(var == "Half-life (days)"), aes(x = val, fill = type)) + 
  geom_histogram(alpha = 1, aes(y = ..density..), binwidth = 5) + 
  geom_density(alpha = 0.3, n = 2^15) + 
  facet_grid(type ~ var, scales = "free") +
  geom_boxplot(aes(x = val, y = 0), width = 0.001, outlier.shape = NA) + 
  coord_cartesian(xlim = c(0, 365)) + theme_bw() + scale_fill_manual(values = mypalette) + 
  theme(legend.position = "none", strip.background = element_rect(fill="#f1e1ab")) + labs(x = "", y = "")

g2 <- ggplot(dfplot %>% filter(var == "Peak S/C"), aes(x = val, fill = type)) + 
  geom_histogram(alpha = 1, aes(y = ..density..), binwidth = 0.5) + 
  geom_density(alpha = 0.3, n = 2^15) + 
  facet_grid(type ~ var, scales = "free") + geom_boxplot(aes(x = val, y = 0), width = 0.01, outlier.shape = NA) + 
  theme_bw() + scale_fill_manual(values = c(mypalette[3], mypalette[1:2])) + 
  theme(legend.position = "none", strip.background = element_rect(fill="#f1e1ab")) + labs(x = "", y = "")

g11 <- ggplot(dfplot %>% filter(var == "Half-life (days)", type == "Repeat Donors"), aes(x = val, fill = type)) + 
  geom_histogram(alpha = 1, aes(y = ..density..), binwidth = 5) + 
  geom_density(alpha = 0.3, n = 2^15) + 
  geom_boxplot(aes(x = val, y = 0), width = 0.001, outlier.shape = NA) + 
  coord_cartesian(xlim = c(0, 365)) + theme_bw() + scale_fill_manual(values = mypalette[1]) + 
  theme(legend.position = "none", strip.background = element_rect(fill="#f1e1ab")) + labs(x = "", y = "") +
  coord_cartesian(xlim = c(0, 365), ylim = c(0, 0.02))

g12 <- ggplot(dfplot %>% filter(var == "Half-life (days)", type == "Plasma Donors"), aes(x = val, fill = type)) + 
  geom_histogram(alpha = 1, aes(y = ..density..), binwidth = 5) + 
  geom_density(alpha = 0.3, n = 2^15) + 
  geom_boxplot(aes(x = val, y = 0), width = 0.001, outlier.shape = NA) + 
  coord_cartesian(xlim = c(0, 365)) + theme_bw() + scale_fill_manual(values = mypalette[2]) + 
  theme(legend.position = "none", strip.background = element_rect(fill="#f1e1ab")) + labs(x = "", y = "") +
  coord_cartesian(xlim = c(0, 365), ylim = c(0, 0.02))


g21 <- ggplot(dfplot %>% filter(var == "Peak S/C", type == "Repeat Donors"), aes(x = val, fill = type)) + 
  geom_histogram(alpha = 1, aes(y = ..density..), binwidth = 0.5) + 
  geom_density(alpha = 0.3, n = 2^15) + geom_boxplot(aes(x = val, y = 0), width = 0.01, outlier.shape = NA) + 
  theme_bw() + scale_fill_manual(values = mypalette[1]) + 
  theme(legend.position = "none", strip.background = element_rect(fill="#f1e1ab")) + labs(x = "", y = "") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 0.2))

g22 <- ggplot(dfplot %>% filter(var == "Peak S/C", type == "Plasma Donors"), aes(x = val, fill = type)) + 
  geom_histogram(alpha = 1, aes(y = ..density..), binwidth = 0.5) + 
  geom_density(alpha = 0.3, n = 2^15)  + geom_boxplot(aes(x = val, y = 0), width = 0.01, outlier.shape = NA) + 
  theme_bw() + scale_fill_manual(values = mypalette[2]) + 
  theme(legend.position = "none", strip.background = element_rect(fill="#f1e1ab")) + labs(x = "", y = "") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 0.2))


g23 <- ggplot(dfplot %>% filter(var == "Peak S/C", type == "Manaus in May 2020"), aes(x = val, fill = type)) + 
  geom_histogram(alpha = 1, aes(y = ..density..), binwidth = 0.5) + 
  geom_density(alpha = 0.3, n = 2^15) + geom_boxplot(aes(x = val, y = 0), width = 0.01, outlier.shape = NA) + 
  theme_bw() + scale_fill_manual(values = mypalette[3]) + 
  theme(legend.position = "none", strip.background = element_rect(fill="#f1e1ab")) + labs(x = "", y = "") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 0.2))

plot_grid(g12, g22, g11, g21, NULL, g23, ncol = 2)

# --- Line plots --- #

df <- read.csv("data/recurrent_donors_hemoam_2021-09-23.csv", sep = ";")
df$donation_date <- as.Date(paste0(df$donyr, "-", df$donmo, "-", df$donda))
df <- df %>% filter(covid_result != "", assay != "COV-2IgGII")
df$result <- as.numeric(sapply(str_split(df$covid_result, pattern = " "), function(x) x[1]))

df <- df %>% group_by(donorid) %>% mutate(ndon = n())
dfplot <- rbind(df %>% filter(ndon >= 2) %>% select(donorid, donation_date, result) %>% mutate(type = "Repeat Blood Donors"),
                hcw %>% select(donorid, donation_date, result) %>% mutate(type = "Convalescent Plasma Donors"))

ggplot(dfplot, aes(x = donation_date, y = result, group = donorid)) +
  geom_line(lwd = 0.5, alpha = 0.1) + geom_hline(yintercept = 0.49) +  geom_point(aes(colour = (result >= 0.49)), size = 1, alpha = 0.5) +
  facet_wrap(~ type, nrow = 2, ncol = 1) + labs(x = "Donation Date", y = "Signal-to-cutoff (S/C)") + theme_bw() + scale_y_log10() +
  scale_colour_manual("", values = mypalette[5:6]) + theme(legend.position = "none",  strip.background = element_rect(fill="#f1e1ab"))

ggsave(last_plot(), file = "figs/lineplots.pdf", width = 10, height = 10)
ggsave(last_plot(), file = "figs/lineplots.png", width = 10, height = 10)


# --- Figure 1 --- #

df <- df %>% group_by(donorid) %>%  mutate(resmax = max(result))
dff <- df %>% filter(resmax >= 0.49) %>% group_by(donorid) %>% arrange(donation_date) %>% mutate(date_first = donation_date[which(result >= 0.49)[1]])
dff <- dff %>% mutate(lag_first = as.numeric(donation_date - date_first))

hcw <- hcw %>% group_by(donorid) %>% arrange(donation_date) %>% mutate(date_first = donation_date[which(result >= 0.49)[1]])
hcw <- hcw %>% mutate(lag_first = as.numeric(donation_date - date_first))
dff <- rbind(dff %>% select(lag_first, result, donorid) %>% mutate(type = "Repeat whole blood donors"),
             hcw %>% select(lag_first, result, donorid) %>% mutate(type = "Convalescent plasma Donors")
             )

lagmax <- dff %>% filter(!is.na(lag_first)) %>% .$lag_first %>% max()

g31 <- ggplot(dff %>% filter(lag_first >= 0, type == "Repeat whole blood donors"), aes(x = lag_first, y = result, group = donorid)) + 
  geom_line(alpha = 0.1, colour = mypalette[1]) + theme_bw() + labs(x = "Days from first positive sample", y = "Signal-to-cutoff (S/C)") +
  scale_y_log10() + coord_cartesian(xlim = c(0, lagmax), ylim = c(0.01, 10)) + geom_hline(yintercept = 0.49, linetype = 'dashed')

g32 <- ggplot(dff %>% filter(lag_first >= 0, type == "Convalescent plasma Donors"), aes(x = lag_first, y = result, group = donorid)) + 
  geom_line(alpha = 0.2, colour = mypalette[2]) + theme_bw() + 
  scale_colour_manual(values = mypalette[2]) + labs(x = "Days from first positive sample", y = "Signal-to-cutoff (S/C)") +
  scale_y_log10() + geom_hline(yintercept = 0.49, linetype = 'dashed') + coord_cartesian(xlim = c(0, lagmax), ylim = c(0.01, 10))

plot_grid(g32, g11, g31, g21, NULL, g23, ncol = 2)
plot_grid(g32, g12, g22, g31, g11, g21, NULL, NULL, g23, ncol = 3)

ggsave(last_plot(), file = "figs/repeat_comparison_may.pdf", width = 15, height = 10)


dfplot <- read.csv("pneg_week.csv")
dfplot <- dfplot %>% group_by(type) %>% mutate(pneg = pneg/sum(pneg))
Ns <- 1E5
dfplot <- dfplot %>% group_by(type) %>% 
  summarise(S = sample(date, Ns, replace = TRUE, prob = pneg))


dfplot <- dfplot %>% filter(type != "Repeat Donors (deaths)")
dfplot$type <- factor(dfplot$type, levels = c("Repeat Donors", "Plasma Donors"))

ggplot(dfplot, aes(y = type, x = S, fill = type)) + 
  geom_boxplot(outlier.shape = NA) +
  theme_bw() + scale_fill_manual(values = mypalette[1:2]) + labs(x = "Interval between seroconversion and seroreversion", y = "") +
  theme(legend.position = "none")
ggsave(last_plot(), file = "figs/tts_repeat_plasma.pdf", width = 5, height = 2)
