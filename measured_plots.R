library(tidyverse)
source("Mamoeiro.R")

codes <- c("HEMOAM" = 130260, "FPS" = 355030, "HEMOCE" = 230440, "HEMOBA" = 292740, 
           "HEMOMINAS" = 310620, "HEMOPE" = 261160, "HEMEPAR" = 410690, "HEMORIO" = 330455)


df <- readRDS("data/Bloodbank.rds")

pop <- readRDS("data/population2020.rds")
pop <- pop %>% mutate(code = floor(code/10)) %>% 
  filter(code %in% codes, age_bin %in% c("15-20", "20-25", "25-30", "30-35", "35-40", "40-45", "45-50", "50-55", "55-60", "60-65"))
agemap <- c("15-20" = "15-24", "20-25" = "15-24", "25-30" = "25-34", "30-35" = "25-34", 
            "35-40" = "35-44", "40-45" = "35-44", "45-50" = "45-54", "50-55" = "45-54", "55-60" = "55-69", "60-65" = "55-69", "65-70" = "55-69")
pop <- pop %>% mutate(new_bin = agemap[age_bin])
pop <- pop %>% group_by(code, new_bin, CS_SEXO) %>% summarise(population = sum(population)) %>% 
  ungroup() %>% rename(age_bin = new_bin) %>% mutate(age_sex = paste0(CS_SEXO, "_", age_bin)) %>% arrange(age_sex)



df <- df %>% filter(!is.na(idade), idade != "", !is.na(sex), sex != "", idade >= 16, idade <= 70, !is.na(result), month <= 15)
df$age_group <- cut(df$idade, c(0, 25, 35, 45, 55, 100), 
                             labels = c("15-24", "25-34", "35-44", "45-54", "55-69"))
df$age_sex <- paste0(df$sex, "_", df$age_group)

df <- df %>% mutate(code = codes[blood_center]) %>% left_join(pop, by = c("code", "age_sex"))
df <- df %>% left_join(pop %>% group_by(code) %>% summarise(poptotal = sum(population)), by = "code")

df <- df %>% group_by(blood_center, month, age_sex) %>% summarise(population = population[1], 
                                                                  poptotal = poptotal[1],
                                                                  npos = sum(result >= 0.49),
                                                                  ntests = n())
df <- df %>% mutate(age_group = substring(age_sex, 3, nchar(age_sex)),
                    sex = substring(age_sex, 1, 1))




dfplot <- df %>% group_by(blood_center, month) %>% summarise(npos.all = sum(npos), ntests.all = sum(ntests), 
                                                             q1.crude = qbeta(0.025, 1+npos.all, 1 + ntests.all - npos.all),
                                                             q2.crude = qbeta(0.25, 1+npos.all, 1 + ntests.all - npos.all),
                                                             q3.crude = qbeta(0.5, 1+npos.all, 1 + ntests.all - npos.all),
                                                             q4.crude = qbeta(0.75, 1+npos.all, 1 + ntests.all - npos.all),
                                                             q5.crude = qbeta(0.975, 1+npos.all, 1 + ntests.all - npos.all),
                                                             q1.corr = NA,
                                                             q2.corr = NA,
                                                             q3.corr = NA,
                                                             q4.corr = NA,
                                                             q5.corr = NA)


Nsamples <- 1E5
se <- rbeta(Nsamples, 1 + 189, 1 + 19)
sp <- rbeta(Nsamples, 1 + 801, 1 + 20)
for(i in 1:nrow(dfplot)) {
  
  I <- which((df$blood_center == dfplot$blood_center[i]) & (df$month == dfplot$month[i]))
  rho <- rep(0, Nsamples)
  for(k in I){
    rho <- rho + df$population[k]*rbeta(Nsamples, 1 + df$npos[k], 1 + df$ntests[k] - df$npos[k])/df$poptotal[k]
  }
  rho <- (rho + sp - 1)/(se + sp - 1)
  dfplot$q1.corr[i] <- quantile(rho, 0.025)
  dfplot$q2.corr[i] <- quantile(rho, 0.25)
  dfplot$q3.corr[i] <- quantile(rho, 0.5)
  dfplot$q4.corr[i] <- quantile(rho, 0.75)
  dfplot$q5.corr[i] <- quantile(rho, 0.975)
}


dfplot <- rbind(dfplot %>% rename(q1 = q1.crude,
                                  q2 = q2.crude,
                                  q3 = q3.crude,
                                  q4 = q4.crude,
                                  q5 = q5.crude) %>% 
                  select(blood_center, month, q1, q2, q3, q4, q5) %>% mutate(type = "Crude"),
                dfplot %>% rename(q1 = q1.corr,
                                  q2 = q2.corr,
                                  q3 = q3.corr,
                                  q4 = q4.corr,
                                  q5 = q5.corr) %>% 
                  select(blood_center, month, q1, q2, q3, q4, q5) %>% mutate(type = "Adjusted"))

dfplot <- dfplot %>% mutate(city = citynames[blood_center])
ggplot(dfplot, aes(x = month, y = 100*q3, fill = type)) + geom_errorbar(aes(ymin = 100*q1, ymax = 100*q5), position = "dodge2", lwd = 0.1) + 
  geom_crossbar(aes(ymin = 100*q2, ymax = 100*q4), position = "dodge2", lwd = 0.1) +
  facet_wrap(~ city, nrow = 4) + scale_fill_manual("", values = mypalette) + theme_bw() + theme(strip.background =element_rect(fill="#f1e1ab")) +
  scale_x_continuous(breaks = 3:15, 
                     labels = c("Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Fev", "Mar")) +
  labs(x = "Month", y = "Prevalence (%)")

ggsave(last_plot(), file = "figs/measured_prevalence_month.png", width = 10, height = 10)
ggsave(last_plot(), file = "figs/measured_prevalence_month.pdf", width = 10, height = 10)


df <- df %>% mutate(city = citynames[blood_center])
library(RColorBrewer)
palM <- colorRampPalette(brewer.pal(5, "Blues"))
palF <- colorRampPalette(brewer.pal(5, "Reds"))

ggplot(df, aes(x = month, y = ntests, fill = age_sex)) + geom_bar(stat = "identity") + 
  facet_wrap(~ city, nrow = 4) + scale_fill_manual("", values = c(palF(5), palM(5))) + theme_bw() + theme(strip.background =element_rect(fill="#f1e1ab")) +
  scale_x_continuous(breaks = 3:15, 
                     labels = c("Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Fev", "Mar")) +
  labs(x = "Month", y = "Number of monthly tests")
ggsave(last_plot(), file = "figs/ntests.png", width = 10, height = 10)
ggsave(last_plot(), file = "figs/ntests.pdf", width = 10, height = 10)

#Table
df.out <- dfplot %>% filter((month %in% c(12, 15)) | ((month == 14) & (city == "Recife"))) %>% 
                              mutate(txt = paste0(100*round(q3, 3), " (", 100*round(q1, 3), " - ", 100*round(q5, 3), ")"))
df.out$month[df.out$month == 14] <- 15
df.out <- df.out %>% arrange(type, month, city) 
write.csv(df.out, "measured_seroprevalence_month.csv")

write.csv(dfplot, "measured_seroprevalence_month_all.csv")

# --- Seroprevalence disaggregated by age/sex --- #

dfplot.agesex.all <- df %>% group_by(blood_center, age_sex) %>% summarise(npos.all = sum(npos), ntests.all = sum(ntests), 
                                                                          q1.crude = qbeta(0.025, 1+npos.all, 1 + ntests.all - npos.all),
                                                                          q2.crude = qbeta(0.25, 1+npos.all, 1 + ntests.all - npos.all),
                                                                          q3.crude = qbeta(0.5, 1+npos.all, 1 + ntests.all - npos.all),
                                                                          q4.crude = qbeta(0.75, 1+npos.all, 1 + ntests.all - npos.all),
                                                                          q5.crude = qbeta(0.975, 1+npos.all, 1 + ntests.all - npos.all))

dfplot.agesex <- df %>% group_by(blood_center, month, age_sex) %>% summarise(npos.all = sum(npos), ntests.all = sum(ntests), 
                                                                        q1.crude = qbeta(0.025, 1+npos.all, 1 + ntests.all - npos.all),
                                                                        q2.crude = qbeta(0.25, 1+npos.all, 1 + ntests.all - npos.all),
                                                                        q3.crude = qbeta(0.5, 1+npos.all, 1 + ntests.all - npos.all),
                                                                        q4.crude = qbeta(0.75, 1+npos.all, 1 + ntests.all - npos.all),
                                                                        q5.crude = qbeta(0.975, 1+npos.all, 1 + ntests.all - npos.all),
                                                                        q1.corr = NA,
                                                                        q2.corr = NA,
                                                                        q3.corr = NA,
                                                                        q4.corr = NA,
                                                                        q5.corr = NA)

Nsamples <- 1E5
se <- rbeta(Nsamples, 1 + 189, 1 + 19)
sp <- rbeta(Nsamples, 1 + 801, 1 + 20)
for(i in 1:nrow(dfplot.agesex)) {
  k <- which((df$blood_center == dfplot.agesex$blood_center[i]) & (df$month == dfplot.agesex$month[i]) & (df$age_sex == dfplot.agesex$age_sex[i]))
  rho <- rbeta(Nsamples, 1 + df$npos[k], 1 + df$ntests[k] - df$npos[k])
  rho <- (rho + sp - 1)/(se + sp - 1)
  dfplot.agesex$q1.corr[i] <- quantile(rho, 0.025)
  dfplot.agesex$q2.corr[i] <- quantile(rho, 0.25)
  dfplot.agesex$q3.corr[i] <- quantile(rho, 0.5)
  dfplot.agesex$q4.corr[i] <- quantile(rho, 0.75)
  dfplot.agesex$q5.corr[i] <- quantile(rho, 0.975)
}


dfplot.agesex <- rbind(dfplot.agesex %>% rename(q1 = q1.crude,
                                  q2 = q2.crude,
                                  q3 = q3.crude,
                                  q4 = q4.crude,
                                  q5 = q5.crude) %>% 
                  select(blood_center, month, age_sex, q1, q2, q3, q4, q5) %>% mutate(type = "Crude"),
                dfplot.agesex %>% rename(q1 = q1.corr,
                                  q2 = q2.corr,
                                  q3 = q3.corr,
                                  q4 = q4.corr,
                                  q5 = q5.corr) %>% 
                  select(blood_center, month, age_sex, q1, q2, q3, q4, q5) %>% mutate(type = "Adjusted"))

dfplot.agesex <- dfplot.agesex %>% mutate(city = citynames[blood_center])

ggplot(dfplot.agesex %>% filter(type == "Crude"), aes(x = month, y = 100*q3, fill = age_sex)) + 
  geom_errorbar(aes(ymin = 100*q1, ymax = 100*q5), position = "dodge2", lwd = 0.1) + 
  geom_crossbar(aes(ymin = 100*q2, ymax = 100*q4), position = "dodge2", lwd = 0.1) +
  facet_wrap(~ city, nrow = 4) + scale_fill_manual("", values = c(palF(5), palM(5))) + 
  theme_bw() + theme(strip.background =element_rect(fill="#f1e1ab")) +
  scale_x_continuous(breaks = 3:15, 
                     labels = c("Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Fev", "Mar")) +
  labs(x = "Month", y = "Prevalence (%)")
write.csv(dfplot.agesex, file = "crudeseroprevalence_agesex.csv")

# --- Compare crude and adjusted seroprevalence by age and sex --- #

adj <- list()
for(center in names(citynames)) {
  adj[[center]] <- read.csv(paste0("outputs/rho_est_", center, "_repeat_week_agesex.csv")) %>% mutate(center = center)
}
adj <- do.call(rbind, adj)

adj <- adj %>% group_by(center, month, age_sex) %>% summarise(q1 = value[quantile == 0.025],
                                                                    q2 = value[quantile == 0.25], 
                                                                    q3 = value[quantile == 0.5], 
                                                                    q4 = value[quantile == 0.75], 
                                                                    q5 = value[quantile == 0.975])
adj <- adj %>% rename(blood_center = center) %>% mutate(city = citynames[blood_center], type = "Seroreversion",
                                                        dt = as.Date("2020-01-15") + 7*(month - 1))

library(lubridate)
dfplot.agesex <- dfplot.agesex %>% mutate(dt = as.Date("2020-01-15") %m+% months((month - 1)))
dfplot.agesex.serorev <- rbind(dfplot.agesex, adj)

dfplot.agesex.serorev <- dfplot.agesex.serorev %>% mutate(sex = substring(age_sex, 1, 1), age = substring(age_sex, 3, nchar(age_sex)),
                                                          sex = ifelse(sex == "M", "Male", "Female"))

dfplot.agesex.serorev <- dfplot.agesex.serorev %>% mutate(facet_group = paste0(city, ": ", sex, " ", age))
dfplot.agesex.serorev$type[dfplot.agesex.serorev$type == "Adjusted"] <- "Adjusted (no correction for seroreversion)"
dfplot.agesex.serorev$type[dfplot.agesex.serorev$type == "Seroreversion"] <- "Adjusted (correted for seroreversion)"

ggplot(dfplot.agesex.serorev %>% filter(type != "Adjusted (correted for seroreversion)"), aes(x = dt, y = 100*q3, fill = type)) + 
  geom_errorbar(aes(ymin = 100*q1, ymax = 100*q5), position = "dodge2", lwd = 0.01) + 
  geom_crossbar(aes(ymin = 100*q2, ymax = 100*q4), position = "dodge2", lwd = 0.01) + 
  geom_line(data = dfplot.agesex.serorev %>% filter(type == "Adjusted (correted for seroreversion)"), aes(x = dt, y = 100*q3)) +
  geom_ribbon(data = dfplot.agesex.serorev %>% filter(type == "Adjusted (correted for seroreversion)"), aes(x = dt, ymin = 100*q1, ymax = 100*q5, fill = type), alpha = 0.5) +
  facet_wrap(~ facet_group, nrow = 16) + scale_fill_manual("", values = mypalette) + 
  theme_bw() + 
  scale_x_continuous(breaks = unique(dfplot.agesex$dt), 
                     labels = c("Mar/20", "Apr/20", "May/20", "Jun/20", "Jul/20", "Aug/20", "Sep/20", "Oct/20", "Nov/20", "Dec/20", "Jan/21", "Fev/21", "Mar/21")) +
  labs(x = "Month", y = "Seroprevalence (%)") + theme(strip.background =element_rect(fill="#f1e1ab"), legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1, size = 5))

ggsave(last_plot(), file = "figs/seroprevalence_all_agesex_mes.png", width = 10, height = 20)
ggsave(last_plot(), file = "figs/seroprevalence_all_agesex_mes.pdf", width = 10, height = 20)

# --- Seroprevalence disaggregated by age --- #


dfplot.age.all <- df %>% group_by(blood_center, age_group) %>% summarise(npos.all = sum(npos), ntests.all = sum(ntests), 
                                                                          q1.crude = qbeta(0.025, 1+npos.all, 1 + ntests.all - npos.all),
                                                                          q2.crude = qbeta(0.25, 1+npos.all, 1 + ntests.all - npos.all),
                                                                          q3.crude = qbeta(0.5, 1+npos.all, 1 + ntests.all - npos.all),
                                                                          q4.crude = qbeta(0.75, 1+npos.all, 1 + ntests.all - npos.all),
                                                                          q5.crude = qbeta(0.975, 1+npos.all, 1 + ntests.all - npos.all))

dfplot.age <- df %>% group_by(blood_center, month, age_group) %>% summarise(npos.all = sum(npos), ntests.all = sum(ntests), 
                                                                  q1.crude = qbeta(0.025, 1+npos.all, 1 + ntests.all - npos.all),
                                                                  q2.crude = qbeta(0.25, 1+npos.all, 1 + ntests.all - npos.all),
                                                                  q3.crude = qbeta(0.5, 1+npos.all, 1 + ntests.all - npos.all),
                                                                  q4.crude = qbeta(0.75, 1+npos.all, 1 + ntests.all - npos.all),
                                                                  q5.crude = qbeta(0.975, 1+npos.all, 1 + ntests.all - npos.all),
                                                                  q1.corr = NA,
                                                                  q2.corr = NA,
                                                                  q3.corr = NA,
                                                                  q4.corr = NA,
                                                                  q5.corr = NA)

Nsamples <- 1E5
se <- rbeta(Nsamples, 1 + 189, 1 + 19)
sp <- rbeta(Nsamples, 1 + 801, 1 + 20)
for(i in 1:nrow(dfplot.age)) {
  k <- which((df$blood_center == dfplot.age$blood_center[i]) & (df$month == dfplot.age$month[i]) & (df$age_group == dfplot.age$age_group[i]))
  rho <- rbeta(Nsamples, 1 + df$npos[k], 1 + df$ntests[k] - df$npos[k])
  rho <- (rho + sp - 1)/(se + sp - 1)
  dfplot.age$q1.corr[i] <- quantile(rho, 0.025)
  dfplot.age$q2.corr[i] <- quantile(rho, 0.25)
  dfplot.age$q3.corr[i] <- quantile(rho, 0.5)
  dfplot.age$q4.corr[i] <- quantile(rho, 0.75)
  dfplot.age$q5.corr[i] <- quantile(rho, 0.975)
}


dfplot.age <- rbind(dfplot.age %>% rename(q1 = q1.crude,
                                  q2 = q2.crude,
                                  q3 = q3.crude,
                                  q4 = q4.crude,
                                  q5 = q5.crude) %>% 
                  select(blood_center, month, age_group, q1, q2, q3, q4, q5) %>% mutate(type = "Crude"),
                dfplot.age %>% rename(q1 = q1.corr,
                                  q2 = q2.corr,
                                  q3 = q3.corr,
                                  q4 = q4.corr,
                                  q5 = q5.corr) %>% 
                  select(blood_center, month, age_group, q1, q2, q3, q4, q5) %>% mutate(type = "Adjusted"))

dfplot.age <- dfplot.age %>% mutate(city = citynames[blood_center])
dfplot.age.all <- dfplot.age.all %>% mutate(city = citynames[blood_center])

ggplot(dfplot.age %>% filter(type == "Crude"), aes(x = month, y = 100*q3, fill = age_group)) + geom_line(aes(colour = age_group)) + 
  geom_ribbon(aes(ymin = 100*q1, ymax = 100*q5), alpha = 0.1) + 
  scale_fill_manual("", values = mypalette) + 
  theme_bw() + theme(strip.background =element_rect(fill="#f1e1ab")) + 
  geom_errorbar(data = dfplot.age.all, aes(x = 16, ymin = 100*q1.crude, ymax = 100*q5.crude, group = age_group), inherit.aes = FALSE, position = "dodge2", lwd = 0.1) + 
  geom_crossbar(data = dfplot.age.all, aes(x = 16, ymin = 100*q2.crude, ymax = 100*q4.crude, y = 100*q3.crude, fill = age_group), inherit.aes = FALSE, position = "dodge2", lwd = 0.1) +
  geom_vline(xintercept = 15, linetype="dotted") + 
  facet_wrap(~ city, nrow = 4) + scale_x_continuous(breaks = 3:16, 
                                                    labels = c("Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Fev", "Mar", "All months")) +
  labs(x = "Month", y = "Seroprevalence (%)")
ggsave(last_plot(), file = "figs/crude_seroprevalence_age.pdf", width = 15, height = 10)

# --- Seroprevalence disaggregated by sex --- #



dfplot.sex.all <- df %>% group_by(blood_center, sex) %>% summarise(npos.all = sum(npos), ntests.all = sum(ntests), 
                                                                         q1.crude = qbeta(0.025, 1+npos.all, 1 + ntests.all - npos.all),
                                                                         q2.crude = qbeta(0.25, 1+npos.all, 1 + ntests.all - npos.all),
                                                                         q3.crude = qbeta(0.5, 1+npos.all, 1 + ntests.all - npos.all),
                                                                         q4.crude = qbeta(0.75, 1+npos.all, 1 + ntests.all - npos.all),
                                                                         q5.crude = qbeta(0.975, 1+npos.all, 1 + ntests.all - npos.all))

dfplot.sex <- df %>% group_by(blood_center, month, sex) %>% summarise(npos.all = sum(npos), ntests.all = sum(ntests), 
                                                                         q1.crude = qbeta(0.025, 1+npos.all, 1 + ntests.all - npos.all),
                                                                         q2.crude = qbeta(0.25, 1+npos.all, 1 + ntests.all - npos.all),
                                                                         q3.crude = qbeta(0.5, 1+npos.all, 1 + ntests.all - npos.all),
                                                                         q4.crude = qbeta(0.75, 1+npos.all, 1 + ntests.all - npos.all),
                                                                         q5.crude = qbeta(0.975, 1+npos.all, 1 + ntests.all - npos.all),
                                                                         q1.corr = NA,
                                                                         q2.corr = NA,
                                                                         q3.corr = NA,
                                                                         q4.corr = NA,
                                                                         q5.corr = NA)

Nsamples <- 1E5
se <- rbeta(Nsamples, 1 + 189, 1 + 19)
sp <- rbeta(Nsamples, 1 + 801, 1 + 20)
for(i in 1:nrow(dfplot.sex)) {
  k <- which((df$blood_center == dfplot.sex$blood_center[i]) & (df$month == dfplot.sex$month[i]) & (df$sex == dfplot.sex$sex[i]))
  rho <- rbeta(Nsamples, 1 + df$npos[k], 1 + df$ntests[k] - df$npos[k])
  rho <- (rho + sp - 1)/(se + sp - 1)
  dfplot.sex$q1.corr[i] <- quantile(rho, 0.025)
  dfplot.sex$q2.corr[i] <- quantile(rho, 0.25)
  dfplot.sex$q3.corr[i] <- quantile(rho, 0.5)
  dfplot.sex$q4.corr[i] <- quantile(rho, 0.75)
  dfplot.sex$q5.corr[i] <- quantile(rho, 0.975)
}


dfplot.sex <- rbind(dfplot.sex %>% rename(q1 = q1.crude,
                                                q2 = q2.crude,
                                                q3 = q3.crude,
                                                q4 = q4.crude,
                                                q5 = q5.crude) %>% 
                         select(blood_center, month, sex, q1, q2, q3, q4, q5) %>% mutate(type = "Crude"),
                       dfplot.sex %>% rename(q1 = q1.corr,
                                                q2 = q2.corr,
                                                q3 = q3.corr,
                                                q4 = q4.corr,
                                                q5 = q5.corr) %>% 
                         select(blood_center, month, sex, q1, q2, q3, q4, q5) %>% mutate(type = "Adjusted"))

dfplot.sex <- dfplot.sex %>% mutate(city = citynames[blood_center])
dfplot.sex.all <- dfplot.sex.all %>% mutate(city = citynames[blood_center])

ggplot(dfplot.sex %>% filter(type == "Crude"), aes(x = month, y = 100*q3, fill = sex)) + geom_line(aes(colour = sex)) + 
  geom_ribbon(aes(ymin = 100*q1, ymax = 100*q5), alpha = 0.5) + 
   scale_fill_manual("", values = c(palF(5)[3], palM(5)[3])) + 
  theme_bw() + theme(strip.background =element_rect(fill="#f1e1ab")) + 
  geom_errorbar(data = dfplot.sex.all, aes(x = 16, ymin = 100*q1.crude, ymax = 100*q5.crude, group = sex), inherit.aes = FALSE, position = "dodge2", lwd = 0.1) + 
  geom_crossbar(data = dfplot.sex.all, aes(x = 16, ymin = 100*q2.crude, ymax = 100*q4.crude, y = 100*q3.crude, fill = sex), inherit.aes = FALSE, position = "dodge2", lwd = 0.1) +
  geom_vline(xintercept = 15, linetype="dotted") + 
  facet_wrap(~ city, nrow = 4) + scale_x_continuous(breaks = 3:16, 
                     labels = c("Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Fev", "Mar", "All months")) +
  labs(x = "Month", y = "Seroprevalence (%)")
ggsave(last_plot(), file = "figs/crude_seroprevalence_sex.pdf", width = 15, height = 10)

