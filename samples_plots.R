library(tidyverse)
library(RColorBrewer)
source("Mamoeiro.R")
cities <- c("HEMOAM", "FPS", "HEMOCE", "HEMOBA", "HEMOMINAS", "HEMOPE", "HEMEPAR", "HEMORIO")

fields <- paste0("drho[42, ", 1:10, "]")
fieldnames <- c("F_15-24", "F_25-34", "F_35-44", "F_45-54", "F_55-69",
                "M_15-24", "M_25-34", "M_35-44", "M_45-54", "M_55-69")

dfplot <- list()
dfplot.sex <- list()
dfplot.age <- list()
i <- 1

pop <- readRDS("data/population2020_eightcities.rds") %>% mutate(code = floor(code/10))
codes <- c("HEMOAM" = 130260, "FPS" = 355030, "HEMOCE" = 230440, "HEMOBA" = 292740, 
             "HEMOMINAS" = 310620, "HEMOPE" = 261160, "HEMEPAR" = 410690, "HEMORIO" = 330455)

for(city in cities) {
  dff_inc <- readRDS(paste0("outputs/samples_", city, "_repeat_week_agesex.rds"))
  dff <- 0*dff_inc
  pop.city <- pop %>% filter(code == codes[city]) %>% arrange(age_sex) %>% .$population
    
  for(i in 1:10) {
    dff[, paste0("drho[", 1, ", ", i, "]")] <- dff_inc[, paste0("drho[", 1, ", ", i, "]")]
    for(j in 2:42) {
      dff[, paste0("drho[", j, ", ", i, "]")] <- dff[, paste0("drho[", j - 1, ", ", i, "]")] + dff_inc[, paste0("drho[", j, ", ", i, "]")]
    }
  }
  
  prev.F <- pop.city[1]*dff[, "drho[42, 1]"] + pop.city[2]*dff[, "drho[42, 2]"] + pop.city[3]*dff[, "drho[42, 3]"] + 
    pop.city[4]*dff[, "drho[42, 4]"] + pop.city[5]*dff[, "drho[42, 5]"] 
  prev.M <- pop.city[6]*dff[, "drho[42, 6]"] + pop.city[7]*dff[, "drho[42, 7]"] + pop.city[8]*dff[, "drho[42, 8]"] + 
    pop.city[9]*dff[, "drho[42, 9]"] + pop.city[10]*dff[, "drho[42, 10]"]
  prev.F <- prev.F/sum(pop.city[1:5])
  prev.M <- prev.M/sum(pop.city[6:10])
  
  
  prev.young <- pop.city[1]*dff[, "drho[42, 1]"] + pop.city[2]*dff[, "drho[42, 2]"] + pop.city[3]*dff[, "drho[42, 3]"] + 
    pop.city[4]*dff[, "drho[42, 4]"] + pop.city[6]*dff[, "drho[42, 6]"] + pop.city[7]*dff[, "drho[42, 7]"] + 
    pop.city[8]*dff[, "drho[42, 8]"] + pop.city[9]*dff[, "drho[42, 9]"]
  
  prev.old <- pop.city[5]*dff[, "drho[42, 5]"] + pop.city[10]*dff[, "drho[42, 10]"]
  
  prev.young <- prev.young/sum(pop.city[c(1,2,3,4,6,7,8,9)])
  prev.old <- prev.old/sum(pop.city[c(5,10)])
  
  
  dff <- dff[, fields]/dff[, "drho[42, 6]"]
  dfplot[[city]] <- rbind(data.frame(value = sapply(1:ncol(dff), function(x) quantile(dff[, x], 0.025)), quantile = 0.025, age_sex = fieldnames),
                          data.frame(value = sapply(1:ncol(dff), function(x) quantile(dff[, x], 0.25)), quantile = 0.25, age_sex = fieldnames),
                          data.frame(value = sapply(1:ncol(dff), function(x) quantile(dff[, x], 0.5)), quantile = 0.5, age_sex = fieldnames),
                          data.frame(value = sapply(1:ncol(dff), function(x) quantile(dff[, x], 0.75)), quantile = 0.75, age_sex = fieldnames),
                          data.frame(value = sapply(1:ncol(dff), function(x) quantile(dff[, x], 0.975)), quantile = 0.975, age_sex = fieldnames)) %>%
    mutate(city = city)
  
  dfplot.sex[[city]] <- data.frame(
    q1.M = quantile(prev.M, 0.025),
    q2.M = quantile(prev.M, 0.25),
    q3.M = quantile(prev.M, 0.5),
    q4.M = quantile(prev.M, 0.75),
    q5.M = quantile(prev.M, 0.975),
    q1.F = quantile(prev.F, 0.025),
    q2.F = quantile(prev.F, 0.25),
    q3.F = quantile(prev.F, 0.5),
    q4.F = quantile(prev.F, 0.75),
    q5.F = quantile(prev.F, 0.975),
    q1 = quantile(prev.M/prev.F, 0.025),
    q2 = quantile(prev.M/prev.F, 0.25),
    q3 = quantile(prev.M/prev.F, 0.5),
    q4 = quantile(prev.M/prev.F, 0.75),
    q5 = quantile(prev.M/prev.F, 0.975),
    city = city)
  
  dfplot.age[[city]] <- data.frame(
    q1.young = quantile(prev.young, 0.025),
    q2.young = quantile(prev.young, 0.25),
    q3.young = quantile(prev.young, 0.5),
    q4.young = quantile(prev.young, 0.75),
    q5.young = quantile(prev.young, 0.975),
    q1.old = quantile(prev.old, 0.025),
    q2.old = quantile(prev.old, 0.25),
    q3.old = quantile(prev.old, 0.5),
    q4.old = quantile(prev.old, 0.75),
    q5.old = quantile(prev.old, 0.975),
    q1 = quantile(prev.young/prev.old, 0.025),
    q2 = quantile(prev.young/prev.old, 0.25),
    q3 = quantile(prev.young/prev.old, 0.5),
    q4 = quantile(prev.young/prev.old, 0.75),
    q5 = quantile(prev.young/prev.old, 0.975),
    city = city)
}


dfplot <- do.call(rbind, dfplot)
dfplot.sex <- do.call(rbind, dfplot.sex)
dfplot.age <- do.call(rbind, dfplot.age)

write.csv(dfplot.sex, file = "prev_december_sex.csv")
write.csv(dfplot.age, file = "prev_december_age.csv")

saveRDS(dfplot, "RR_prev_december.rds")


dfplot <- dfplot %>% group_by(age_sex, city) %>% summarise(q1 = value[which(quantile == 0.025)[1]],
                                                           q2 = value[which(quantile == 0.25)[1]],
                                                           q3 = value[which(quantile == 0.5)[1]],
                                                           q4 = value[which(quantile == 0.75)[1]],
                                                           q5 = value[which(quantile == 0.975)[1]])
dfplot <- dfplot %>% mutate(cityname = citynames[city])

palM <- colorRampPalette(brewer.pal(9, "Blues"))
palF <- colorRampPalette(brewer.pal(9, "Reds"))
palM <- palM(9)[c(3, 5, 6, 8, 9)]
palF <- palF(9)[c(3, 5, 6, 8, 9)]

ggplot(dfplot, aes(x = cityname, y = q3, fill = age_sex)) +  geom_errorbar(aes(ymin = q1, ymax = q5), position = "dodge2", width = 0.8, lwd = 0.1) + 
  geom_crossbar(aes(ymin = q2, ymax = q4), position = "dodge2", width = 0.8, lwd = 0.1) +
  theme_bw() + scale_fill_manual("", values = c(palF, palM)) + labs(x = "", y = "Relative Risk") + geom_hline(yintercept = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave(last_plot(), file = "figs/Prevalence_RR2.png", width = 11, height = 3, dpi = 600)
ggsave(last_plot(), file = "figs/Prevalence_RR2.pdf", width = 11, height = 3, dpi = 600)


# --- By sex only --- #
