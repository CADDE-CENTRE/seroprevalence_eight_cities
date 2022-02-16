library(tidyverse)
library(lubridate)


compute.decay.rate <- function(df) #Estimates the decay rate for each donation from a repeat blood donor
{
  df$row_id <- 1:nrow(df)
  df$next_id <- NA
  for(i in 1:nrow(df))
  {
    dfs <- df %>% filter(donorid == df$donorid[i], donation_date > df$donation_date[i]) %>% arrange(donation_date)
    if(nrow(dfs) > 0)
      df$next_id[i] <- dfs$row_id[1]
  }
  
  # --- decay rate --- #
  df$decay_rate <- NA
  for(i in 1:nrow(df))
  {
    if(!is.na(df$next_id[i]))
      df$decay_rate[i] <- (log(df$result[df$next_id[i]]) - log(df$result[i]))/as.numeric(df$donation_date[df$next_id[i]] - df$donation_date[i])
  }
  return(df)
}



calculate.tts.coefs <- function(y, y0, tau, thr) #Calculates time to seroreversion
{
  a = (y-y0)/tau
  return((thr - y0)/a)
}


calculate.seroreversion.date <- function(df, thr = 0.49, median.halflife = FALSE)
  #Calculates the date of seroreversion
  #If median.halflife = TRUE, the median halflife for each donor will be used 
  #in the exponential interpolation. In the paper, we use median.halflife = FALSE.
{
  df <- df %>% group_by(donorid) %>% mutate(result_max = max(result))
  df <- df %>% filter(result_max >= thr)
  
  if(median.halflife)
  {
    dff <- df %>% mutate(halflife = -log(2)/decay_rate) %>% filter(!is.na(halflife), halflife > 0, result >= thr)
    HLmed <- median(dff$halflife)
    print(HLmed)
  }
  ids <- unique(df$donorid)
  HL <- rep(NA, length(ids))
  serorev_date <- rep(NA, length(ids))
  serorev_date_lin <- rep(NA, length(ids))
  halflife <- rep(NA, length(ids))
  serorev_obs <- rep(NA, length(ids))
  for(i in 1:length(ids))
  {
    print(paste0("Processing patient ", i, "/", length(ids)))
    pat <- df %>% filter(donorid == ids[i]) %>% arrange(donation_date)
    
    kfirst <- which(pat$result >= thr)[1]
    if(kfirst < nrow(pat)) {
      klast <- tail(which(pat$result >= thr), 1)
      kmax <- which(pat$result == max(pat$result))[1]
      if(klast < nrow(pat))  #Seroreversion was observed!
      {
        ka <- klast
        kb <- klast + 1
        serorev_obs[i] <- TRUE
      }
      else {
        ka <- klast - 1      #Seroreversion didn't happen yet
        kb <- klast
        serorev_obs[i] <- FALSE
      }
      
      Dt <- as.numeric(pat$donation_date[kb] - pat$donation_date[ka])
      Dy.exp <- (log(pat$result[kb]) - log(pat$result[ka]))
      Dy.lin <- (pat$result[kb] - pat$result[ka])
      
      if(Dy.lin < 0) #In this case Dy.exp is also < 0 
      {
        if(median.halflife & (serorev_obs[i] == FALSE))
        {
          serorev_date[i] <- as.character(pat$donation_date[klast] + HLmed*(log(pat$result[klast]) - log(thr)))
        } else {
          serorev_date[i] <- as.character(pat$donation_date[ka] + (log(thr) - log(pat$result[ka]))*Dt/Dy.exp)
        }
        serorev_date_lin[i] <- as.character(pat$donation_date[ka] + (thr - pat$result[ka])*Dt/Dy.lin)
        halflife[i] <- -log(2)*Dt/Dy.exp
      }
    }
  
  }
  df <- left_join(df, data.frame(donorid = ids, serorev_date = as.Date(serorev_date), 
                                 serorev_date_lin = as.Date(serorev_date_lin), 
                                 halflife = halflife, serorev_obs = serorev_obs), 
                  by = "donorid")
  return(df)
}
