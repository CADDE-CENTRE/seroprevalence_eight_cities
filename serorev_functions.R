## load the nimble library and set seed
library('nimble')
library('nimbleSMC')
library(tidyverse)
library(lubridate)
library(readxl)
library(survey)

modelCode.BSC <- function()
{
 
  code <- nimbleCode({
    se ~ dbeta(1+npos.plasma, 1+nneg.plasma)
    spb ~ dbeta(1+npos.prepandemic, 1+nneg.prepandemic) #spb = 1 - specificity
    
    for(m in 1:ncompartments)
    {
      drhon[1:nmonths, m] ~ ddirch(alpha[1:nmonths])
      rho_max[m] ~ dunif(rho_max.lb, rho_max.ub) 
      for(n in 1:nmonths)
        drho[n, m] <- drhon[n,m]*rho_max[m] 
    }
    
    for(n in 1:nmonths) {
      for(m in 1:ncompartments) {
        mu[n,m] <- se*sum(A[n, 1:n]*drho[1:n,m]) + spb*(1 - sum(drho[1:n,m]))
      }
    }
    
    for(n in 1:nmonths) {
      for(m in 1:ncompartments) {
        npos[n,m] ~ dbinom(size = ntests[n,m], prob = mu[n,m])
      }
    }
  })
  
  
  return(code)
}


prevalence.blood_donors <- function(prev, ppos, population, cases_max = 1, cases_prior = NA,
                                    nburnin = 10000, niter = 100000, alpha = NA)
{
  npos.plasma <- 189
  nneg.plasma <- 19
  npos.prepandemic <- 20
  nneg.prepandemic <- 801
  
  ntests = prev.to.matrix(prev, "ntests")
  npos = prev.to.matrix(prev, "npos")
  
  nmonths = nrow(npos)
  ncompartments = ncol(npos)
  
  if(is.na(cases_prior)){
    rho_max.lb <- 0
    rho_max.ub <- cases_max
  }
  else{
    rho_max.lb <- 0
    rho_max.ub <- min(cases_prior, cases_max)
  }
  if(is.na(alpha))
  {
    alpha <- rep(1, nmonths)
  }
  
  A <- ppos.to.matrix(ppos, nmonths)
  
  modelCode <- modelCode.BSC()
  
  constants <- list(
    nmonths = nmonths,
    ntests = as.matrix(ntests),
    ncompartments = ncompartments,
    A = A,
    alpha = alpha,
    rho_max.lb = rho_max.lb,
    rho_max.ub = rho_max.ub,
    npos.plasma = npos.plasma,
    nneg.plasma = nneg.plasma,
    npos.prepandemic = npos.prepandemic,
    nneg.prepandemic = nneg.prepandemic
  )
  
  data <- list(npos = as.matrix(npos))
  inits <- list(rho_max = rep(0.5*(rho_max.lb + rho_max.ub), ncompartments), 
                drhon = matrix(data = 1/nmonths, nrow = nmonths, ncol = ncompartments),
                se = 0.9,
                spb = 0.01
                )
  
  
  nimbleMCMC_samples <- nimbleMCMC(code = modelCode, 
                                   constants = constants, 
                                   data = data, 
                                   inits = inits,
                                   nburnin = nburnin, niter = niter, monitors = c("drho", "se"))

  prev_est <- matrix(data = 0, nrow = nrow(nimbleMCMC_samples), ncol = nmonths)
  drho_est <- matrix(data = 0, nrow = nrow(nimbleMCMC_samples), ncol = nmonths)
  
  for(m in 1:ncompartments) {
    drho_samples = nimbleMCMC_samples[, sapply(1:nmonths, function(x) sprintf(paste0("drho[%d, ", m, "]"), x))]
    rho_samples = apply(drho_samples, 1, cumsum) %>% t()
    prev_est <- prev_est + population[m]*rho_samples/sum(population)
    drho_est <- drho_est + population[m]*drho_samples/sum(population)
  }
  
  drho <- array(0, dim = c(nmonths, ncompartments, nrow(nimbleMCMC_samples)))
  for(i in 1:dim(drho)[1])
    for(j in 1:dim(drho)[2])
        drho[i,j,] <- nimbleMCMC_samples[, sprintf("drho[%d, %d]", i, j)]
  

  colnames(drho) <- prev$age_sex %>% unique() %>% sort()  
  rownames(drho) <- prev$month %>% unique() %>% sort()  
          
  return(list(prev_est = prev_est, cases_est = drho_est, drho = drho, samples = nimbleMCMC_samples))

}


prevalence.blood_donors.hom <- function(prev, ppos, population, nburnin = 10000, niter = 100000, alpha = NA, cases_max = 1)
{
  prev_hom <- prev %>% group_by(month) %>% summarise(ntests = sum(ntests), npos = sum(npos), age_sex = "MF_0-100")  #Homogeneous prevalence
  lst <- prevalence.blood_donors(prev_hom, ppos, sum(population), cases_max, nburnin = nburnin, niter = niter, alpha = alpha)
  lst <- prevalence.blood_donors(prev, ppos, population, cases_max, 2*quantile(lst$prev_est[, ncol(lst$prev_est)], 0.95), 
                                 nburnin = nburnin, niter = niter, alpha = alpha)
  
  
  
  return(lst)
}



mergeAges.blood_donors <- function(pop, age_bin = 5)
{
  age_sex_exclude <- c("M_0-4","F_0-4","M_5-9", "F_5-9", "M_10-14", "F_10-14", "M_70-74", "F_70-74", 
                       "M_75-79", "F_75-79", "M_80-84", "F_80-84", "M_85-89", "F_85-89", "M_90-94", "F_90-94", "M_95-99", "F_95-99",
                       "M_100+", "F_100+")

  if(age_bin == 5)
  {
    age_sex_merge_M <- c("M_55-59", "M_60-64", "M_65-69")
    age_sex_merge_F <- c("F_55-59", "F_60-64", "F_65-69")
    age_sex_merge <- c(age_sex_merge_M, age_sex_merge_F)
    
    pop_merged <- data.frame(age_sex = c("M_55-69", "F_55-69"), 
                             population = c( sum((pop %>% filter(age_sex %in% age_sex_merge_M))$population), 
                                             sum((pop %>% filter(age_sex %in% age_sex_merge_F))$population)))
  }
  else #age_bin = 10
  {
    age_sex_merge <- unique(pop$age_sex) #all age-sex groups will be merged with another group

    pop_merged <- data.frame(age_sex = c("M_15-24", "M_25-34", "M_35-44", "M_45-54", "M_55-69", 
                                         "F_15-24", "F_25-34", "F_35-44", "F_45-54", "F_55-69"), 
                             population = c( sum((pop %>% filter(age_sex %in% c("M_15-19", "M_20-24")))$population), 
                                             sum((pop %>% filter(age_sex %in% c("M_25-29", "M_30-34")))$population), 
                                             sum((pop %>% filter(age_sex %in% c("M_35-39", "M_40-44")))$population), 
                                             sum((pop %>% filter(age_sex %in% c("M_45-49", "M_50-54")))$population), 
                                             sum((pop %>% filter(age_sex %in% c("M_55-59", "M_60-64", "M_65-69")))$population), 
                                             sum((pop %>% filter(age_sex %in% c("F_15-19", "F_20-24")))$population), 
                                             sum((pop %>% filter(age_sex %in% c("F_25-29", "F_30-34")))$population), 
                                             sum((pop %>% filter(age_sex %in% c("F_35-39", "F_40-44")))$population), 
                                             sum((pop %>% filter(age_sex %in% c("F_45-49", "F_50-54")))$population), 
                                             sum((pop %>% filter(age_sex %in% c("F_55-59", "F_60-64", "F_65-69")))$population)
                                          )
                             )
  }
  
  pop2 <- pop %>% filter(! (age_sex %in% c(age_sex_exclude, age_sex_merge)))
  pop2 <- rbind(pop2, pop_merged) %>% arrange(age_sex)
  return(pop2)
}

generate.drho <- function(prev, rho_max, rho_var)
{
  prev <- prev %>% arrange(month, age_sex)
  drho <- matrix(data = 0.0, nrow = length(unique(prev$month)), ncol = length(unique(prev$age_sex)))
  rho <- matrix(data = 0.0, nrow = length(unique(prev$month)), ncol = length(unique(prev$age_sex)))
  
  for(i in 1:ncol(drho))
  {
    drhon <- rdirch(n=1, alpha = rep(1, nrow(drho)))
    rho_max_i <- rho_max*(1 + rho_var*(runif(1) - 0.5))
    drho[,i] <- drhon*rho_max_i
    rho[,i] <- cumsum(drho[,i])
  }
  
  return(list(drho = drho, rho = rho))
}


ppos.to.matrix <- function(ppos, nmonths)
{
  if(length(ppos) < nmonths)
    ppos <- c(ppos, rep(0, nmonths - length(ppos)))
  if(length(ppos) > nmonths)
    ppos <- ppos[1:nmonths]
  
  A <- matrix(data = 0, nrow = length(ppos), ncol = length(ppos))
  for(i in 1:(length(ppos)))
    A[,i] <- c(rep(0, i-1), ppos[1:(length(ppos)-i+1)])
  return(A)
}

prev.to.matrix <- function(prev, field)  #Groups by age_sex and month and converts to matrix
{
  M <- dcast(prev %>% select(month, age_sex, field) %>% arrange(month, age_sex), month ~ age_sex, value.var = field, fill = 0)
  if(ncol(M) == 2)
    M <- as.matrix(M[, 2])
  else
    M <- M[, 2:ncol(M)]
  
  colnames(M) <- NULL
  
  return(M)
}

normalize.drho <- function(drho)
{
  rho_max <- rep(0, nrow(drho))
  drhon <- matrix(data = 0, nrow = nrow(drho), ncol = ncol(drho))
  for(i in 1:ncol(drho)) {
    rho_max[i] <- sum(drho[,i])
    drhon[,i] <- drho[,i]/rho_max[i]
  }
  return(list(rho_max = rho_max, drhon = drhon))
}

simulate.epidemic <- function(prev, drho, ppos, cases_max = 1)
{
  npos.plasma <- 189
  nneg.plasma <- 19
  npos.prepandemic <- 20
  nneg.prepandemic <- 801
  
  if(length(dim(drho)) == 2) {
    dim(drho) <- c(dim(drho), 1)
  }
  
  if("npos" %in% colnames(prev))
    prev <- prev %>% select(-npos)
  
  ntests = prev.to.matrix(prev, "ntests")
  
  nmonths = nrow(ntests)
  ncompartments = ncol(ntests)
 
  A <- ppos.to.matrix(ppos, nmonths)
  
  
  modelCode <- modelCode.BSC()
  
  constants <- list(
    nmonths = nmonths,
    ntests = as.matrix(ntests),
    ncompartments = ncompartments,
    A = A,
    #alpha.pneg = c(0, diff(ppos)),
    alpha = rep(1, nmonths),
    rho_max.lb = 0,
    rho_max.ub = cases_max,
    npos.plasma = npos.plasma,
    nneg.plasma = nneg.plasma,
    npos.prepandemic = npos.prepandemic,
    nneg.prepandemic = nneg.prepandemic
  )
  
  model <- nimbleModel(modelCode, constants = constants, check = FALSE)
  
  previ <- array(0, dim = dim(drho))
  for(i in 1:dim(previ)[3]) {
    if(i %% 100 == 0)
      print(paste0(i, " / ", dim(drho)[3]))
    inits <- normalize.drho(drho[,,i])
    model$setInits(inits)
    model$simulate(c("mu", "npos", "drho", "se", "spb"))
    previ[,,i] <- model$npos
  }
  
  colnames(previ) <- prev$age_sex %>% unique() %>% sort()  
  rownames(previ) <- prev$month %>% unique() %>% sort()  
  df_out <- melt(previ) %>% 
    rename(month = Var1, age_sex = Var2, id_sample = Var3, npos = value) %>%
    left_join(prev, by = c("month", "age_sex"))
  
  return(df_out)
}

summarise.samples <- function(samples)
{
  if(length(dim(samples)) == 3) {
    mg <- c(1,2)
  } else {
    mg <- 2
  }
  
  q1 <- apply(samples, MARGIN = mg, FUN = function(x) quantile(x, 0.025))
  q2 <- apply(samples, MARGIN = mg, FUN = function(x) quantile(x, 0.25))
  q3 <- apply(samples, MARGIN = mg, FUN = function(x) quantile(x, 0.5))
  q4 <- apply(samples, MARGIN = mg, FUN = function(x) quantile(x, 0.75))
  q5 <- apply(samples, MARGIN = mg, FUN = function(x) quantile(x, 0.975))
  
  dfout <- rbind(melt(as.matrix(q1)) %>% mutate(quantile = 0.025), 
                 melt(as.matrix(q2)) %>% mutate(quantile = 0.25), 
                 melt(as.matrix(q3)) %>% mutate(quantile = 0.5), 
                 melt(as.matrix(q4)) %>% mutate(quantile = 0.75), 
                 melt(as.matrix(q5)) %>% mutate(quantile = 0.975)) %>% 
    rename(month = Var1, age_sex = Var2)

  return(dfout)
}

prevalence.blood_donors.measured <- function(prev, population, cases_prior = NA,
                                      nburnin = 10000, niter = 100000)
{
  npos.plasma <- 189
  nneg.plasma <- 19
  npos.prepandemic <- 20
  nneg.prepandemic <- 801
  
  ntests = prev.to.matrix(prev, "ntests")
  npos = prev.to.matrix(prev, "npos")
  
  nmonths = nrow(npos)
  ncompartments = ncol(npos)
  
  if(any(is.na(cases_prior)))
  {
    cases_prior <- rep(1, nmonths)
  }
  else {
    for(i in 1:length(cases_prior))
    {
      cases_prior[i] <- min(1, cases_prior[i])
    }
  }
  
  
  modelCode <- nimbleCode({
    se ~ dbeta(1+npos.plasma, 1+nneg.plasma)
    spb ~ dbeta(1+npos.prepandemic, 1+nneg.prepandemic) #spb = 1 - specificity
    
    for(n in 1:nmonths) {
      for(m in 1:ncompartments) {
        rho[n,m] ~ dunif(0,cases_prior[n])
        mu[n,m] <- se*rho[n,m] + spb*(1 - rho[n,m])
        npos[n,m] ~ dbinom(size = ntests[n,m], prob = mu[n,m])
      }
    }
  })
  
  constants <- list(
    nmonths = nmonths,
    ntests = as.matrix(ntests),
    ncompartments = ncompartments,
    npos.plasma = npos.plasma,
    nneg.plasma = nneg.plasma,
    npos.prepandemic = npos.prepandemic,
    nneg.prepandemic = nneg.prepandemic,
    cases_prior = cases_prior
  )
  
  data <- list(npos = as.matrix(npos))
  inits <- list(rho = as.matrix(npos)/as.matrix(ntests),
                se = 0.9,
                spb = 0.01
  )
  
  
  nimbleMCMC_samples <- nimbleMCMC(code = modelCode, 
                                   constants = constants, 
                                   data = data, 
                                   inits = inits,
                                   nburnin = nburnin, niter = niter, monitors = c("rho", "se", "spb"))
  
  prev_est <- matrix(data = 0, nrow = nrow(nimbleMCMC_samples), ncol = nmonths)
  drho_est <- matrix(data = 0, nrow = nrow(nimbleMCMC_samples), ncol = nmonths)
  
  rho_samples <- list()
  for(m in 1:ncompartments) {
    rho_samples[[m]] = nimbleMCMC_samples[, sapply(1:nmonths, function(x) sprintf(paste0("rho[%d, ", m, "]"), x))]
    prev_est <- prev_est + population[m]*rho_samples[[m]]/sum(population)
  }

  return(list(prev_est = prev_est, rho = rho_samples))
}
