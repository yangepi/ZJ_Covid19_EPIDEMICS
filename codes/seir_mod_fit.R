# ======================================================================
# ======================================================================
library(here)
library(tidyverse)
library(patchwork)
library(nloptr)
library(lubridate)
library(deSolve)
library(pracma)
library(snow)

rm(list = ls())

source(here("plot_theme.R"))
myseed = 525014
# import data
zhejiang_data <- readRDS(here("Zjchn_data")) 
case_data <- readRDS(here("zhejiang_cases")) 
num_split = 5

# ======================================================================
n_samples = 1000

# ======================================================================
temp_pars <- list()
set.seed(seed = myseed)

for (j in 1:n_samples) {
  
  R0s = runif(1, min = 3, max = 4)
  gammas = rgamma(n = 1, scale = 0.1, shape = 20)
  deltas = rgamma(n = 1, scale = 0.2, shape = 20)
  
  ctm_a = runif(1, min = 0.5, max = 1.5)
  ctm_b = runif(1, min = 17,  max = 20)
  
  temp_pars[[j]] <- data.frame(ctm_a = ctm_a, ctm_b = ctm_b, R0s = R0s, gammas = gammas, deltas = deltas, samples = j)
}

all_pars <- bind_rows(temp_pars) %>%
  distinct() %>%
  bind_rows(data.frame(ctm_a = 1, ctm_b = 19, R0s = 3.5, gammas = 2, deltas = 4, samples = 0))

summary(all_pars)

saveRDS(object = all_pars, file = here("all_pars.rds"))
all_pars <- readRDS(file = here("all_pars.rds"))

set.seed(myseed)
all_pars$split <- sample(1:num_split,nrow(all_pars),replace = T)
all_pars$ids <- 1:nrow(all_pars)
table(all_pars$split)

# ======================================================================
cluster_fun <- function(i = NA){
  
  # import data
  zhejiang_data <- readRDS(here("Zjchn_data")) 
  # ode format
  age_labels <<-  c("[0,5)","[5,10)","[10,15)","[15,20)","[20,25)","[25,30)","[30,35)","[35,40)","[40,45)","[45,50)","[50,55)","[55,60)","[60,65)","65+")
  zhejiang_data$pop <- zhejiang_data$pop %>%
    mutate(p_pop = round(population/sum(population),4))
  N <- zhejiang_data$pop[,2] %>% as.matrix()  
  zhejiang_cases <- zhejiang_data$cases 
  # CTM age label
  rownames(zhejiang_data$ctm_base) <- age_labels
  colnames(zhejiang_data$ctm_base) <- age_labels
  rownames(zhejiang_data$ctm_outbreak) <- age_labels
  colnames(zhejiang_data$ctm_outbreak) <- age_labels
  zhejiang_data <<- zhejiang_data 
  # ===============
  source(here("ode_functions.R"))
  source(here("get_beta.R"))
  
  # i=1 # from 1 to all_pars 
  beta_lpi <- get_beta(C = zhejiang_data$ctm_base,
                       propns = as.matrix(zhejiang_data$pop[,"p_pop"]),
                       inf_period = all_pars$gammas[i] + all_pars$deltas[i],
                       r0 = all_pars$R0s[i]) %>%
    log()
  
  # creat ith fitting parameter setting
  ranglb = c(ctm_a = all_pars$ctm_a[i], ctm_b = all_pars$ctm_b[i], beta_I = beta_lpi, par_epsilon = 1,  fp = 0.1, C0 = 0.9, C1 = 0.9, par_gamma = all_pars$gammas[i], par_delta = all_pars$deltas[i])
  x0     = c(ctm_a = all_pars$ctm_a[i], ctm_b = all_pars$ctm_b[i], beta_I = beta_lpi, par_epsilon = 5,  fp = 0.4, C0 = 1,   C1 = 1,   par_gamma = all_pars$gammas[i], par_delta = all_pars$deltas[i])
  rangub = c(ctm_a = all_pars$ctm_a[i], ctm_b = all_pars$ctm_b[i], beta_I = beta_lpi, par_epsilon = 15, fp = 0.9, C0 = 1.3, C1 = 1.3, par_gamma = all_pars$gammas[i], par_delta = all_pars$deltas[i])
  
  mod_fit_data = zhejiang_cases
  test_dur = max(mod_fit_data$days_post)
  
  # fit the model with SEIR
  set.seed(myseed)
  mod_fit <- zj_fit(dat = mod_fit_data, ranglb=ranglb, x0=x0, rangub=rangub, test_dur = test_dur, N = N)
  
  return(list(mod = mod_fit, pars = x0, lps = i))
}

# ====================================================================== 
# parallel
cl <- makeCluster(25, type = "SOCK") 
clusterEvalQ(cl, c(library(tidyverse), library(patchwork), library(nloptr), library(deSolve), library(pracma), library(here), library(lubridate)))
clusterExport(cl, c("zhejiang_data","all_pars","myseed"))

for (k in 1:num_split) {
  inum <- all_pars %>% filter(split == k) %>% pull(ids)
  savepath <- paste0("mod_fit_",k,".rds")
  a1 <- Sys.time()
  cat("\nStart to fit ",k, "(",a1,")")
  
  # split the workload by each time 
  zj_fit_all <- clusterApplyLB(cl, inum, cluster_fun)
  saveRDS(zj_fit_all,file = here(savepath))
  gc()
  
  cat("; Finished",k," run, ", k*100/num_split, "%; took ", Sys.time()-a1)
}

stopCluster(cl)

