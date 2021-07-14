# SEIR model
rm(list = ls())

library(here)
library(tidyverse)
library(patchwork)
library(nloptr)
library(lubridate)
library(deSolve)
library(pracma)
library(snow)

rm(list = ls())
zhejiang_data <- readRDS(here("Zjchn_data")) 
combine_zj1kfit <- readRDS(here("combine_zj1kfit.rds"))
combine_zjfxfit <- readRDS(here("combine_zjfxfit.rds"))
all1k1_zjfits <- append(combine_zjfxfit,combine_zj1kfit)
pars1k_tab <- list()

for (i in 1:length(combine_zj1kfit)) {
  model <- combine_zj1kfit[[i]]
  model_pars <- model$mod$nloptr$solution
  
  C0 = as.numeric(model_pars[6])
  C1 = as.numeric(model_pars[7])
  
  pars1k_tab[[i]] <- data.frame(mod = i,
                                par_epsilon = as.numeric(model_pars[4]),
                                fp = as.numeric(model_pars[5]),
                                C0scale = C0,
                                C1scale = C1,
                                ctm_reduced = ctm_reduced)
}


pars1k_tab <- bind_rows(pars1k_tab)
cip <- rethinking::PI(pars1k_tab$fp, prob = 0.95)
cie <- rethinking::PI(pars1k_tab$par_epsilon, prob = 0.95)
cic0 <- rethinking::PI(pars1k_tab$C0scale, prob = 0.95)
cic1 <- rethinking::PI(pars1k_tab$C1scale, prob = 0.95)

pars_ci <- data.frame(pars = c("fp","par_epsilon", "C0scale", "C1scale"),
                      mean = c(mean(pars1k_tab$fp),mean(pars1k_tab$par_epsilon),mean(pars1k_tab$C0scale),mean(pars1k_tab$C1scale)),
                      hdilo     = c(cip[1],cie[1],cic0[1],cic1[1]),
                      hdiup     = c(cip[2],cie[2],cic0[2],cic1[2])) %>% 
  mutate_if(is.numeric,round,4) %>% 
  mutate(labels = paste0(round(mean,2),"\n95%CI (",round(hdilo,2),", ",round(hdiup,2),")"))
