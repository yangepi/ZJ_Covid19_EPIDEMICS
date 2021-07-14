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

source(here("plot_theme.R"))

myseed = 525014

combine_zj1kfit <- readRDS(here("combine_zj1kfit.rds"))

# import data
zhejiang_data <- readRDS(here("Zjchn_data")) 
dur_dat = max(zhejiang_data$cases$days_post)

# ====================================================================== 
cluster_fun <- function(i = NA){
  
  test_dur = 365
  
  source(here("reopen_odes_functions.R"))
  zhejiang_data$pop <- zhejiang_data$pop %>% 
    mutate(p_pop = round(population/sum(population),4))
  N <- zhejiang_data$pop[,2] %>% as.matrix() 
  zhejiang_cases <- zhejiang_data$cases 
  
  solutions <- combine_zj1kfit[[i]]$mod$nloptr$solution
  
  ctm_a_fit = as.numeric(solutions[1])
  ctm_b_fit = as.numeric(solutions[2])
  beta_I_fit  = as.numeric(solutions[3])
  epsilon_fit  = as.numeric(solutions[4])
  fp_fit = as.numeric(solutions[5])
  ctm_lo_fit = as.numeric(solutions[6])
  ctm_up_fit = as.numeric(solutions[7])
  par_gamma_fit = as.numeric(solutions[8])
  par_delta_fit = as.numeric(solutions[9])
  
  matrix_pool <- data.frame(mod_id = paste0("mod",rep(i,3)),
                            scenarios = paste0("scenarios",1:3),
                            speed0 = c(epsilon_fit, epsilon_fit, epsilon_fit),
                            speed1 = c(epsilon_fit,3,2),
                            speed3 = c(epsilon_fit,1.5,1),
                            tp1 = c(fp_fit,fp_fit,fp_fit),
                            tp2 = c(fp_fit,0.4,0.5),
                            tp3 = c(fp_fit,0.8,0.8),
                            addcasesp = NA,
                            addcases = NA,
                            ifc_p = NA)
  
  ts <- list()
  
  for (j in 1:3) {
    epsilon_list_i <- matrix_pool[j, 3:5] %>% as.numeric()
    fp_list_i <- matrix_pool[j, 6:8] %>% as.numeric()
    
    # reopening setting
    res <- epi_ode_seir_reopen(C0 = zhejiang_data$ctm_base*ctm_lo_fit, 
                               C1 = zhejiang_data$ctm_outbreak*ctm_up_fit, 
                               Ns = N, 
                               ts=seq(0,test_dur,by=1),
                               ctm_a = ctm_a_fit, 
                               ctm_b = ctm_b_fit,
                               beta_I = exp(beta_I_fit), 
                               par_gamma = par_gamma_fit, 
                               par_delta = par_delta_fit, 
                               par_epsilon_list = epsilon_list_i, 
                               fp_list = fp_list_i)
    
    res$daily_newcases$date <- min(zhejiang_cases$date) + c(0:test_dur)
    matrix_pool$addcasesp[j] <- (res$daily_newcases$H[test_dur+1] - res$daily_newcases$H[47])/res$daily_newcases$H[47]
    matrix_pool$addcases[j]  <- (res$daily_newcases$H[test_dur+1] - res$daily_newcases$H[47])
    matrix_pool$ifc_p[j]  <- round(res$daily_newcases$H[test_dur+1] / sum(zhejiang_data$pop$population),5)
    
    ts[[j]] <- res$daily_newcases %>% mutate(mod_id = paste0("mod",i),scenarios = paste0("Scenarios ",j),times = paste0("days",c(0:test_dur)))
  }
  
  op_smr <- matrix_pool
  op_tsd <- bind_rows(ts)
  
  return(list(smr = op_smr, ts = op_tsd))
}

# ====================================================================== 
# parallel
cl <- makeCluster(20, type = "SOCK") 
clusterEvalQ(cl, c(library(tidyverse), library(nloptr), library(deSolve), library(pracma), library(here), library(lubridate)))
clusterExport(cl, c("zhejiang_data","combine_zj1kfit"))
reopens <- clusterApplyLB(cl, 1:1000, cluster_fun)
stopCluster(cl)
# saveRDS(reopens,here("reopen_simulations.rds"))
# reopens <- readRDS(here("reopen_simulations.rds"))
smr <- list()
ts <- list()

for (i in 1:length(reopens)) {
  ts[[i]] <- reopens[[i]]$ts
  smr[[i]] <- reopens[[i]]$smr
}

ts_tab <- bind_rows(ts)
smr_tab <- bind_rows(smr)

ts_tab_smr <- ts_tab %>% 
  group_by(date,scenarios) %>% 
  summarise(case_mean = mean(dailynew,na.rm=T),
            case_lo = quantile(dailynew,0.025,na.rm=T),
            case_up = quantile(dailynew,0.975,na.rm=T)) %>% 
  mutate(case_mean2 = if_else(case_mean<=1,1,case_mean),
         case_lo2 = if_else(case_lo<=1,1,case_lo),
         case_up2 = if_else(case_up<=1,1,case_up))

facet_text <- data.frame(scenarios = paste0("Scenarios ",1:3),
                         lkde = c("\U003B5","\U003B5","\U003B5"), 
                         phase1e = c("\U003B5","3 days","2 days"), 
                         phase2e = c("\U003B5","1.5 days","1 day"),
                         lkdp = c("\U00070","\U00070","\U00070"),
                         phase1p = c("\U00070","40%","50%"), 
                         phase2p = c("\U00070","80%","80%")) %>% 
  mutate(labels = paste0(lkde," \U02192 ",phase1e," \U02192 ",phase2e,"\n",lkdp," \U02192 ",phase1p," \U02192 ",phase2p))

reopen_plot <- ggplot()+ 
  geom_line(data = ts_tab_smr, aes(x = date, y = case_mean), size = 0.5) +
  geom_ribbon(data = ts_tab_smr, aes(x=date, ymin = case_lo, ymax = case_up), alpha = 0.2) +
  facet_grid(scenarios~.,scales = "free", labeller = labeller(scenarios = c("Scenarios 1" = "Scenario 1", "Scenarios 2" = "Scenario 2", "Scenarios 3" = "Scenario 3"))) +
  scale_x_date(date_labels = "%b %d", date_breaks = "30 day", expand = c(0, 0)) +
  geom_vline(xintercept = min(ts_tab_smr$date) + 90, linetype = 2, size = 0.3) +
  geom_vline(xintercept = min(ts_tab_smr$date) + 180, linetype = 2, size = 0.3) +
  annotate("rect", xmin = min(ts_tab_smr$date) + 90, xmax = min(ts_tab_smr$date) + 180, ymin = -Inf, ymax = Inf, fill = "#afe6bf", alpha = .15) +
  annotate("rect", xmin = min(ts_tab_smr$date) + 180, xmax = min(ts_tab_smr$date) + 365, ymin = -Inf, ymax = Inf, fill = "#bf7fbf", alpha = .15) +
  geom_text(data = facet_text,aes(x =  min(ts_tab_smr$date) + 27, y = Inf, label = labels), hjust = 0, vjust = 1.5, size = 3) +
  theme_bw() +
  plot_theme +
  # scale_y_log10() +
  theme(strip.background =element_rect(fill="white", colour="white"),
        legend.position = "bottom",
        legend.margin = margin(-10,0,0,0),
        plot.title = element_text(size = 10, margin=margin(0,0,-10,0), hjust = 0.5),
        axis.text.x=element_text(angle=60, hjust=1),
        legend.key.width = unit(1.2,"cm"),
        panel.grid.major = element_blank(),
        panel.border = element_blank())
  
reopen_plot
ggsave(file = here("figure5-normal-y.png"), reopen_plot, dpi=dpiset, units="in", width=8, height=5)

scenario1 <- ts_tab_smr %>% filter(scenarios == "Scenarios 1")
sum(scenario1$case_mean,na.rm = T) / sum(zhejiang_data$pop$population)

scenario2 <- ts_tab_smr %>% filter(scenarios == "Scenarios 2")
sum(scenario2$case_mean,na.rm = T)

scenario3 <- ts_tab_smr %>% filter(scenarios == "Scenarios 3")
sum(scenario3$case_mean,na.rm = T)

reopen_plot <- ggplot()+ 
  geom_line(data = ts_tab_smr, aes(x = date, y = case_mean2), size = 0.5) +
  geom_ribbon(data = ts_tab_smr, aes(x=date, ymin = case_lo2, ymax = case_up2), alpha = 0.2) +
  facet_grid(scenarios~.,labeller = labeller(scenarios = c("Scenarios 1" = "Scenario 1", "Scenarios 2" = "Scenario 2", "Scenarios 3" = "Scenario 3"))) +
  scale_x_date(date_labels = "%b %d", date_breaks = "30 day", expand = c(0, 0)) +
  geom_vline(xintercept = min(ts_tab_smr$date) + 90, linetype = 2, size = 0.3) +
  geom_vline(xintercept = min(ts_tab_smr$date) + 180, linetype = 2, size = 0.3) +
  theme_bw() +
  scale_y_log10(labels = function(x) format(x, scientific = FALSE), limits = c(1,5e6)) +
  labs(x = "Date", y = "Counts of cases (log scale)") +
  annotate("rect", xmin = min(ts_tab_smr$date) + 90, xmax = min(ts_tab_smr$date) + 180, ymin = 0, ymax = Inf, fill = "#afe6bf", alpha = .15) +
  annotate("rect", xmin = min(ts_tab_smr$date) + 180, xmax = min(ts_tab_smr$date) + 365, ymin = 0, ymax = Inf, fill = "#bf7fbf", alpha = .15) +
  geom_text(data = facet_text,aes(x =  min(ts_tab_smr$date) + 35, y = Inf, label = labels), hjust = 0, vjust = 1.5, size = 3) +
  plot_theme +
  theme(strip.background =element_rect(fill="white", colour="white"),
        legend.position = "bottom",
        legend.margin = margin(-10,0,0,0),
        plot.title = element_text(size = 10, margin=margin(0,0,-10,0), hjust = 0.5),
        axis.text.x=element_text(angle=60, hjust=1),
        legend.key.width = unit(1.2,"cm"),
        panel.grid.major = element_blank(),
        panel.border = element_blank())

reopen_plot
ggsave(file = here("figure5-log-y.png"), reopen_plot, dpi=dpiset, units="in", width=8, height=5)
