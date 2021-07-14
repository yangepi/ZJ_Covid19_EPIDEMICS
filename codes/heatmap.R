library(here)
library(tidyverse)
library(patchwork)
library(nloptr)
library(lubridate)
library(deSolve)
library(pracma)
library(snow)
library(wesanderson)
library(viridis)

rm(list = ls())

source(here("plot_theme.R"))
source(here("ode_functions.R"))
myseed = 525014

# import data
zhejiang_data <- readRDS(here("Zjchn_data")) 
combine_zjfxfit <- readRDS(here("combine_zjfxfit.rds"))
N <- zhejiang_data$pop[,2] %>% as.matrix()  

# ====================================================================== 
cluster_fun <- function(i = NA){
  # browser()
  source(here("ode_functions.R"))
  
  age_labels <<- c("[0,5)","[5,10)","[10,15)","[15,20)","[20,25)","[25,30)","[30,35)","[35,40)","[40,45)","[45,50)","[50,55)","[55,60)","[60,65)","65+")
  seir_op <- list()
  scenarios_op <- list()
  
  model <- combine_zjfxfit[[1]]
  model_pars <- model$mod$nloptr$solution
  model_id <-  model$lps
  model_preidct_cases <- model$mod$fit_ts$cumulative_cases_sum
  model_preidct_cases <- 664
  
  ctm_a = as.numeric(model_pars[1])
  ctm_b = as.numeric(model_pars[2])
  beta_I  = as.numeric(model_pars[3])
  C0scale = as.numeric(model_pars[6])
  C1scale = as.numeric(model_pars[7])
  par_gamma = as.numeric(model_pars[8])
  par_delta = as.numeric(model_pars[9])
  
  ctm_reduced <- max(Re(eigen(C1scale*zhejiang_data$ctm_outbreak)$values))/max(Re(eigen(C0scale*zhejiang_data$ctm_base)$values))
  ctm_scale_pool <- c(0.1,0.3,0.5,0.7)
  epsilon_p_pool = crossing(ctm_scale = ctm_scale_pool, epsilon_pool = seq(2,8,0.1), fp_pool = seq(0.2,0.6,0.01))

  outbreak_ctm <- C0scale*zhejiang_data$ctm_base * epsilon_p_pool$ctm_scale[i]
  
  temp <- epi_ode_seir(C0 = zhejiang_data$ctm_base*C0scale, 
                       C1 = outbreak_ctm, 
                       Ns = N, 
                       ts=seq(0,46,by=1),
                       ctm_a = ctm_a, 
                       ctm_b = ctm_b,
                       beta_I = exp(beta_I),
                       par_gamma = par_gamma, 
                       par_delta = par_delta, 
                       par_epsilon = epsilon_p_pool$epsilon_pool[i], 
                       fp = epsilon_p_pool$fp_pool[i])
  
  ttcases <- tibble(fp = epsilon_p_pool$fp_pool[i], epsilon = epsilon_p_pool$epsilon_pool[i], ctm_scale = epsilon_p_pool$ctm_scale[i], 
                    ttcases = temp$cumulative_cases_sum, 
                    casesratio = round(ttcases/model_preidct_cases,2)) 
  
  return(ttcases)
}

# ====================================================================== 
# parallel
cl <- makeCluster(20, type = "SOCK") 
clusterEvalQ(cl, c(library(tidyverse), library(nloptr), library(deSolve), library(pracma), library(here), library(lubridate)))
clusterExport(cl, c("zhejiang_data","combine_zjfxfit","N"))
heatmap <- clusterApplyLB(cl, 1:10004, cluster_fun)
stopCluster(cl)
# saveRDS(heatmap,here("heatmap.rds"))

heatmap_dat <- bind_rows(heatmap) %>% 
  mutate(ctm_labels = paste0("Outbreak CNT = ", ctm_scale*100, "% of baseline")) %>% 
  mutate(ctm_labels = factor(ctm_labels,c("Outbreak CNT = 10% of baseline",
                                          "Outbreak CNT = 30% of baseline",
                                          "Outbreak CNT = 50% of baseline",
                                          "Outbreak CNT = 70% of baseline")))

heatmap_plot <- ggplot(heatmap_dat, aes(x = epsilon, y = fp, z = casesratio))+
  geom_contour_filled(breaks = c(0,2,10,100,500,1000), alpha = 0.7) +
  scale_fill_brewer(palette = "Spectral", direction = -1) +
  facet_wrap(~ctm_labels) +
  scale_y_continuous(breaks = seq(0.2,0.6,0.1), labels = paste0(100*seq(0.2,0.6,0.1),"%")) +
  scale_x_continuous(breaks = seq(2,8), labels = paste0(seq(2,8))) +
  labs(x = expression(paste(epsilon, ": case isolation speed (days)")),  
       y = expression(paste("p: contact tracing proportion")),
       fill = "Ratio (predicted / observed # of cases)") +
  plot_theme + 
  theme_bw() +
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5)) +
  theme(strip.background =element_rect(fill="white", colour="white"),
        legend.key.width = unit(1, "cm"),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
  
heatmap_plot
ggsave(file = here("results/figure4.png"), heatmap_plot, dpi=dpiset, units="in", width=6, height=6) 
 