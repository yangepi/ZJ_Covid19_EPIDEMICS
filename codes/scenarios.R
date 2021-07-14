
library(here)
library(tidyverse)
library(patchwork)
library(nloptr)
library(lubridate)
library(deSolve)
library(pracma)
library(snow)
library(ggthemes)
selfcolors <- colorblind_pal()(8)

rm(list = ls())

source(here("plot_theme.R"))
source(here("ode_functions.R"))
myseed = 525014

# import data
zhejiang_data <- readRDS(here("Zjchn_data")) 
combine_zj1kfit <- readRDS(here("combine_zj1kfit.rds"))
N <- zhejiang_data$pop[,2] %>% as.matrix()  

# ====================================================================== 
cluster_fun <- function(i = NA){
  
  source(here("ode_functions.R"))
  
  age_labels <<- c("[0,5)","[5,10)","[10,15)","[15,20)","[20,25)","[25,30)","[30,35)","[35,40)","[40,45)","[45,50)","[50,55)","[55,60)","[60,65)","65+")
  seir_op <- list()
  scenarios_op <- list()
  
  model <- combine_zj1kfit[[i]]
  model_pars <- model$mod$nloptr$solution
  model_id <-  model$lps
  
  ctm_a = as.numeric(model_pars[1])
  ctm_b = as.numeric(model_pars[2])
  beta_I  = as.numeric(model_pars[3])
  C0scale = as.numeric(model_pars[6])
  C1scale = as.numeric(model_pars[7])
  par_gamma = as.numeric(model_pars[8])
  par_delta = as.numeric(model_pars[9])
  
  epsilon_value = as.numeric(model_pars[4]) * c(1,1.2,0.8)
  fp_value = as.numeric(model_pars[5]) * c(1,1.2,0.8)
  
  scenarios = data.frame(epsilon_names = c("Fitted value","Fitted value","Fitted value","Increased by 20%","Increased by 20%","Increased by 20%","Decreased by 20%","Decreased by 20%","Decreased by 20%"),
                         epsilon = c(epsilon_value[1],epsilon_value[1],epsilon_value[1],epsilon_value[2],epsilon_value[2],epsilon_value[2],epsilon_value[3],epsilon_value[3],epsilon_value[3]),
                         fp_names = c("Fitted value","Increased by 20%","Decreased by 20%","Fitted value","Increased by 20%","Decreased by 20%","Fitted value","Increased by 20%","Decreased by 20%"),
                         fp = c(fp_value[1],fp_value[2],fp_value[3],fp_value[1],fp_value[2],fp_value[3],fp_value[1],fp_value[2],fp_value[3]))
  
  for (j in 1:nrow(scenarios)) {
    temp <- epi_ode_seir(C0 = zhejiang_data$ctm_base*C0scale, 
                         C1 = zhejiang_data$ctm_outbreak*C1scale, 
                         Ns = N, 
                         # because here is the prediction, so we used the total time of real data 46 days
                         ts=seq(0,46,by=1),
                         ctm_a = ctm_a, 
                         ctm_b = ctm_b,
                         beta_I = exp(beta_I),
                         par_gamma = par_gamma, 
                         par_delta = par_delta, 
                         par_epsilon = scenarios$epsilon[j], 
                         fp = scenarios$fp[j])
    
    seir_op[[j]] <- temp$daily_newcases %>% 
      mutate(days = min(zhejiang_data$cases$date) + c(0:46),
             model_id = model_id, 
             epsilon_names = scenarios$epsilon_names[j], par_epsilon = scenarios$epsilon[j],
             fp_names = scenarios$fp_names[j], fp = scenarios$fp[j])
    
  }
  op <- bind_rows(seir_op)
  return(op)
}

# ====================================================================== 
# parallel
cl <- makeCluster(25, type = "SOCK") 
clusterEvalQ(cl, c(library(tidyverse), library(nloptr), library(deSolve), library(pracma), library(here), library(lubridate)))
clusterExport(cl, c("zhejiang_data","combine_zj1kfit","N"))
scenarios <- clusterApplyLB(cl, 1:1000, cluster_fun)
stopCluster(cl)
# saveRDS(scenarios,here("scenarios.rds"))
# scenarios <- readRDS(here("scenarios.rds"))
scenarios_tab <- bind_rows(scenarios)

scenarios_smr <- scenarios_tab %>% 
  group_by(days,epsilon_names,fp_names) %>% 
  summarise(case_mean = mean(dailynew,na.rm=T),
            case_lo = quantile(dailynew,0.025,na.rm=T),
            case_up = quantile(dailynew,0.975,na.rm=T)) %>% 
  mutate(epsilon_names = factor(epsilon_names,levels = c("Decreased by 20%","Fitted value","Increased by 20%")),
         fp_names = factor(fp_names,levels = c("Decreased by 20%","Fitted value","Increased by 20%"))) %>% 
  ungroup()

fig3 <- ggplot() + 
  geom_line(data = scenarios_smr, aes(x = days, y = case_mean, linetype = factor(fp_names)), size = 0.5) +
  geom_ribbon(data = scenarios_smr, aes(x=days, ymin = case_lo, ymax = case_up, group = factor(fp_names), fill = factor(fp_names)), 
              linetype=2, alpha = 0.2, show.legend = F) +
  facet_grid(~epsilon_names) +
  scale_linetype_manual(values = c("Decreased by 20%" = 2, "Fitted value" = 1, "Increased by 20%" = 3),
                        name = "Contacts tracing: p",
                        guide = guide_legend(override.aes = list(size = 2))) +
  labs(x = "Date", y = "Counts of cases", title = bquote(atop("Isolation speed: " ~ epsilon))) +
  scale_x_date(date_labels = "%b %d", date_breaks = "5 day") +
  coord_cartesian(ylim = c(0,150)) +
  scale_fill_colorblind() +
  theme_bw() +
  plot_theme +
  guides(linetype = guide_legend(override.aes = list(size = 1))) +
  theme(strip.background =element_rect(fill="white", colour="white"),
        strip.text = element_text(size = 9),
        legend.position = "bottom",
        legend.margin = margin(-10,0,0,0),
        plot.title = element_text(size = 10, margin=margin(0,0,-10,0), hjust = 0.5, face = "bold"),
        axis.text.x=element_text(angle=60, hjust=1),
        legend.key.width = unit(1.2,"cm"),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank())  
  
fig3
ggsave(file = here("figure3.png"), fig3, dpi=dpiset, units="in", width=7, height=5)
