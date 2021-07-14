# ======================================================================
# ======================================================================
library(here)
library(tidyverse)
library(patchwork)
library(lubridate)
library(deSolve)
library(pracma)

rm(list = ls())

source(here("plot_theme.R"))
myseed = 525014
# import data
zhejiang_data <- readRDS(here("Zjchn_data")) 
case_data <- readRDS(here("zhejiang_cases")) 
all_pars <- readRDS(file = here("all_pars.rds"))

num_split = 5

# ====================================================================== 
# combine all fit objects
combine_zj_fit <- list()
pars_x0_list <- list()
pars_xfit_list <- list()

# read fitted models
for (k in 1:num_split) {
  savepath <- paste0("mod_fit_",k,".rds")
  lista <- readRDS(file = here(savepath))
  if(k==1){combine_zj_fit <- lista}
  if(k> 1){combine_zj_fit <- append(combine_zj_fit, lista)}
  gc()
}

combine_zj1kfit <- combine_zj_fit[-791]
combine_zjfxfit <- combine_zj_fit[791]
saveRDS(combine_zj1kfit,here("combine_zj1kfit.rds"))
saveRDS(combine_zjfxfit,here("combine_zjfxfit.rds"))

for (k in 1:length(combine_zj1kfit)) {
  pars_x0_list[[k]]   <- as.data.frame(t(combine_zj1kfit[[k]]$pars)) 
  pars_xfit_list[[k]] <- as.data.frame(t(combine_zj1kfit[[k]]$mod$nloptr$solution)) %>% 
    mutate(ttcase = combine_zj1kfit[[k]]$mod$fit_ts$cumulative_cases_sum)
  # cat("\nFinished",k," run, ", k*100/length(combine_zj1kfit), "%")
}

pars_x0 <- bind_rows(pars_x0_list)
pars_xfit <- bind_rows(pars_xfit_list)
colnames(pars_xfit)[1:9] <- colnames(pars_x0)

cip <- rethinking::PI(pars_xfit$fp, prob = 0.95)
cie <- rethinking::PI(pars_xfit$par_epsilon, prob = 0.95)
cic0 <- rethinking::PI(pars_xfit$C0, prob = 0.95)
cic1 <- rethinking::PI(pars_xfit$C1, prob = 0.95)

pars_ci <- data.frame(pars = c("fp","par_epsilon", "C0scale", "C1scale"),
                      mean = c(mean(pars_xfit$fp),mean(pars_xfit$par_epsilon),mean(pars_xfit$C0),mean(pars_xfit$C1)),
                      hdilo     = c(cip[1],cie[1],cic0[1],cic1[1]),
                      hdiup     = c(cip[2],cie[2],cic0[2],cic1[2])) %>% 
  mutate_if(is.numeric,round,4) %>% 
  mutate(notes = c("\n\U00070 = ","\n\U003B5 = ","\n\U003BC = ","\n\U003BC = "),
         labels = paste0(notes,round(mean,2),"\n95%CI (",round(hdilo,2),", ",round(hdiup,2),")"))

pars_ci


# ====================================================================== 
# create prediction intervals
# =================
# social contact interval
test_dur = max(zhejiang_data$cases$days_post)
changes <- matrix(nrow = length(combine_zj1kfit), ncol = test_dur+1, data = NA)

for (i in 1:length(combine_zj1kfit)) {
  changes[i,] <- sigmoid(0:test_dur-pars_xfit$ctm_b[i], a = pars_xfit$ctm_a[i], b = 0) %>% as.numeric()
}

# plot(x = , y = )
ctm_tab <- tibble(date = min(zhejiang_data$cases$date) + c(0:test_dur), 
                  ctm_mean = 1-apply(changes,2,mean,na.rm=T),
                  ctm_lo =   1-apply(changes,2,quantile,0.025,na.rm=T),
                  ctm_up =   1-apply(changes,2,quantile,0.975,na.rm=T))

ctm_plot <- ggplot(ctm_tab,aes(x = date, y = ctm_mean)) +
  geom_ribbon(aes( ymin = ctm_lo, ymax = ctm_up), fill = "grey", alpha = 0.8) +
  geom_line(size = 0.5) +
  geom_vline(xintercept = as_date("2020-01-23"), linetype = 2, col = "red", size = 0.5) +
  geom_vline(xintercept = as_date("2020-02-01"), linetype = 2, col = "red", size = 0.5) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1),labels = c("Outbreak\nstrength","Middle","Baseline\nstrength")) +
  scale_x_date(date_labels = "%b %d", date_breaks = "3 day") +
  labs(x = "Date", y="Contact matrix") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=60, hjust=1),
        panel.border = element_blank(),
        panel.grid.major = element_blank()) +
  plot_theme

ctm_plot

# =================
# daily new cases interval
dnewcases <- matrix(nrow = length(combine_zj1kfit), ncol = test_dur+1, data = NA)

for (i in 1:length(combine_zj1kfit)) {
  dnewcases[i,] <- combine_zj1kfit[[i]]$mod$fit_ts$daily_newcases$dailynew
}

case_tab <- tibble(date = min(zhejiang_data$cases$date) + c(0:test_dur), 
                   case_mean = apply(dnewcases,2,mean,na.rm=T),
                   case_lo = apply(dnewcases,2,quantile,0.025,na.rm=T),
                   case_up = apply(dnewcases,2,quantile,0.975,na.rm=T))

p_plot <- data.frame(pars = "p", value = pars_xfit$fp)  %>% 
  ggplot(aes(x=value)) + 
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity", binwidth = 0.05, color="black", fill="white") +
  geom_density(adjust=2,size = 0.5,col=c("blue")) +
  labs(x = paste0("Contact tracing \nproportion",pars_ci$labels[1]), y = "Density") + # , title = expression(paste("Density plot of p"))) +
  scale_x_continuous(breaks = seq(0.1,0.9,0.2)) +
  plot_theme +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_text(size=9),
        axis.title.x.bottom = element_text(vjust = 12))

e_plot <- data.frame(pars = "e", value = pars_xfit$par_epsilon)  %>% 
  ggplot(aes(x=value)) + 
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity", binwidth = 1, color="black", fill="white") +
  geom_density(adjust=2,size = 0.5,col=c("blue")) +
  labs(x = paste0("Isolation speed (days)",pars_ci$labels[2]), y = "Density") + 
  plot_theme +
  scale_x_continuous(breaks = seq(1,15,3)) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_text(size=9),
        axis.title.x.bottom = element_text(vjust = 12))

ode_plot <- ggplot()+ 
  geom_vline(xintercept = as_date("2020-01-23"), linetype = 2, col = "red", size = 0.5) +
  geom_vline(xintercept = as_date("2020-02-01"), linetype = 2, col = "red", size = 0.5) +
  geom_col(data = zhejiang_data$cases, mapping=aes(x=date, y=dailynewcases, fill = "data"), alpha = 0.5) +
  geom_line(data = case_tab, aes(x=date,y=case_mean,col = "blue"), size = 0.8) +
  geom_ribbon(data = case_tab, aes(x=date,y=case_mean, ymin = case_lo, ymax = case_up), fill = "blue", linetype=2, alpha = 0.2) +
  scale_color_manual(values = c("blue"), labels = c("Model"), name = " ") +
  scale_fill_manual(values = c("grey"), labels = c("Observed cases"), name = " ") +
  labs(x = "Date", y = "Counts of cases") +
  scale_x_date(date_labels = "%b %d", date_breaks = "3 day") +
  theme_bw() +
  plot_theme +
  theme(legend.position = c(0.85,0.75),
        legend.key.size = unit(0.5, "cm"),
        legend.spacing.y = unit(-0.1, "cm"),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.text.x=element_text(angle=60, hjust=1))

ode_plot

fig2 <- (ctm_plot + theme(axis.title.x = element_blank()) + e_plot + p_plot +  plot_layout(widths = c(3, 1.1, 1.1))) /
  (ode_plot)  + plot_layout(heights = c(3, 6)) + plot_annotation(tag_levels = 'A')

fig2
pars_ci

ggsave(file = here("figure2.png"), fig2, dpi=dpiset, units="in", width=8, height=5)