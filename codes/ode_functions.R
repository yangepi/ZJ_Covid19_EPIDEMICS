
setup_C <- function(C1, Ns, beta_scales=NA){
  # browser()
  Nimmunity <- ncol(Ns)
  Nage <- nrow(Ns)
  M <-  kron(C1, ones(Nimmunity,Nimmunity)) 
  propns <- Ns/rowSums(Ns) 
  propns_1 <- repmat(matrix(rep(propns, each=Nimmunity),ncol = Nimmunity,byrow=FALSE),n=1,m=Nage)
  C <- M*propns_1
  beta_scales_C <- repmat(matrix(rep(beta_scales, each=Nimmunity),ncol=Nimmunity*Nage,byrow=TRUE),n=Nage*Nimmunity,m=1)
  C <- C*beta_scales_C
  N_long <- c(t(Ns))
  C <- t(apply(C, 1, function(x) x/N_long))
  C[!is.finite(C)] <- 0
  return(C)
}

# basic ode function
ode_seir <- function(t, y, par, Nage, N, 
                     C1st, C2st, ctm_a, ctm_b, 
                     beta_I, par_gamma, par_delta, par_epsilon, fp){
  
  changes <- sigmoid(t, a = ctm_a, b = ctm_b) 
  ctm_fun <- C1st - (C1st - C2st)*changes
  beta_scales <- c(0.4,0.4,0.4,0.4,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8)   
  C = setup_C(ctm_fun, N, beta_scales)
  
  # 8 compartments: SEPCR + Qp Qc +H
  seir <- matrix(y, ncol=Nage, nrow=8) 
  # P and C: pre-clinic & clinic cases had transmissibility
  
  dS <- - (seir[1,]* seir[3,] %*% t(C)) * beta_I - (seir[1,]* seir[4,] %*% t(C)) * beta_I
  # ---
  dE <- - dS - seir[2,]/par_gamma
  # ---
  dP <- (1-fp)*seir[2,]/par_gamma - seir[3,]/par_delta
  # ---
  dC <- seir[3,]/par_delta - seir[4,]/par_epsilon
  
  # ---
  dQp <- (fp)*seir[2,]/par_gamma - seir[5,]/par_delta
  dQc <- seir[5,]/par_delta - seir[6,]/20
  
  # ---
  dR <- seir[4,]/par_epsilon + seir[6,]/20
  
  # ---
  # for daily symptom cases fitting (Ci)
  dH <- seir[3,]/par_delta + (fp)*seir[2,]/par_gamma # cumulate case having symptom   
  
  tmp <- as.vector(rbind(dS,dE,dP,dC,dR,dQp,dQc,dH))
  return(list(c(tmp)))
}

# ode function + start condition
epi_ode_seir <- function(C0, C1, ctm_a, ctm_b, Ns = N, ts=seq(1,365,by=1),
                         beta_I, 
                         par_gamma, par_delta, par_epsilon, fp){
  
  # contact matrix
  C.initial <- C0
  C.phase1  <- C1
  
  Nage <- nrow(Ns)
  long_Ns <- as.numeric(t(Ns))
  start <- matrix(0,ncol = Nage, nrow = 8)
  start[1,] <- long_Ns
  
  # 2020-01-07: 3 inital cases
  index_agegroup <- which(age_labels %in% c("[30,35)","[45,50)","[55,60)"))
  start[1,index_agegroup] = start[1,index_agegroup] - 1
  start[4,index_agegroup] = start[4,index_agegroup] + 1

  y <- ode(y=start,t=ts, func=ode_seir, N = Ns, Nage = Nage,
           ctm_a = ctm_a, ctm_b = ctm_b,
           C1st = C.initial, C2st = C.phase1, 
           beta_I = beta_I, par_gamma = par_gamma, par_delta = par_delta, par_epsilon = par_epsilon, fp = fp)
  
  y <- as.data.frame(y)
  
  colnames_dat <- data.frame(id = rep(1:14,8)) %>% 
    arrange(id) %>% # dS,dE,dP,dC,dR,dQp,dQc,dH)
        bind_cols(cell =  rep(c("S","E","P","C","R","Qp","Qc","H"), Nage)) %>% 
    mutate(names = paste0(cell,id))
  
  colnames(y) <- c("time",colnames_dat$names)
  
  cumulative_cases <- y[, str_detect(colnames(y),"H")] # case having symptom
  cumulative_cases_sum <- apply(cumulative_cases,1,sum,na.rm = T) # total case having symptom
  daily_newcases <- tibble(H = cumulative_cases_sum, Hlag = lag(cumulative_cases_sum,1)) %>% 
    mutate(dailynew = H-Hlag)
  
  return(list(cumulative_cases_sum = last(cumulative_cases_sum), daily_newcases = daily_newcases))
}

zj_fit <- function(dat = zhejiang_cases, ranglb = ranglb, x0 = x0, rangub = rangub, test_dur = test_dur, N = N){
  
  ssr_eq <- function(par, data) {
    y_base <- epi_ode_seir(C0 = zhejiang_data$ctm_base*par[6], 
                           C1 = zhejiang_data$ctm_outbreak*par[7], 
                           Ns = N, ts=seq(0,test_dur,by=1),
                           ctm_a = par[1], ctm_b = par[2], 
                           beta_I = exp(par[3]),
                           par_gamma = par[8], par_delta = par[9],
                           par_epsilon = par[4], fp = par[5])
    
    mod_case_ts <- as.numeric(y_base[["daily_newcases"]]$dailynew)
    
    ssrdat <- tibble(date = min(data$date) + 0:(max(data$date) - min(data$date)),
                     mod_ts = mod_case_ts) %>% 
      left_join(data,by = "date") %>% 
      mutate(dailynewcases =  replace(dailynewcases,is.na(dailynewcases),0)) %>% 
      mutate(dif = mod_ts - dailynewcases)
    
    
    ssr = sum(ssrdat$dif^2, na.rm = T)
    return(ssr)
  }
  
  # nloptr optimize
  zhejiang_fit <- nloptr(eval_f = ssr_eq, 
                         x0=x0, lb=ranglb, ub=rangub,
                         opts=list("algorithm" = "NLOPT_LN_SBPLX", 
                                   "local_opts" = list( "algorithm" = "NLOPT_LN_NELDERMEAD", "xtol_rel"  = 1e-5), 
                                   "xtol_rel"= 1e-5,
                                   "print_level" = 3, # 3 
                                   "maxeval" = 5e3),
                         data = dat)
  
  # ***************************************************************
  # best fit scenario
  # ***************************************************************
  solutions <- round(zhejiang_fit$solution,2)
  
  ctm_a_fit = as.numeric(solutions[1])
  ctm_b_fit = as.numeric(solutions[2])
  beta_I_fit  = as.numeric(solutions[3])
  epsilon_fit  = as.numeric(solutions[4])
  fp_fit = as.numeric(solutions[5])
  
  ctm_lo_fit = as.numeric(solutions[6])
  ctm_up_fit = as.numeric(solutions[7])
  
  par_gamma_fit = as.numeric(solutions[8])
  par_delta_fit = as.numeric(solutions[9])

  y_seir_fit <- epi_ode_seir(C0 = zhejiang_data$ctm_base*ctm_lo_fit, 
                             C1 = zhejiang_data$ctm_outbreak*ctm_up_fit, 
                             Ns = N, 
                             # because here is the prediction, so we used the total time of real data 46 days
                             ts=seq(0,46,by=1),
                             ctm_a = ctm_a_fit, ctm_b = ctm_b_fit,
                             beta_I = exp(beta_I_fit),
                             par_gamma = par_gamma_fit, par_delta = par_delta_fit, 
                             par_epsilon = epsilon_fit, fp = fp_fit)
  
  
  
  return(list(nloptr = zhejiang_fit, fit_ts = y_seir_fit))
}