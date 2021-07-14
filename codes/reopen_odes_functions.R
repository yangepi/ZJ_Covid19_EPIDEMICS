# model
age_labels <-  c("[0,5)","[5,10)","[10,15)","[15,20)","[20,25)","[25,30)","[30,35)","[35,40)","[40,45)","[45,50)","[50,55)","[55,60)","[60,65)","65+")

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

ode_seir_reopen <- function(t, y, par, Nage, N, 
                            C.initial, C.outbreak, ctm_a, ctm_b,  
                            beta_I, par_gamma, par_delta, par_epsilon_list, fp_list){
  
  
  # **********************************************************************
  beta_scales <- c(0.4,0.4,0.4,0.4,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8)   # rep(1,14)
  
  if (t<90) { # outbreak peak period 3 months
    changes <- sigmoid(t, a = ctm_a, b = ctm_b) 
    ctm_fun <- C.initial - (C.initial - C.outbreak)*changes
    par_epsilon <- par_epsilon_list[1]
    fp <- fp_list[1]
    C = setup_C(ctm_fun, N, beta_scales)
    
  }
  if (t>=90 & t<180) { # phase I reopen + 3 monthes
    ctm_fun <- C.initial*0.5
    par_epsilon <- par_epsilon_list[2]
    fp <- fp_list[2]
    C = setup_C(ctm_fun, N, beta_scales)
    
  }
  if (t>=180) { # phase II reopen
    ctm_fun <- C.initial*0.8
    par_epsilon <- par_epsilon_list[3]
    fp <- fp_list[3]
    C = setup_C(ctm_fun, N, beta_scales)
    
  }
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


epi_ode_seir_reopen <- function(C0, C1, ctm_a, ctm_b,  Ns = N, ts=seq(1,365,by=1),
                                beta_I, 
                                par_gamma, par_delta, par_epsilon_list, fp_list){
  
  Nage <- nrow(Ns)
  long_Ns <- as.numeric(t(Ns))
  start <- matrix(0,ncol = Nage, nrow = 8)
  start[1,] <- long_Ns
  
  # 2020-01-07: 3 inital cases
  index_agegroup <- which(age_labels %in% c("[30,35)","[45,50)","[55,60)"))
  start[1,index_agegroup] = start[1,index_agegroup] - 1
  start[4,index_agegroup] = start[4,index_agegroup] + 1
  # browser()
  
  y <- ode(y=start,t=ts, func=ode_seir_reopen, N = Ns, Nage = Nage,
           ctm_a = ctm_a, ctm_b = ctm_b, 
           C.initial = C0, C.outbreak = C1, 
           beta_I = beta_I, par_gamma = par_gamma, par_delta = par_delta, par_epsilon_list = par_epsilon_list, fp_list = fp_list)
  
  y <- as.data.frame(y)
  
  colnames_dat <- data.frame(id = rep(1:14,8)) %>% 
    arrange(id) %>% # dS,dE,dP,dC,dR,dQp,dQc,dH)
    bind_cols(cell =  rep(c("S","E","P","C","R","Qp","Qc","H"), Nage)) %>% 
    mutate(names = paste0(cell,id))
  
  colnames(y) <- c("time",colnames_dat$names)
  
  cumulative_cases <- y[, str_detect(colnames(y),"H")] # case having symptom
  cumulative_cases_sum <- apply(cumulative_cases,1,sum,na.rm = T) # case having symptom
  daily_newcases <- tibble(H = cumulative_cases_sum, Hlag = lag(cumulative_cases_sum,1)) %>% 
    mutate(dailynew = H-Hlag)
  
  return(list(cumulative_cases_sum = last(cumulative_cases_sum), daily_newcases = daily_newcases))
}