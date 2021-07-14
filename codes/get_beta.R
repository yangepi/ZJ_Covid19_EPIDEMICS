get_beta <- function(C, propns, inf_period=5, r0){
  beta_scales <- c(0.4,0.4,0.4,0.4,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8)
  fi <- kron(propns, ones(1, length(propns)))
  fj <- kron(t(propns), ones(length(propns),1))
  M <- C * beta_scales * fi/fj
  Rx <- max(Re(eigen(M)$values))
  beta <- r0/(Rx*inf_period)
  beta
}
