###############################################################################################
### estimating efficiency gain of covariate adjusted analysis: analytical, survival outcomes ##
###############################################################################################
library(MASS)
library(survival)
library(survSuperLearner)

### library of learning algorithms in super learner
#   can include additional learners if desired
event.SL.library <- cens.SL.library <- c("survSL.km", "survSL.coxph", "survSL.gam")#, "survSL.rfsrc")

##############################################
### survival outcomes: estimating gain + CI ##
##############################################

### estimate nuisance functions --- S: conditional survival function, h: conditional hazard function, C: conditional distribution of censoring time
### dat is a data.frame with the following columns: time, censoring status (event, 1=event, 0=censoring), covariate W, 1-censoring status
### tt: a sequence of time points on which we estimate the nuisance function
### method: method to estimate the nuisance functions
###         coxph = sequence of coxph models + BIC
###         super = super learning without cross-fitting (used for our simulations)
###         supercf = super learning with cross-fitting
### output matrix of dimension n by length(tt)
calc_shc = function(dat,tt,method){
  k = length(tt)
  n = dim(dat)[1]
  S = matrix(0,n,k)
  h = S
  C = S
  H = S
  
  # using a sequence of coxph model with polynomials of W with increasing order, 1,3,5,7
  # choose the model with best BIC
  if (method == "coxph"){
    mod_lm = coxph(Surv(time,event)~W,data = dat)
    mod_poly3 = coxph(Surv(time,event)~W + I(W^2) + I(W^3),data = dat)
    mod_poly5 = coxph(Surv(time,event)~W + I(W^2) + I(W^3) + I(W^4) + I(W^5),data = dat)
    mod_poly7 = coxph(Surv(time,event)~W + I(W^2) + I(W^3) + I(W^4) + I(W^5)+ I(W^6) + I(W^7),data = dat)
    mods = list(mod_lm,mod_poly3,mod_poly5,mod_poly7)
    index = which.min(BIC(mod_lm,mod_poly3,mod_poly5,mod_poly7)$BIC)
    
    mod_lm_c = coxph(Surv(time,censor)~W,data = dat)
    mod_poly3_c = coxph(Surv(time,censor)~W + I(W^2) + I(W^3),data = dat)
    mod_poly5_c = coxph(Surv(time,censor)~W + I(W^2) + I(W^3) + I(W^4) + I(W^5),data = dat)
    mod_poly7_c = coxph(Surv(time,censor)~W + I(W^2) + I(W^3) + I(W^4) + I(W^5)+ I(W^6) + I(W^7),data = dat)
    mods_c = list(mod_lm_c,mod_poly3_c,mod_poly5_c,mod_poly7_c)
    index_c = which.min(BIC(mod_lm_c,mod_poly3_c,mod_poly5_c,mod_poly7_c)$BIC)
    
    for (i in 1:k){
      newdata = data.frame(cbind(rep(tt[i],n),rep(1,n),dat$W,rep(1,n)))
      colnames(newdata) = c("time","event","W","censor")
      pred = predict(mods[[index]],type="expected",newdata = newdata)
      S[,i] = exp(-pred)
      H[,i] = pred
      pred_c = predict(mods_c[[index_c]],type="expected",newdata = newdata)
      C[,i] = exp(-pred_c)
    }
    C = cbind(rep(1,n),C[,-k])
    h = H - cbind(rep(0,n),H[,-k])
  } else if (method == "super"){
    # estimate nuisance with survSuperLearner
    mod_super = survSuperLearner(time = dat$time, event = dat$event, X = data.frame("W"=dat$W), newX = data.frame("W"=dat$W), 
                                 new.times = tt, event.SL.library = event.SL.library, cens.SL.library = cens.SL.library, verbose = FALSE)
    C = mod_super$cens.SL.predict
    S = mod_super$event.SL.predict
    h = 1 - S / cbind(rep(1,n),S[,-k])
  } else if (method == "supercf"){
    # estimate nuisance with survSuperLearner using cross-fitting
    nfold = 10
    fold.size = ceiling(n/nfold)
    fold.id = NULL
    for (i in 1:(nfold-1)){
      fold.id = c(fold.id,rep(i,fold.size))
    }
    fold.id = c(fold.id,rep(nfold,n-(nfold-1)*fold.size))
    fold.id = fold.id[sample(c(1:length(fold.id)),length(fold.id))]
    
    for (i in 1:nfold){
      outfoldx = data.frame("W" = dat$W[fold.id!=i])
      infoldx  = data.frame("W" = dat$W[fold.id==i])
      mod_super = survSuperLearner(time = dat$time[fold.id!=i], event = dat$event[fold.id!=i], X = outfoldx, newX = infoldx, 
                                   new.times = tt, event.SL.library = event.SL.library, cens.SL.library = cens.SL.library, verbose = FALSE)
      C[fold.id==i,] = mod_super$cens.SL.predict
      S[fold.id==i,] = mod_super$event.SL.predict
      h[fold.id==i,] = 1 - S[fold.id==i] / cbind(rep(1,n),S[fold.id==i,-k])
    }
  }
  return(list(S,h,C))
}

########################################################
### estimate relative efficiency phi_f for RD and RR
########################################################

# margg: marginal distribution of C, censoring time, in the future trial, user-specified
# g: conditional distribution of C|W in the future trial, user-specified
# dat is a data.frame with the following columns: time, censoring status (event, 1=event, 0=censoring), covariate W, 1-censoring status
# tt: all unique event times before some time of interest t0, the sequence of time points on which we estimate the nuisance function
# S, h and C are the nuisance function estimates from the previous function
est_gain_RD = function(dat,tt,S,h,C,margg,g){
  k = length(tt)
  marginalS = rep(0,k)
  n = dim(dat)[1]
  
  ### unadjusted variance
  tau.mat =matrix(0,n,k)
  tau.mat[,1] = (-1)/(S[,1]*C[,1])*((dat$time == tt[1] & dat$event == 1) - h[,1]*(dat$time >= tt[1]))
  for (u in 2:k){
    tau.mat[,u] = (-1)/(S[,u]*C[,u])*((dat$time == tt[u] & dat$event == 1) - h[,u]*(dat$time >= tt[u])) + tau.mat[,(u-1)]
  }
  for (u in 1:k){
    tau.mat[,u] = tau.mat[,u]*S[,u]
  }
  marginalS = apply(tau.mat + S,2,mean)
  phit = sweep(tau.mat + S,2,marginalS)
  sigma_u = sum((1/marginalS - 1/c(1,marginalS[-k]))/margg)*marginalS[k]^2
  firstpiece = 2*marginalS[k]*phit[,k]*(sum(1/(marginalS*margg))-sum(1/(c(1,marginalS[-k])*margg)))
  secondpiece = 1/margg[1]*(-phit[,1]/marginalS[1]^2)
  for (i in 2:k){
    secondpiece = secondpiece + 1/margg[i]*(-phit[,i]/marginalS[i]^2 + phit[,(i-1)]/marginalS[(i-1)]^2)
  } 
  secondpiece = marginalS[k]^2*secondpiece
  var_u = mean((firstpiece + secondpiece)^2)
  
  ### adjusted variance
  EE = (2*S[,k]*tau.mat[,k]*(1/S[,1]-1) - S[,k]^2/S[,1]^2*tau.mat[,1] + S[,k]^2*(1/S[,1]-1))/g[,1]
  for (u in 2:k){
    EE = EE + (2*S[,k]*tau.mat[,k]*(1/S[,u]-1/S[,u-1]) - S[,k]^2/S[,u]^2*tau.mat[,u] + S[,k]^2/S[,u-1]^2*tau.mat[,u-1] + S[,k]^2*(1/S[,u]-1/S[,u-1]))/g[,u]
  }
  
  sigma_a = mean(EE)
  var_a = mean((EE-sigma_a)^2)
  
  releff = sigma_a/sigma_u
  var_releff = mean(((EE-sigma_a)/sigma_u - sigma_a/sigma_u^2*(firstpiece + secondpiece))^2)
  
  return(c(sigma_u,var_u,sigma_a,var_a,releff,var_releff))
}

########################################################
### estimate relative efficiency phi_f for RMST
########################################################

# t0 is the time of interest
est_gain_RMST = function(dat,tt,S,h,C,margg,g,t0){
  k = length(tt)
  n = dim(dat)[1]
  
  # unadjusted variance
  marginalS = rep(0,k)
  tau.mat =matrix(0,n,k)
  tau.mat[,1] = (-1)/(S[,1]*C[,1])*((dat$time == tt[1] & dat$event == 1) - h[,1]*(dat$time >= tt[1]))
  for (u in 2:k){
    tau.mat[,u] = (-1)/(S[,u]*C[,u])*((dat$time == tt[u] & dat$event == 1) - h[,u]*(dat$time >= tt[u])) + tau.mat[,(u-1)]
  }
  for (u in 1:k){
    tau.mat[,u] = tau.mat[,u]*S[,u]
  }
  marg.S.adjusted = apply(tau.mat + S,2,mean)
  phit = sweep(tau.mat + S,2,marg.S.adjusted)
  phit.sum = matrix(0,n,k)
  phit.sum[,k] = phit[,k]*(t0 - tt[k])
  marg.S.sum = rep(0,k)
  marg.S.sum[k] = marg.S.adjusted[k]*(t0 - tt[k])
  for (u in 1:(k-1)){
    phit.sum[,k-u] = phit[,k-u]*(tt[k-u+1] - tt[k-u]) + phit.sum[,k-u+1]
    marg.S.sum[k-u] = marg.S.adjusted[k-u]*(tt[k-u+1] - tt[k-u]) + marg.S.sum[k-u+1]
  }
  
  sigma_u = sum(marg.S.sum^2*(1/marg.S.adjusted - 1/c(1,marg.S.adjusted[-k]))/margg)
  IF_u = 2*(1/marg.S.adjusted[1] - 1)/margg[1]*marg.S.sum[1]*phit.sum[,1] - marg.S.sum[1]^2/marg.S.adjusted[1]^2*phit[,1]
  for (u in 2:k){
    IF_u = IF_u + 2*(1/marg.S.adjusted[u] - 1/marg.S.adjusted[u-1])/margg[u]*marg.S.sum[u]*phit.sum[,u] - 
      marg.S.sum[u]^2/marg.S.adjusted[u]^2*phit[,u]/margg[u] + 
      marg.S.sum[u]^2/marg.S.adjusted[u-1]^2*phit[,u-1]/margg[u]
  }
  var_u = mean(IF_u^2)
  
  # adjusted variance
  EE = 0
  tau.mat.sum = matrix(0,n,k)
  S.sum = matrix(0,n,k)
  tau.mat.sum[,k] = tau.mat[,k]*(t0 - tt[k])
  S.sum[,k] = S[,k]*(t0 - tt[k])
  for (u in 1:(k-1)){
    tau.mat.sum[,(k-u)] = tau.mat[,(k-u)]*(tt[k-u+1] - tt[k-u])  + tau.mat.sum[,(k-u+1)]
    S.sum[,(k-u)] = S[,(k-u)]*(tt[k-u+1] - tt[k-u]) + S.sum[,(k-u+1)]
  }
  EE = EE + S.sum[,1]^2*(1/S[,1] - 1)/g[,1] +2*(1/S[,1] - 1)/g[,1]*tau.mat.sum[,1]*S.sum[,1]-S.sum[,1]^2*tau.mat[,1]/(S[,1]^2)/g[,1]
  for (u in 2:k){
    term.1 = 2*(1/S[,u] - 1/S[,(u-1)])*tau.mat.sum[,u]*S.sum[,u]/g[,u]
    term.2 = S.sum[,u]^2*tau.mat[,u]/(S[,u]^2)/g[,u]
    term.3 = S.sum[,u]^2*tau.mat[,u-1]/(S[,u-1]^2)/g[,u]
    term.4 = S.sum[,u]^2*(1/S[,u] - 1/S[,u-1])/g[,u]
    EE = EE + term.1 - term.2 + term.3 + term.4
  }
  
  sigma_a = mean(EE)
  var_a = mean((EE-sigma_a)^2)
  
  releff = sigma_a/sigma_u
  var_releff = mean(((EE-sigma_a)/sigma_u - sigma_a/sigma_u^2*IF_u)^2)
  
  return(c(sigma_u,var_u,sigma_a,var_a,releff,var_releff))
}

### example
n=1000
t0 = 3 #time point of interest
W = runif(n)
tevent = rexp(n,rate = 0.1*(1+9*W))
tcensor = rexp(n,rate = 0.1)
time = pmin(tevent,tcensor)
event = (time == tevent)
tt = sort(unique(time[event ==1 & time<=t0]))

dat = data.frame(cbind(time,event,W,1-event))
colnames(dat) = c("time","event","W","censor")
g = matrix(rep(1-pexp(tt,rate=0.1),n),byrow=TRUE,ncol=length(tt))
margg = g[1,]
estimates = calc_shc(dat,tt,method = "super")
S = estimates[[1]];h = estimates[[2]]; C = estimates[[3]];rm(estimates)
result_super = est_gain_RD(dat,tt,S,h,C,margg,g)

### CI: log scale then transform
CIlower = exp(log(result_super[5]) - qnorm(0.975)*sqrt(result_super[6]/n)/result_super[5])
CIupper = exp(log(result_super[5]) + qnorm(0.975)*sqrt(result_super[6]/n)/result_super[5])

CIlower = exp(log(result_super[11]) -qnorm(0.975)*sqrt(result_super[12]/n)/(result_super[11]))
CIupper = exp(log(result_super[11]) +qnorm(0.975)*sqrt(result_super[12]/n)/(result_super[11]))

### CI: original scale
CIlower = result_super[5] - qnorm(0.975)*sqrt(result_super[6]/n)
CIupper = result_super[5] + qnorm(0.975)*sqrt(result_super[6]/n)

##########################
### simulation studies ###
##########################
# setting = 1, discrete time
#           2, discrete time, dependent censoring
#           3, continuous time

# one simulation
# n = sample size
# t0 = time of interest

do.one = function(n,t0, setting){
  # generate data
  if (setting == 1){
    W = runif(n)
    tevent = rexp(n,rate = 0.1*(1+9*W))
    tcensor = rexp(n,rate = 0.1)
    tevent = floor(tevent*5)+1
    tcensor = floor(tcensor*5)+1
    tcensor[tcensor>50] = 50
    
    time = pmin(tevent,tcensor)
    event = (time == tevent)
    tt = sort(unique(time[event ==1 & time<=t0]))
  } else if (setting == 2){
    W = runif(n)
    tevent = rexp(n,rate = 0.1*(1+9*W))
    tcensor = rexp(n,rate = 1-0.8*W)
    tevent = floor(tevent*5)+1
    tcensor = floor(tcensor*5)+1
    tcensor[tcensor>50] = 50
    
    time = pmin(tevent,tcensor)
    event = (time == tevent)
    tt = sort(unique(time[event ==1 & time<=t0]))
  } else if (setting ==3){
    W = runif(n)
    tevent = rexp(n,rate = 0.1*(1+9*W))
    tcensor = rexp(n,rate = 0.1)
    time = pmin(tevent,tcensor)
    event = (time == tevent)
    tt = sort(unique(time[event ==1 & time<=t0]))
  }
  
  
  dat = data.frame(cbind(time,event,W,1-event))
  colnames(dat) = c("time","event","W","censor")
  if (setting == 1 | setting ==2){
    g = matrix(rep(1-pexp((tt-1)*0.2,rate=0.1),n),byrow=TRUE,ncol=length(tt))
  } else if (setting == 3){
    g = matrix(rep(1-pexp(tt,rate=0.1),n),byrow=TRUE,ncol=length(tt))
  }
  margg = g[1,]
  
  estimates = calc_shc(dat,tt,method = "super")
  S = estimates[[1]];h = estimates[[2]]; C = estimates[[3]];rm(estimates)
  result_super = est_gain_RD(dat,tt,S,h,C,margg,g)
  #result_super = est_gain_RMST(dat,tt,S,h,C,margg,g)
  
  estimates = calc_shc(dat,tt,method = "supercf")
  S = estimates[[1]];h = estimates[[2]]; C = estimates[[3]];rm(estimates)
  result_supercf = est_gain_RD(dat,tt,S,h,C,margg,g)
  #result_super = est_gain_RMST(dat,tt,S,h,C,margg,g)
  
  return(c(result_super,result_supercf))
}

# many simulations, run on cluster, each job run one simulation replication
args <- commandArgs(TRUE)
print(args)
if(length(args) == 0) {
  print("No arguments supplied.")
  seed = 1
} else {
  for (i in 1:length(args)) {
    eval(parse(text = args[i]))
  }
}
seed = S
nn = N
t0 = Q

### for discrete time, set t0 to either 6, 11, or 16
### for continuous time, set t0 to 1, 2 or 3
set.seed(9*seed)
filename = paste("seed",seed,"n",nn,"t",t0,".RData",sep="")
suppressWarnings({results = replicate(1,do.one(nn,t0,setting=3))})
save(results,file=filename)


### functions to compute true RE ###
# discrete time, RD
margsurv = function(a,b,t){
  ss = (exp(-a*t)-exp(-a*(1+b)*t))/(a*b*t)
  return(ss)
}
dis_time = seq(1,50,1)
margsurv_dis = margsurv(0.1,9,dis_time*0.2)

calc_truth_RD = function(margsurv_dis,t0){
  g = 1-pexp((seq(1,t0,1)-1)*0.2,rate=0.1)
  var_u = margsurv_dis[t0]^2*sum((1/margsurv_dis[1:t0]- 1/c(1,margsurv_dis[1:t0-1]))/g[1:t0])
  return(var_u)
}
integrand_RD = function(w,a0,b0,k){
  value = 0
  for (mm in 1:k){
    value = value + exp(-a0*(1+b0*w)*0.2*(k+k))*(exp(a0*(1+b0*w)*0.2*mm)-exp(a0*(1+b0*w)*0.2*(mm-1)))*exp(0.1*0.2*(mm-1))
  }
  return(value)
}

sigma_u_rd = calc_truth_RD(margsurv_dis,tt)
sigma_a_rd = integrate(f=integrand_RD,lower=0,upper=1,a0=0.1,b0=9,k=tt)$value
re = sigma_a_rd / sigma_u_rd

# discrete time, RMST
calc_truth_RMST = function(margsurv_dis,t0){
  g = 1-pexp((seq(1,t0,1)-1)*0.2,rate=0.1)
  var_u = 0
  for (j in 1:t0){
    for (l in 1:t0){
      newterm = sum((1/margsurv_dis[1:min(j,l)]-1/c(1,margsurv_dis[1:(min(j,l)-1)]))/g[1:min(j,l)])
      var_u = var_u + margsurv_dis[j]*margsurv_dis[l]*newterm
    }
  }
  return(var_u)
}
integrand_RMST = function(w,a0,b0,k){
  value = 0
  for (jj in 1:k){
    for (ll in 1:k){
      for (mm in 1:min(jj,ll)){
        value = value + exp(-a0*(1+b0*w)*0.2*(jj+ll))*(exp(a0*(1+b0*w)*0.2*mm)-exp(a0*(1+b0*w)*0.2*(mm-1)))*exp(0.1*0.2*(mm-1))
      }
    }
  }
  return(value)
}
re = integrate(f=integrand_RMST,lower=0,upper=1,a0=0.1,b0=9,k=15)$value/calc_truth_RMST(margsurv_dis,15)

# continuous time, RD
integrand = function(t,a,b){
  value = ((a*exp(-a*t)-a*(1+b)*exp(-(a+a*b)*t))*a*b*t+a*b*(exp(-a*t)-exp(-a*(1+b)*t)))/(exp(-a*t)-exp(-a*(1+b)*t))^2*exp(0.1*t)
  return(value)
}
integrand_a = function(t00,a0,b0,w){
  value = exp(-2*(a0*(1+b0*w)*t00))*(a0+a0*b0*w)/(a0+a0*b0*w+0.1)*(exp((a0+a0*b0*w+0.1)*t00)-1)
  return(value)
}

calc_truth_analytical = function(t0,a,b){
  var_u = ((exp(-a*t0)-exp(-a*(1+b)*t0))/(a*b*t0))^2*integrate(f=integrand,lower=0,upper=t0,a=a,b=b)$value
  var_a = integrate(f=integrand_a,lower=0,upper=1,t00=t0,a0=a,b0=b)$value
  return(c(var_u,var_a))
}
truth = calc_truth_analytical(3,0.1,9)
re = truth[2]/truth[1]

# continuous time, RMST
integrand.rmst.a = function(x,a,b){
  w=x[1]
  t=x[2]
  val = (a*(1+b*w))^(-1)*exp(a*(1+b*w)*t + 0.1*t)*(exp(-a*(1+b*w)*t)-exp(-a*(1+b*w)*3))^2
  return(val)
}
sigma_a_true = hcubature(f=integrand.rmst.a,lowerLimit=c(0,0),upperLimit=c(1,3),a=0.1,b=9)$integral

inner.integrand.rmst = function(t,a,b){
  inner.integral = integrate(function(y){exp(-(a*y))*(1-exp(-(a*b*y)))/(a*b*y)},lower=t,upper=3)$value
  return(inner.integral)
}
integrand.rmst = function(t,a,b){
  outer.integrand = ((a*exp(-a*t)-a*(1+b)*exp(-(a+a*b)*t))*a*b*t+a*b*(exp(-a*t)-exp(-a*(1+b)*t)))/(exp(-a*t)-exp(-a*(1+b)*t))^2*exp(0.1*t)
  value = inner.integrand.rmst(t,a,b)^2*outer.integrand
  return(value)
}
t.grid = seq(0.00001,3,by=0.00001)
t.grid.value = sapply(t.grid,integrand.rmst,a=0.1,b=9)
sigma_u_true = sum(t.grid.value)*0.00001

re.rmst.cont = sigma_a_true/sigma_u_true