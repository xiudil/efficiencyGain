#################################################################
### estimating efficiency gain of covariate adjusted analysis ###
#################################################################
library(MASS)
library(dplyr)
library(nleqslv)
library(mvtnorm)

source("helper_fns.R")

page = c(0.01,0.09,0.12,0.13,0.18,0.22,0.25)
pevent = matrix(c(0,0.01,0.03,0.08,0.11,0.17,0.37,0,0.18,0.32,0.31,0.37,0.47,0.35,
                  rep(0,7)),ncol=3,byrow=FALSE)
pevent[,3] = 1-pevent[,1]-pevent[,2]

p.cond = function(w1,w2){
  # p1.cond = expit(-x)*expit(1.5*x+w2^2)
  # p2.cond = (1-expit(-x))*expit(1.5*x+w2^2)
  p1.cond = expit(1.5*w1-2*w2^2)
  p2.cond = expit(1.5*w1+w2^2)-expit(1.5*x-2*w2^2)
  p3.cond = 1-expit(1.5*w1+w2^2)
  return(c(p1.cond, p2.cond, p3.cond))
}





#############################
### this is for bootstrap ###
#############################

### draw a second layer resample and compute treatment effect estimate on second layer resample
# Yboot is first layer outcome, vector
# Wboot is first layer covariate, data.frame
# n: first layer sample size
# oversample: n^oversample is second layer sample size, taken to be 1.3 in the sims
# pi: treatment probability
gen_theta_boot = function(Yboot,Wboot,n,oversample,pi,u=I){
  resample2 = sample(x=n,size = n^oversample,replace=TRUE)
  Yboot2 = Yboot[resample2]
  Wboot2 = Wboot[resample2,]
  trt = rbinom(n^oversample,size=1,prob = pi)
  
  
  DIM_u = mean(u(Yboot2[trt == 1])) - mean(u(Yboot2[trt == 0]))
  modtrt = ord_work_param_clean(Yboot2[trt==1], Wboot2[trt==1,], Wboot2)
  pevent_mod_trt = modtrt[[1]];cdf_mod_trt = modtrt[[2]]; rm(modtrt)
  modctr = ord_work_param_clean(Yboot2[trt==0], Wboot2[trt==0,], Wboot2)
  pevent_mod_ctr = modctr[[1]];cdf_mod_ctr = modctr[[2]]; rm(modctr)
  condmean_trt = u(c(1:3)) %*% t(pevent_mod_trt)
  condmean_ctr = u(c(1:3)) %*% t(pevent_mod_ctr)
  DIM_a = mean(condmean_trt) - mean(condmean_ctr)
  
  marginal_trt = rep(0,3);marginal_ctr = marginal_trt;h = marginal_trt
  for (i in 1:3){
    marginal_trt[i] = sum(Yboot2==i & trt ==1)/sum(trt==1)
    marginal_ctr[i] = sum(Yboot2==i & trt ==0)/sum(trt==0)
  }
  for (i in 1:3){
    h[i] = sum(marginal_trt[i:3])-0.5*marginal_trt[i]
  }
  MW_u = sum(h*marginal_ctr)
  
  trt_mod = apply(pevent_mod_trt,2,mean)
  ctr_mod = apply(pevent_mod_ctr,2,mean)
  h_mod = rep(0,3)
  for (i in 1:3){
    h_mod[i] = sum(trt_mod[i:3])-0.5*trt_mod[i]
  }
  MW_a = sum(h_mod*ctr_mod)
  
  LOR_u = 0
  LOR_a = 0
  for (i in 1:2){
    LOR_u = LOR_u + logit(sum(marginal_trt[1:i])) - logit(sum(marginal_ctr[1:i]))
    LOR_a = LOR_a + logit(sum(trt_mod[1:i])) - logit(sum(ctr_mod[1:i]))
  }
  LOR_u = LOR_u/2
  LOR_a = LOR_a/2
  
  result = c(DIM_u,DIM_a,MW_u,MW_a,LOR_u,LOR_a)
  names(result) = c("DIM_u","DIM_a","MW_u","MW_a","LOR_u","LOR_a")
  return(result)
}

# link specifies the scale on which the RE and standard error are reported, can be identity, log or logit
# Bone: number of first layer resample, 100 in the simulations reported
# Btwo: number of second layer resample, 2000 in the simulations reported
estimate_gain_boot = function(y,w,oversample,pi,link,Bone,Btwo){
  thetas_est = replicate(Btwo,gen_theta_boot(y,w,n,oversample,pi))
  vars_est = apply(thetas_est,1,var)*n^oversample
  releff_DIM = vars_est[2]/vars_est[1]
  releff_MW = vars_est[4]/vars_est[3]
  releff_LOR = vars_est[6]/vars_est[5]
  est = c(releff_DIM,releff_MW,releff_LOR)
  
  releff = matrix(0,3,Bone)
  for (i in 1:Bone){
    resample = sample(x=n,size=n,replace = TRUE)
    yboot = y[resample]
    wboot = w[resample,]
    thetas = replicate(Btwo,gen_theta_boot(yboot,wboot,n,oversample,pi))
    vars = apply(thetas,1,var)*n^oversample
    if (link == "logit"){
      releff[1,i] = logit(vars[2]/vars[1]);releff[2,i] = logit(vars[4]/vars[3]);releff[3,i] = logit(vars[6]/vars[5])
    } else if (link == "log"){
      releff[1,i] = log(vars[2]/vars[1]);releff[2,i] = log(vars[4]/vars[3]);releff[3,i] = log(vars[6]/vars[5])
    } else if (link == "identity"){
      releff[1,i] = vars[2]/vars[1];releff[2,i] = vars[4]/vars[3];releff[3,i] = vars[6]/vars[5]
    }
  }
  sds = apply(releff,1,sd,na.rm=T)
  return(c(est,sds))
}

### example
n = 500
x = rmvnorm(n, mean=c(0,0), sigma = matrix(c(1,0.3,0.3,1),ncol=2))
x2 = runif(n,min=-1,max=1)
x1 = x[,1]
x3 = x[,2]
y = rep(0,n)
for (i in 1:n){
  y[i] = sample(x=3, size=1, prob=p.cond(x1[i],x2[i]))
}

w = data.frame("X1" = x1, "X2"=x2, "X3"=x3)

estimate_gain_boot(y=y,w=w, oversample=1.3,pi=0.5,link="logit",Bone=100,Btwo=2000)

###################
### simulations ###
###################

do.one.boot = function(n,oversample,pi,link,Bone,Btwo,setting){
  if (setting == 1){
    dat = gen_ordinal_data(n,page,pevent)
    y = dat$Y
    w = data.frame("W" = dat$W)
  } else if (setting == 2){
    x = rmvnorm(n, mean=c(0,0), sigma = matrix(c(1,0.3,0.3,1),ncol=2))
    x2 = runif(n,min=-1,max=1)
    x1 = x[,1]
    x3 = x[,2]
    y = rep(0,n)
    for (i in 1:n){
      y[i] = sample(x=3, size=1, prob=p.cond(x1[i],x2[i]))
    }
    w = data.frame("X1" = x1, "X2"=x2, "X3"=x3)
  }
  
  result = estimate_gain_boot(y,w,oversample,pi,link,Bone,Btwo)
  return(result)
}

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
set.seed(68*seed)
filename = paste("n",nn,"seed",seed,"largeB.RData",sep="")
system.time({rr = replicate(2,do.one.boot(nn,oversample=1.3,pi=0.5,link="logit",Bone=100,Btwo=2000,setting=1))})
save(rr,file=filename)
