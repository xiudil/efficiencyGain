####################
# COVID-19 Data    #
# 05/07/2022       #
# compare with existing methods # 
####################

# computationally intensive, so we can run cluster
# P is the argument specifying seed
# need to fit adjusted and unadjusted estimators B times on simulated data

args <- commandArgs(TRUE)
if(length(args) == 0) {
  print("No arguments supplied.")
} else {
  for (i in 1:length(args)) {
    eval(parse(text = args[i]))
  }
}
args = P

library(MASS)
library(SuperLearner)
library(tidyverse)
library(nleqslv)
library(survival)
library(glmnet)
library(MASS)
library(survival)
library(survSuperLearner)


data = read.csv("covid_analysis/mock_data.csv")

logit = function(x){
  return(log(x/(1-x)))
}
expit = function(x){
  return(exp(x)/(1+exp(x)))
}

### --------------------ordinal outcome--------------------
# fit adjusted estimator based on working models
my_polr_clean = function(y,w){
  ycat = length(unique(y))
  n = dim(w)[1]
  binary.y = NULL
  x = matrix(0,(ycat-1)*n,ycat-1)
  ww = NULL
  for (i in 1:(ycat-1)){
    binary.y = c(binary.y,1*(y<=i))
    x[(n*(i-1)+1):(n*i),i] = 1
    ww = rbind(ww,as.matrix(w))
  }
  fit = glm(binary.y~cbind(x,ww)-1,family="binomial")
  
  coef = fit$coefficients
  return(coef)
}
ord_work_param_clean = function(y,w,neww=w){
  coef = my_polr_clean(y,w)
  ycat = length(unique(y))
  
  pevent_mod = matrix(0,dim(neww)[1],ycat)
  for (i in 1:ycat){
    if (i == 1){
      pevent_mod[,i] = expit(coef[i]+coef[ycat:length(coef)] %*% t(as.matrix(neww)))
    } else if (i < ycat & i>1){
      pevent_mod[,i] = expit(coef[i]+coef[ycat:length(coef)] %*% t(as.matrix(neww))) - pevent_mod[,i-1]
    } else {
      pevent_mod[,i] = 1 - rowSums(pevent_mod[,1:i-1])
    }
  }
  
  cdf_mod = pevent_mod[,1:(ycat-1)]
  for (i in 1:(ycat-1)){
    if (i == 1){
      cdf_mod[,i] = pevent_mod[,i]
    } else{
      cdf_mod[,i] = rowSums(pevent_mod[,1:i])
    }
  }
  return(list(pevent_mod,cdf_mod))
}


#Lib = c("SL.mean","SL.glm","SL.glm.interaction","SL.gam","SL.ranger")
learners = create.Learner("SL.gam",tune=list(deg.gam = c(5,10,20)))
Lib = c("SL.gam","SL.glm","SL.glm.interaction",learners$names,"SL.ranger")

### this function estimate the treatment effect (difference in mean)
###   on a dataset simulated from the original data
###   using AIPW estimator and working-model-based estimators
###   adjusting for covariates in which.W
###   as well as using unadjusted estimator
#          nuisance functions in AIPW estimated via super learner with library = Lib
one <- function(which.W,P){
  set.seed(P)
  # simulate data by resampling
  data = data[sample(1:nrow(data),nrow(data),replace = T),]
  
  # generate treatment variables
  A <- rbinom(nrow(data),1,0.5)
  W <- data[,which.W,drop  = FALSE]
  Y <- data$y_ord
  
  sl  <- SuperLearner(Y = Y, X = data.frame(A,W),family = gaussian,SL.library = Lib)
  cond1 <- predict(sl,newdata = data.frame(A=1,W))$pred
  cond0 <- predict(sl,newdata = data.frame(A=0,W))$pred
  
  unest <- mean(Y[A==1]) - mean(Y[A==0])
  est <- mean((A==1)/0.5*(Y - cond1) + cond1) - mean((A==0)/0.5*(Y - cond0) + cond0)
  ### working model
  modtrt = ord_work_param_clean(Y[A==1],W[A==1,,drop=FALSE],W)
  pevent_mod_trt = modtrt[[1]];cdf_mod_trt = modtrt[[2]]; rm(modtrt)
  modctr = ord_work_param_clean(Y[A==0],W[A==0,,drop=FALSE],W)
  pevent_mod_ctr = modctr[[1]];cdf_mod_ctr = modctr[[2]]; rm(modctr)
  condmean_trt = c(1:3) %*% t(pevent_mod_trt)
  condmean_ctr = c(1:3) %*% t(pevent_mod_ctr)
  DIM_m = mean(condmean_trt - condmean_ctr)
  
  return(list(adj = est, unadj = unest,model = DIM_m))
  
}
allP <-  c("Gender","race","cardiovascular_disease","hypertension","chronic_kidney_disease", "cholesterol_meds","hypertension_meds","dm","age")
individual = sapply(allP, one,P=P)
all = one(allP,P=P)
save(individual,all,file = paste0("./",P,"_resultsRD.rda"))


#### survival outcome
event.SL.library <- cens.SL.library <- c("survSL.km", "survSL.coxph", "survSL.gam","survSL.expreg","survSL.weibreg","survSL.loglogreg")#, "survSL.rfsrc")

### this function estimate the treatment effect (rd & rmst)
###   on a dataset simulated from the original data
###   using AIPCW estimator adjusting for covariates (implemented in the function) in which.W
###   as well as using unadjusted estimator
#          nuisance functions in AIPCW estimated via super learner with library in event.SL.library
one.surv = function(which.W,P){
  set.seed(P)
  data = data[sample(1:nrow(data),nrow(data),replace = T),]
  n = nrow(data)
  t0 = 21000
  A = rbinom(n,1,0.5)
  W = data[,which.W,drop  = FALSE]
  
  tt = c(sort(unique(data$minutes_to_discharge_or_censor[data$event == 1 & data$minutes_to_discharge_or_censor < t0])),t0)
  k = length(tt)
  mod_super = survSuperLearner(time = data$minutes_to_discharge_or_censor, event = data$event, X = data.frame(A,W), newX=rbind(data.frame(A=1,W),data.frame(A=0,W)),
                               new.times = tt, event.SL.library = event.SL.library, cens.SL.library = cens.SL.library, verbose = FALSE)
  C1 = mod_super$cens.SL.predict[1:n,]
  S1 = mod_super$event.SL.predict[1:n,]
  h1 = 1 - S1 / cbind(rep(1,n),S1[,-k])
  
  C0 = mod_super$cens.SL.predict[-c(1:n),]
  S0 = mod_super$event.SL.predict[-c(1:n),]
  h0 = 1 - S0 / cbind(rep(1,n),S0[,-k])
  
  tau.mat1 =matrix(0,n,k)
  tau.mat1[,1] = (-1)/(S1[,1]*C1[,1])*((data$minutes_to_discharge_or_censor == tt[1] & data$event == 1) - h1[,1]*(data$minutes_to_discharge_or_censor >= tt[1]))
  for (u in 2:k){
    tau.mat1[,u] = (-1)/(S1[,u]*C1[,u])*((data$minutes_to_discharge_or_censor == tt[u] & data$event == 1) - h1[,u]*(data$minutes_to_discharge_or_censor >= tt[u])) + tau.mat1[,(u-1)]
  }
  for (u in 1:k){
    tau.mat1[,u] = tau.mat1[,u]*S1[,u]
  }
  marginalS1 = apply((A==1)/0.5*tau.mat1 + S1,2,mean)
  
  tau.mat0 =matrix(0,n,k)
  tau.mat0[,1] = (-1)/(S0[,1]*C0[,1])*((data$minutes_to_discharge_or_censor == tt[1] & data$event == 1) - h0[,1]*(data$minutes_to_discharge_or_censor >= tt[1]))
  for (u in 2:k){
    tau.mat0[,u] = (-1)/(S0[,u]*C0[,u])*((data$minutes_to_discharge_or_censor == tt[u] & data$event == 1) - h0[,u]*(data$minutes_to_discharge_or_censor >= tt[u])) + tau.mat0[,(u-1)]
  }
  for (u in 1:k){
    tau.mat0[,u] = tau.mat0[,u]*S0[,u]
  }
  marginalS0 = apply((A==0)/0.5*tau.mat0 + S0,2,mean)
  
  rd = marginalS1[k] - marginalS0[k]
  time.diff = c(tt[-1],t0) - tt
  rmst = sum((marginalS1 - marginalS0)*time.diff)
  
  # plot(tt,marginalS1)
  # points(tt,marginalS0,col="blue")
  
  km = survfit(Surv(data$minutes_to_discharge_or_censor,data$event)~A)
  
  kmpreds = matrix(summary(km,time = tt)$surv,ncol = 2)
  colnames(kmpreds) <- c("A0","A1")
  
  km.rd = kmpreds[k,2]-kmpreds[k,1]
  km.rmst = sum((kmpreds[,2]-kmpreds[,1])*time.diff)
  
  return(list(rd = rd,km.rd = km.rd, rmst = rmst, km.rmst = km.rmst  ))
}
individual.surv = sapply(allP, one.surv,P=P)
all.surv = one.surv(allP,P=P)
save(individual.surv,all.surv,file = paste0("./",P,"_resultsSurvival.rda"))

### relative efficiency is estimated by doing 10,000 simulations
###   and computing the empirical variance of the adjusted and unadjusted estimators
