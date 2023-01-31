#################################################################
### estimating efficiency gain of covariate adjusted analysis ###
#################################################################
library(MASS)
library(SuperLearner)
library(dplyr)
library(nleqslv)
library(mvtnorm)

source("helper_fns.R")

#################################################
### ordinal simulations: estimating gain + CI ###
#################################################

### estimate nuisance functions using super learner with the following library of learning algorithms
###    SL.gam (df=default), SL.glm, SL.glm.interaction, SL.ranger (random forest)
### we include additional gam learners in the super learner with df=5,10,and 20
###    can include other learning algorithms
learners = create.Learner("SL.gam",tune=list(deg.gam = c(5,10,20)))

### w: a data.frame of the covariates, n by p
### y: a vector of the outcome of length n
### u: monotone transform on the outcome in difference in mean, u=I for the sims
### pevent_mod, cdf_mod: n by K, conditional pdf and cdf of Y|W, estimated through proportional odds model
### link: scale on which RE estimate and CI are reported; identity, logit or log, see the example below
### cf: if true, nuisance functions are estimated using 10-fold cross fitting; default to true
### for single categorical covariate, nuisance functions can also be estimated using within group sample mean

est_var_DIM_clean = function(y,w,u,pevent_mod,cdf_mod,link, cf = TRUE, nfold=10){
  n = dim(w)[1]
  
  ### unadjusted variance
  sigma_u = var(u(y))
  var_u = mean(((u(y)-mean(u(y)))^2-sigma_u)^2)
  
  ### fully adjusted variance
  # use super learner to estimate the conditional mean of u(Y)
  if (!cf){
    mod_super = SuperLearner(Y=u(y), X=w, family=gaussian(),
                             SL.library=c("SL.gam","SL.glm","SL.glm.interaction",learners$names,"SL.ranger"))
    mean_super= as.vector(predict(mod_super, w)$pred)
  } else {
    fold.size = ceiling(n/nfold)
    fold.id = NULL
    for (i in 1:(nfold-1)){
      fold.id = c(fold.id,rep(i,fold.size))
    }
    fold.id = c(fold.id,rep(nfold,n-(nfold-1)*fold.size))
    fold.id = fold.id[sample(c(1:length(fold.id)),length(fold.id))]
    
    mean_super = rep(0,n)
    for (i in 1:nfold){
      outfoldx = w[fold.id != i,,drop=FALSE]
      infoldx  = w[fold.id == i,,drop=FALSE]
      mod_super = SuperLearner(Y=u(y[fold.id != i]), X=outfoldx, family=gaussian(),
                               SL.library=c("SL.gam","SL.glm","SL.glm.interaction",learners$names,"SL.ranger"))
      mean_super[fold.id==i] = predict(mod_super, infoldx)$pred
    }
  }
  sigma_a2 = mean((u(y)-mean_super)^2)
  var_a2 = mean(((u(y)-mean_super)^2-sigma_a2)^2)
  
  if (link == "logit"){
    releff_a2 = logit(sigma_a2/sigma_u)
    var_re_a2 = mean((((u(y)-mean_super)^2-sigma_a2)/sigma_u-((u(y)-mean(u(y)))^2-sigma_u)*sigma_a2/sigma_u^2)^2)
    var_re_a2 = var_re_a2/(expit(releff_a2)*(1-expit(releff_a2)))^2
  } else if (link == "log"){
    releff_a2 = log(sigma_a2/sigma_u)
    var_re_a2 = mean((((u(y)-mean_super)^2-sigma_a2)/sigma_u-((u(y)-mean(u(y)))^2-sigma_u)*sigma_a2/sigma_u^2)^2)
    var_re_a2 = var_re_a2/(exp(releff_a2))^2
  } else if (link == "identity"){
    releff_a2 = sigma_a2/sigma_u
    var_re_a2 = mean((((u(y)-mean_super)^2-sigma_a2)/sigma_u-((u(y)-mean(u(y)))^2-sigma_u)*sigma_a2/sigma_u^2)^2)
  }
  
  
  ### parametric adjusted variance
  mean_mod = as.numeric(u(c(1:dim(pevent_mod)[2])) %*% t(pevent_mod))
  sigma_m = mean((u(y)-mean_mod)^2)
  
  # influence function calculation
  ycat = dim(pevent_mod)[2]
  
  # partial derivative of the adjusted variance with respect to the proportional odds model coefficients
  dlayer1 = rep(0,ycat-1+dim(w)[2])
  for (i in 1:(ycat-1)){
    dlayer1[i] = mean(2*(u(y)-mean_mod)*(u(i+1)-u(i))*cdf_mod[,i]*(1-cdf_mod[,i]))
  }
  dsum = matrix(0,n,dim(w)[2])
  for (j in 1:(ycat-1)){
    dsum = dsum + (u(j+1)-u(j))*as.matrix(w)*cdf_mod[,j]*(1-cdf_mod[,j])
  }
  dlayer1[ycat:(ycat-1+dim(w)[2])] = apply(2*(u(y)-mean_mod)*dsum,2,mean)
  
  # information matrix for the proportional odds model coefficients
  dlayer2 = matrix(0,ycat-1+dim(w)[2],ycat-1+dim(w)[2])
  for (i in 1:(ycat-1)){
    dlayer2[i,i] = mean(cdf_mod[,i]*(1-cdf_mod[,i]))
    dlayer2[ycat:(ycat-1+dim(w)[2]),i] = apply(w*cdf_mod[,i]*(1-cdf_mod[,i]),2,mean)
    dlayer2[i,ycat:(ycat-1+dim(w)[2])] = dlayer2[ycat:(ycat-1+dim(w)[2]),i]
    dlayer2[ycat:(ycat-1+dim(w)[2]),ycat:(ycat-1+dim(w)[2])] = dlayer2[ycat:(ycat-1+dim(w)[2]),ycat:(ycat-1+dim(w)[2])] +  t(as.matrix(w)) %*% diag(cdf_mod[,i]*(1-cdf_mod[,i])) %*% as.matrix(w)/n #mean(w^2*cdf_mod[,i]*(1-cdf_mod[,i]))
  }
  
  # estimating functions for the proportional odds model coefficients
  IFab = matrix(0,ycat-1+dim(w)[2],n)
  for (i in 1:(ycat-1)){
    IFab[i,] = (y <= i) - cdf_mod[,i]
    IFab[ycat:(ycat-1+dim(w)[2]),] = IFab[ycat:(ycat-1+dim(w)[2]),] + t(w*((y <= i) - cdf_mod[,i]))
  }
  
  var_m = mean(((u(y)-mean_mod)^2-sigma_m + t(IFab) %*% solve(dlayer2) %*% dlayer1)^2)
  if (link == "logit"){
    releff_m = logit(sigma_m/sigma_u)
    var_re_m = mean((((u(y)-mean_mod)^2-sigma_m + t(IFab) %*% solve(dlayer2) %*% dlayer1)/sigma_u - ((u(y)-mean(u(y)))^2-sigma_u)*sigma_m/sigma_u^2)^2)
    var_re_m = var_re_m/(expit(releff_m)*(1-expit(releff_m)))^2
  } else if (link == "log"){
    releff_m = log(sigma_m/sigma_u)
    var_re_m = mean((((u(y)-mean_mod)^2-sigma_m + t(IFab) %*% solve(dlayer2) %*% dlayer1)/sigma_u - ((u(y)-mean(u(y)))^2-sigma_u)*sigma_m/sigma_u^2)^2)
    var_re_m = var_re_m/(exp(releff_m))^2
  } else if (link == "identity"){
    releff_m = sigma_m/sigma_u
    var_re_m = mean((((u(y)-mean_mod)^2-sigma_m + t(IFab) %*% solve(dlayer2) %*% dlayer1)/sigma_u - ((u(y)-mean(u(y)))^2-sigma_u)*sigma_m/sigma_u^2)^2)
  }
  
  variances = c(releff_a2,var_re_a2,releff_m,var_re_m)
  names(variances) = c("releff_a2","var_re_a2","releff_m","var_re_m")
  return(variances)
}

### Mann-Whitney
est_var_MW_clean = function(y,w,pevent_mod,cdf_mod,link, cf=TRUE, nfold=10){
  ycat = dim(pevent_mod)[2]
  n = dim(w)[1]
  
  ### unadjusted variance, first need to estimate the eta function (h function here)
  marginal = rep(0,ycat)
  h = rep(0,ycat)
  for (i in 1:ycat){
    marginal[i] = sum(y == i)/n
    h[i] = (sum(y<i) + sum(y==i)/2)/n
  }
  sigma_u = var(h[y])
  IFs = 0
  for (i in 1:(ycat-1)){
    IFs = IFs + (marginal[ycat]^2-marginal[i]^2)*((y==i)-marginal[i])/4
  }
  var_u = mean(IFs^2)
  
  ### fully adjusted variance
  # estimate the conditional mean function of eta given covariates using super learner
  hh = h[y]
  if (!cf){
    mod_super = SuperLearner(Y=hh, X=w, family=gaussian(),
                             SL.library=c("SL.gam","SL.glm","SL.glm.interaction",learners$names,"SL.ranger"))
    mean_super= as.vector(predict(mod_super, w)$pred)
  } else {
    fold.size = ceiling(n/nfold)
    fold.id = NULL
    for (i in 1:(nfold-1)){
      fold.id = c(fold.id,rep(i,fold.size))
    }
    fold.id = c(fold.id,rep(nfold,n-(nfold-1)*fold.size))
    fold.id = fold.id[sample(c(1:length(fold.id)),length(fold.id))]
    
    mean_super = rep(0,n)
    for (i in 1:nfold){
      outfoldx = w[fold.id != i,,drop=FALSE]
      infoldx  = w[fold.id == i,,drop=FALSE]
      mod_super = SuperLearner(Y=hh[fold.id != i], X=outfoldx, family=gaussian(),
                               SL.library=c("SL.gam","SL.glm","SL.glm.interaction",learners$names,"SL.ranger"))
      mean_super[fold.id==i] = predict(mod_super, infoldx)$pred
    }
  }
  sigma_a2 = mean((hh-mean_super)^2)
  
  # estimate the variance of sigma_a2 using close form expression of its influence function
  IFint = rep(0,ycat)
  for (i in 1:ycat){
    IFint[i] = mean((hh-mean_super)*((y>i)+0.5*(y==i)))
  }
  var_a2 = mean(((hh-mean_super)^2 - 3*sigma_a2 + 2*IFint[y])^2)
  if (link == "logit"){
    releff_a2 = logit(sigma_a2/sigma_u)
    var_re_a2 = mean((((hh-mean_super)^2 - 3*sigma_a2 + 2*IFint[y])/sigma_u-IFs*sigma_a2/sigma_u^2)^2)
    var_re_a2 = var_re_a2/(expit(releff_a2)*(1-expit(releff_a2)))^2
  } else if (link == "log"){
    releff_a2 = log(sigma_a2/sigma_u)
    var_re_a2 = mean((((hh-mean_super)^2 - 3*sigma_a2 + 2*IFint[y])/sigma_u-IFs*sigma_a2/sigma_u^2)^2)
    var_re_a2 = var_re_a2/(exp(releff_a2))^2
  } else if (link == "identity"){
    releff_a2 = sigma_a2/sigma_u
    var_re_a2 = mean((((hh-mean_super)^2 - 3*sigma_a2 + 2*IFint[y])/sigma_u-IFs*sigma_a2/sigma_u^2)^2)
  }
  
  ### parametrically adjusted
  hmean_mod = as.numeric(h %*% t(pevent_mod))
  sigma_ap = mean((hh-hmean_mod)^2)
  
  # partial derivative of the adjusted variance with respect to the proportional odds model coefficients
  dlayer1 = rep(0,ycat-1+dim(w)[2])
  for (i in 1:(ycat-1)){
    dlayer1[i] = mean(2*(hh-hmean_mod)*(h[i+1]-h[i])*cdf_mod[,i]*(1-cdf_mod[,i]))
  }
  dsum = 0
  for (j in 1:(ycat-1)){
    dsum = dsum + (h[j+1]-h[j])*as.matrix(w)*cdf_mod[,j]*(1-cdf_mod[,j])
  }
  dlayer1[ycat:(ycat-1+dim(w)[2])] = apply(2*(hh-hmean_mod)*dsum,2,mean)
  
  # information matrix for the proportional odds model coefficients
  dlayer2 = matrix(0,ycat-1+dim(w)[2],ycat-1+dim(w)[2])
  for (i in 1:(ycat-1)){
    dlayer2[i,i] = mean(cdf_mod[,i] * (1-cdf_mod[,i]))
    dlayer2[ycat:(ycat-1+dim(w)[2]),i] = apply(as.matrix(w) * cdf_mod[,i] * (1-cdf_mod[,i]),2,mean)
    dlayer2[i,ycat:(ycat-1+dim(w)[2])] = dlayer2[ycat:(ycat-1+dim(w)[2]),i]
    dlayer2[ycat:(ycat-1+dim(w)[2]),ycat:(ycat-1+dim(w)[2])] = dlayer2[ycat:(ycat-1+dim(w)[2]),ycat:(ycat-1+dim(w)[2])] + 
      t(as.matrix(w)) %*% diag(cdf_mod[,i]*(1-cdf_mod[,i])) %*% as.matrix(w)/n
  }
  
  # estimating functions for the proportional odds model coefficients
  IFab = matrix(0,ycat-1+dim(w)[2],n)
  for (i in 1:(ycat-1)){
    IFab[i,] = (y <= i) - cdf_mod[,i]
    IFab[ycat:(ycat-1+dim(w)[2]),] = IFab[ycat:(ycat-1+dim(w)[2]),] + t(w*((y <= i) - cdf_mod[,i]))
  }
  
  # partial derivative of the adjusted variance with respect to the marginal pdf of Y
  dp = rep(0,ycat)
  for (i in 1:ycat){
    if (i==1){
      dp[i] = mean((-1)*(hh-hmean_mod)*((y <=1)-cdf_mod[,1]))
    } else if (i == ycat){
      dp[i] = mean((-1)*(hh-hmean_mod)*((y <= (ycat-1))-cdf_mod[,ycat-1]))
    } else {
      dp[i] = mean((-1)*(hh-hmean_mod)*((y <=i)-cdf_mod[,i]+(y<=(i-1))-cdf_mod[,i-1]))
    }
  }
  
  # influence function of the marginal pdf of Y
  IFp = matrix(0,ycat,n)
  for (i in 1:ycat){
    IFp[i,] = (y == i) - marginal[i] 
  }
  
  # estimate the variance of sigma_ap using close form expression of its influence function
  var_ap = mean(((hh-hmean_mod)^2-sigma_ap + t(IFab) %*% solve(dlayer2) %*% dlayer1 + t(IFp) %*% dp)^2)
  if (link == "logit"){
    releff_ap = logit(sigma_ap/sigma_u)
    var_re_ap = mean((((hh-hmean_mod)^2-sigma_ap + t(IFab) %*% solve(dlayer2) %*% dlayer1 + t(IFp) %*% dp)/sigma_u - IFs*sigma_ap/sigma_u^2)^2)
    var_re_ap = var_re_ap/(expit(releff_ap)*(1-expit(releff_ap)))^2
  } else if (link == "log"){
    releff_ap = log(sigma_ap/sigma_u)
    var_re_ap = mean((((hh-hmean_mod)^2-sigma_ap + t(IFab) %*% solve(dlayer2) %*% dlayer1 + t(IFp) %*% dp)/sigma_u - IFs*sigma_ap/sigma_u^2)^2)
    var_re_ap = var_re_ap/(exp(releff_ap))^2
  } else if (link == "identity"){
    releff_ap = sigma_ap/sigma_u
    var_re_ap = mean((((hh-hmean_mod)^2-sigma_ap + t(IFab) %*% solve(dlayer2) %*% dlayer1 + t(IFp) %*% dp)/sigma_u - IFs*sigma_ap/sigma_u^2)^2)
  }
  
  variances = c(releff_a2,var_re_a2,releff_ap,var_re_ap)
  names(variances) = c("releff_a2","var_re_a2","releff_ap","var_re_ap")
  return(variances)
}

### Log odds ratio
est_var_LOR_clean = function(y,w,pevent_mod,cdf_mod,link, cf=TRUE, nfold=10){
  ycat = dim(pevent_mod)[2]
  n = dim(w)[1]
  
  ### unadjusted variance
  # first estimate the marginal cdf of Y
  cdf = rep(0, ycat-1)
  for (j in 1:(ycat-1)){
    cdf[j] = mean(y <=j)
  }
  sums = 0
  for (j in 1:(ycat-1)){
    sums = sums + ((y <= j)-cdf[j])/(cdf[j]*(1-cdf[j]))
  }
  sums = sums/(ycat-1)
  sigma_u = mean(sums^2)
  
  # partial derivative of sigma_u with respect to and influence function of the marginal cdf of Y
  dtheta = rep(0,ycat-1)
  IFtheta = matrix(0,n,ycat-1)
  for (i in 1:(ycat-1)){
    dtheta[i] = mean((-2)/(ycat-1)*sums*(((y<=i)-cdf[i])/(cdf[i]*(1-cdf[i])))^2)
    IFtheta[,i] = (y<=i) - cdf[i]
  }
  # close form expression of the influence function of sigma_u
  var_u = mean((sums^2-sigma_u + IFtheta %*% dtheta)^2)
  
  ### fully adjusted
  # estimate the conditional distribution function of Y given covariates using super learner
  cdf_super = matrix(0, n, ycat-1)
  if (!cf){
    for (j in 1:(ycat - 1)){
      mod_super = SuperLearner(Y=1*(y <= j), X=w, family=binomial(),
                               SL.library=c("SL.gam","SL.glm","SL.glm.interaction",learners$names,"SL.ranger"))
      cdf_super[,j] = as.vector(predict(mod_super, w)$pred)
    }
  } else {
    for (j in 1:(ycat-1)){
      fold.size = ceiling(n/nfold)
      fold.id = NULL
      for (i in 1:(nfold-1)){
        fold.id = c(fold.id,rep(i,fold.size))
      }
      fold.id = c(fold.id,rep(nfold,n-(nfold-1)*fold.size))
      fold.id = fold.id[sample(c(1:length(fold.id)),length(fold.id))]
      
      for (i in 1:nfold){
        outfoldx = w[fold.id != i,,drop=FALSE]
        infoldx  = w[fold.id == i,,drop=FALSE]
        mean.mod = SuperLearner(Y=1*(y[fold.id != i] <= j), X=outfoldx, family=binomial(),
                                SL.library=c("SL.gam","SL.glm","SL.glm.interaction",learners$names,"SL.ranger"))
        cdf_super[fold.id==i,j] = predict(mean.mod, infoldx)$pred
      }
    }
  }
  sums_super = 0
  for (j in 1:(ycat-1)){
    sums_super = sums_super + ((y <= j)-cdf_super[,j])/(cdf[j]*(1-cdf[j]))
  }
  sums_super = sums_super/(ycat-1)
  sigma_a2 = mean(sums_super^2)
  dtheta_super = rep(0,ycat-1)
  for (i in 1:(ycat-1)){
    dtheta_super[i] = mean((-2)/(ycat-1)*sums_super*((y <= i)-cdf_super[,i])*(1-2*cdf[i])/(cdf[i]*(1-cdf[i]))^2)
  }
  # estimate variance of sigma_a2 using its influence function
  var_a2 = mean((sums_super^2-sigma_a2 + IFtheta %*% dtheta_super)^2)
  if (link == "logit"){
    releff_a2 = logit(sigma_a2/sigma_u)
    var_re_a2 = mean(((sums_super^2-sigma_a2 + IFtheta %*% dtheta_super)/sigma_u-(sums^2-sigma_u + IFtheta %*% dtheta)*sigma_a2/sigma_u^2)^2)
    var_re_a2 = var_re_a2/(expit(releff_a2)*(1-expit(releff_a2)))^2
  } else if (link == "log"){
    releff_a2 = log(sigma_a2/sigma_u)
    var_re_a2 = mean(((sums_super^2-sigma_a2 + IFtheta %*% dtheta_super)/sigma_u-(sums^2-sigma_u + IFtheta %*% dtheta)*sigma_a2/sigma_u^2)^2)
    var_re_a2 = var_re_a2/(exp(releff_a2))^2
  } else if (link == "identity"){
    releff_a2 = sigma_a2/sigma_u
    var_re_a2 = mean(((sums_super^2-sigma_a2 + IFtheta %*% dtheta_super)/sigma_u-(sums^2-sigma_u + IFtheta %*% dtheta)*sigma_a2/sigma_u^2)^2)
  }
  
  ### parametric adjusted
  sums_p = 0
  for (j in 1:(ycat-1)){
    sums_p = sums_p + ((y <= j)-cdf_mod[,j])/(cdf[j]*(1-cdf[j]))
  }
  sums_p = sums_p/(ycat-1)
  sigma_ap = mean(sums_p^2)
  dlayer1 = rep(0,ycat-1+dim(w)[2])
  for (i in 1:(ycat-1)){
    dlayer1[i] = mean((-2)*sums_p/(ycat-1)/(cdf[i]*(1-cdf[i]))*cdf_mod[,i]*(1-cdf_mod[,i]))
    dlayer1[ycat:(ycat-1+dim(w)[2])] = dlayer1[ycat:(ycat-1+dim(w)[2])] + apply((-2)*sums_p/(ycat-1)/(cdf[i]*(1-cdf[i]))*as.matrix(w)*cdf_mod[,i]*(1-cdf_mod[,i]),2,mean)
  }
  
  dlayer2 = matrix(0,ycat-1+dim(w)[2],ycat-1+dim(w)[2])
  for (i in 1:(ycat-1)){
    dlayer2[i,i] = mean(cdf_mod[,i]*(1-cdf_mod[,i]))
    dlayer2[ycat:(ycat-1+dim(w)[2]),i] = apply(w*cdf_mod[,i]*(1-cdf_mod[,i]),2,mean)
    dlayer2[i,ycat:(ycat-1+dim(w)[2])] = dlayer2[ycat:(ycat-1+dim(w)[2]),i]
    dlayer2[ycat:(ycat-1+dim(w)[2]),ycat:(ycat-1+dim(w)[2])] = dlayer2[ycat:(ycat-1+dim(w)[2]),ycat:(ycat-1+dim(w)[2])] + 
      t(as.matrix(w)) %*% diag(cdf_mod[,i]*(1-cdf_mod[,i])) %*% as.matrix(w)/n
  }
  IFab = matrix(0,ycat-1+dim(w)[2],n)
  for (i in 1:(ycat-1)){
    IFab[i,] = (y <= i) - cdf_mod[,i]
    IFab[ycat:(ycat-1+dim(w)[2]),] = IFab[ycat:(ycat-1+dim(w)[2]),] + t(w*((y <= i) - cdf_mod[,i]))
  }
  dtheta_p = rep(0,ycat-1)
  for (i in 1:(ycat-1)){
    dtheta_p[i] = mean((-2)/(ycat-1)*sums_p*((y <= i)-cdf_mod[,i])*(1-2*cdf[i])/(cdf[i]*(1-cdf[i]))^2)
  }
  # estimate variance using close form expression for its influence function
  var_ap = mean((sums_p^2-sigma_ap + IFtheta %*% dtheta_p + t(IFab) %*% solve(dlayer2) %*% dlayer1)^2)
  if (link == "logit"){
    releff_ap = logit(sigma_ap/sigma_u)
    var_re_ap = mean(((sums_p^2-sigma_ap + IFtheta %*% dtheta_p + t(IFab) %*% solve(dlayer2) %*% dlayer1)/sigma_u - (sums^2-sigma_u + IFtheta %*% dtheta)*sigma_ap/sigma_u^2)^2)
    var_re_ap = var_re_ap/(expit(releff_ap)*(1-expit(releff_ap)))^2
  } else if (link == "log"){
    releff_ap = log(sigma_ap/sigma_u)
    var_re_ap = mean(((sums_p^2-sigma_ap + IFtheta %*% dtheta_p + t(IFab) %*% solve(dlayer2) %*% dlayer1)/sigma_u - (sums^2-sigma_u + IFtheta %*% dtheta)*sigma_ap/sigma_u^2)^2)
    var_re_ap = var_re_ap/(exp(releff_ap))^2
  } else if (link == "identity"){
    releff_ap = sigma_ap/sigma_u
    var_re_ap = mean(((sums_p^2-sigma_ap + IFtheta %*% dtheta_p + t(IFab) %*% solve(dlayer2) %*% dlayer1)/sigma_u - (sums^2-sigma_u + IFtheta %*% dtheta)*sigma_ap/sigma_u^2)^2)
  }
  
  variances = c(releff_a2,var_re_a2,releff_ap,var_re_ap)
  #names(variances) = c("unadjusted","adjusted_1","adjusted_2","adjusted_param")
  return(variances)
}

### example with a univariate covariate, age category with 7 categories
n=500
page = c(0.01,0.09,0.12,0.13,0.18,0.22,0.25)
pevent = matrix(c(0,0.01,0.03,0.08,0.11,0.17,0.37,0,0.18,0.32,0.31,0.37,0.47,0.35,
                  rep(0,7)),ncol=3,byrow=FALSE)
pevent[,3] = 1-pevent[,1]-pevent[,2]
dat = gen_ordinal_data(n,page,pevent)
y = dat$Y
w = data.frame("W" = dat$W)
param_mod = ord_work_param_clean(y,w)
pevent_mod = param_mod[[1]]
cdf_mod = param_mod[[2]]
result = est_var_DIM_clean(y=y, w=w, u=I, pevent_mod=pevent_mod,
                           cdf_mod=cdf_mod, link="logit",cf=TRUE,nfold=10)
# build CI on logit scale then transform
CIlow = expit(result[1] -qnorm(0.975)*sqrt(result[2]/n))
CIup = expit(result[1] +qnorm(0.975)*sqrt(result[2]/n))

expit(result[3] -qnorm(0.975)*sqrt(result[4]/n))
expit(result[3] +qnorm(0.975)*sqrt(result[4]/n))


###################
### simulations ###
###################
# the following code reproduce the simulation results in the paper
# the results were obtained by running on a cluster computer

### data generating mechanism, setting I: a univariate covariate, age category with 7 categories
page = c(0.01,0.09,0.12,0.13,0.18,0.22,0.25)
pevent = matrix(c(0,0.01,0.03,0.08,0.11,0.17,0.37,0,0.18,0.32,0.31,0.37,0.47,0.35,
                  rep(0,7)),ncol=3,byrow=FALSE)
pevent[,3] = 1-pevent[,1]-pevent[,2]

### data generating mechanism, setting II: continuous covariates, with the following Y|W distribution
p.cond = function(w1,w2){
  p1.cond = expit(1.5*w1-2*w2^2)
  p2.cond = expit(1.5*w1+w2^2)-expit(1.5*w1-2*w2^2)
  p3.cond = 1-expit(1.5*w1+w2^2)
  return(c(p1.cond, p2.cond, p3.cond))
}

### run one simulation replication
###   n = sample size, 500, 1000 or 2000
###   param = treatment effect estimand, DIM or MW or LOR
###   pcov = marginal distribution of univariate categorical covariate, only used when setting = 1
###   pout = conditional distribution of Y given univariate categorical covariate, only used when setting = 1
###   link = identity or logit or log 
###   u = monotone transformation on the outcome for DIM
###   cf = whether cross fitting is used for super learner
###   setting = 1 or 2
do.one.clean = function(n,param,pcov,pout,link,u=I,cf,nfold=10, setting){
  if (setting == 1){
    dat = gen_ordinal_data(n,pcov,pout)
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
  
  # fit the working proportional odds model, see the helper functions for details
  # output two n by K matrices
  param_mod = ord_work_param_clean(y,w)
  pevent_mod = param_mod[[1]]
  cdf_mod = param_mod[[2]]
  
  if (param == "DIM"){
    result = est_var_DIM_clean(y,w,u,pevent_mod,cdf_mod,link,cf,nfold)
  } else if (param == "MW"){
    result = est_var_MW_clean(y,w,pevent_mod,cdf_mod,link,cf,nfold)
  } else if (param == "LOR"){
    result = est_var_LOR_clean(y,w,pevent_mod,cdf_mod,link,cf,nfold)
  }
  return(result)
}

### many sims run using cluster computer, each job run one simulation
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

set.seed(9*seed)
res = replicate(1,do.one.clean(nn,"LOR",page,pevent,"logit",I,TRUE,setting=1))
filename = paste("seed",seed,"n",nn,"DIM.RData",sep="")
save(res,file=filename)

### compute the true RE based on true data generating mechanism, Setting I
calc_true_gain = function(pcov,pout,param,u=I){
  wcat = length(pcov)
  ycat = dim(pout)[2]
  marginal = pcov %*% pout
  cdf = c(marginal[1],rep(0,ycat-1))
  for (i in 2:ycat){
    cdf[i] = sum(marginal[1:i])
  }
  
  mod_coef = nleqslv(x=rep(0,dim(pout)[2]),fn=prop_odd_EE,pcov=pcov,pout=pout)$x
  pout_mod = matrix(0,wcat,ycat)
  for (i in 1:ycat){
    if (i == 1){
      pout_mod[,i] = expit(mod_coef[i]+mod_coef[ycat]*c(1:wcat))
    } else if (i < ycat & i>1){
      pout_mod[,i] = expit(mod_coef[i]+mod_coef[ycat]*c(1:wcat)) - pout_mod[,i-1]
    } else {
      pout_mod[,i] = 1 - rowSums(pout_mod[,1:i-1])
    }
  }
  
  if (param == "DIM"){
    m = sum(u(1:ycat)*marginal)
    sigma_u = sum((u(1:ycat)-m)^2*marginal)
    
    condvar = rep(0,wcat)
    condm = rep(0,wcat)
    for (i in 1:wcat){
      condm[i] = sum(u(1:ycat)*pout[i,])
      condvar[i] = sum((u(1:ycat)-condm[i])^2*pout[i,])
    }
    sigma_a = sum(condvar*pcov)
    
    #print(t(pout_mod)%*%pcov)
    condm_mod = u(1:ycat)%*% t(pout_mod)
    sigma_ap = sigma_a + sum((condm-condm_mod)^2*pcov)
  }
  if (param == "MW"){
    h = cdf - 0.5*marginal
    hmean = sum(h*marginal)
    sigma_u = sum((h-hmean)^2*marginal)
    
    cond_hmean = rep(0,wcat)
    cond_hvar = rep(0,wcat)
    for (i in 1:wcat){
      cond_hmean[i] = sum(h*pout[i,])
      cond_hvar[i] = sum((h-cond_hmean[i])^2*pout[i,])
    }
    sigma_a = sum(cond_hvar*pcov)
    
    cond_hmean_mod = h %*% t(pout_mod)
    sigma_ap = sigma_a + sum((cond_hmean-cond_hmean_mod)^2*pcov)
  }
  if (param == "LOR"){
    f = 1/cdf[-ycat]/(1-cdf[-ycat])
    h = rep(0,ycat)
    for (i in 1:(ycat-1)){
      h[i] = sum(f[i:(ycat-1)])
    }
    hmean = sum(h*marginal)
    sigma_u = sum((h-hmean)^2*marginal)/(ycat-1)^2
    cond_hmean = rep(0,wcat)
    cond_hvar = rep(0,wcat)
    for (i in 1:wcat){
      cond_hmean[i] = sum(h*pout[i,])
      cond_hvar[i] = sum((h-cond_hmean[i])^2*pout[i,])
    }
    sigma_a = sum(cond_hvar*pcov)/(ycat-1)^2
    
    cond_hmean_mod = h %*% t(pout_mod)
    sigma_ap = sigma_a + sum((cond_hmean-cond_hmean_mod)^2*pcov)/(ycat-1)^2
  }
  return(c(sigma_u,sigma_a,sigma_ap,sigma_a/sigma_u,sigma_ap/sigma_u))
}
truth = calc_true_gain(page,pevent,"DIM")

### compute the truth based on true data generating mechanism, Setting II
### here population quantities involve integrals with respect to a 3-d gaussian density
### such integrals are computed using the cubature package, hcubature function
library(cubature)
p1marg = hcubature(function(v){p.cond(v[1],v[2])[1]*dnorm(v[1])*0.5},c(-100,-1),c(100,1))$integral
p2marg = hcubature(function(v){p.cond(v[1],v[2])[2]*dnorm(v[1])*0.5},c(-100,-1),c(100,1))$integral
p3marg = hcubature(function(v){p.cond(v[1],v[2])[3]*dnorm(v[1])*0.5},c(-100,-1),c(100,1))$integral
u=I
mean.marg = u(1)*p1marg + u(2)*p2marg + u(3)*p3marg
var.marg = (u(1)-mean.marg)^2*p1marg + (u(2)-mean.marg)^2*p2marg + (u(3)-mean.marg)^2*p3marg
rm(u)

cond.var = function(x,y, u=I){
  mean.cond = sum(c(u(1),u(2),u(3))*p.cond(x,y))
  var.cond = sum((c(u(1),u(2),u(3)) - mean.cond)^2*p.cond(x,y))
  return(var.cond)
}
mean.cond.var = hcubature(function(v){cond.var(v[1],v[2])*dnorm(v[1])*0.5},c(-100,-1),c(100,1))$integral

#true.coef = optim(par=c(-0.627,  0.32376,  1.44255,  0, 0),fn=prop.odd.obj, ycat=3, cutoff=20)$par
true.coef = c(-0.6269274017, 0.3237308977, 1.4408499574, -0.0002779427, 0.0001601257)


cond.var.mod = function(x,y,z,u=I){
  mean.cond.mod = u(1)*expit(true.coef[1]+sum(true.coef[3:5]*c(x,y,z))) +
    u(2)*(expit(true.coef[2]+sum(true.coef[3:5]*c(x,y,z)))-expit(true.coef[1]+sum(true.coef[3:5]*c(x,y,z))))+
    u(3)*(1-expit(true.coef[2]+sum(true.coef[3:5]*c(x,y,z))))
  var.cond.mod = sum(c((u(1) - mean.cond.mod)^2,(u(2)-mean.cond.mod)^2,(u(3)-mean.cond.mod)^2)*p.cond(x,y))
  return(var.cond.mod)
}

mean.mod.var = hcubature(function(v){cond.var.mod(v[1],v[2],v[3])*dmvnorm(c(v[1],v[3]),mean=c(0,0), sigma=matrix(c(1,0.3,0.3,1),ncol=2))*0.5},
                         c(-100,-1,-100),c(100,1,100))$integral

re.DIM = mean.cond.var/var.marg
re.DIM.mod = mean.mod.var/var.marg

eta1 = p1marg/2
eta2 = p1marg + p2marg/2
eta3 = p1marg + p2marg + p3marg/2
mean.eta = sum(c(eta1,eta2,eta3)*c(p1marg,p2marg,p3marg))
var.eta = sum((c(eta1,eta2,eta3)-mean.eta)^2*c(p1marg,p2marg,p3marg))

var.marg.eta = (1 - sum(c(p1marg,p2marg,p3marg)^3))/12

cond.var.eta = function(x,y){
  mean.cond = sum(c(eta1,eta2,eta3)*p.cond(x,y))
  var.cond = sum((c(eta1,eta2,eta3)-mean.cond)^2*p.cond(x,y))
  return(var.cond)
}
mean.cond.var.eta = hcubature(function(v){cond.var.eta(v[1],v[2])*dnorm(v[1])*0.5},c(-100,-1),c(100,1))$integral

cond.var.mod.eta = function(x,y,z){
  mean.cond.mod = eta1*expit(true.coef[1]+sum(true.coef[3:5]*c(x,y,z))) +
    eta2*(expit(true.coef[2]+sum(true.coef[3:5]*c(x,y,z)))-expit(true.coef[1]+sum(true.coef[3:5]*c(x,y,z))))+
    eta3*(1-expit(true.coef[2]+sum(true.coef[3:5]*c(x,y,z))))
  var.cond.mod = sum(c((eta1 - mean.cond.mod)^2,(eta2-mean.cond.mod)^2,(eta3-mean.cond.mod)^2)*p.cond(x,y))
  return(var.cond.mod)
}

mean.mod.var.eta = hcubature(function(v){cond.var.mod.eta(v[1],v[2],v[3])*dmvnorm(c(v[1],v[3]),mean=c(0,0), sigma=matrix(c(1,0.3,0.3,1),ncol=2))*0.5},
                             c(-100,-1,-100),c(100,1,100))$integral
re.MW = mean.cond.var.eta / var.marg.eta
re.MW.mod = mean.mod.var.eta / var.marg.eta

F.marg = cumsum(c(p1marg,p2marg,p3marg))
zeta1 = 0.5*(1/(F.marg[1]*(1-F.marg[1])) + 1/(F.marg[2]*(1-F.marg[2])))
zeta2 = 0.5*(1/(F.marg[2]*(1-F.marg[2])))
zeta3 = 0
mean.zeta = sum(c(zeta1,zeta2,zeta3)*c(p1marg,p2marg,p3marg))
var.zeta = sum((c(zeta1,zeta2,zeta3) - mean.zeta)^2*c(p1marg,p2marg,p3marg))

cond.var.zeta = function(x,y){
  mean.cond = sum(c(zeta1,zeta2,zeta3)*p.cond(x,y))
  var.cond = sum((c(zeta1,zeta2,zeta3)-mean.cond)^2*p.cond(x,y))
  return(var.cond)
}
mean.cond.var.zeta = hcubature(function(v){cond.var.zeta(v[1],v[2])*dnorm(v[1])*0.5},c(-100,-1),c(100,1))$integral

cond.var.mod.zeta = function(x,y,z){
  mean.cond.mod = zeta1*expit(true.coef[1]+sum(true.coef[3:5]*c(x,y,z))) +
    zeta2*(expit(true.coef[2]+sum(true.coef[3:5]*c(x,y,z)))-expit(true.coef[1]+sum(true.coef[3:5]*c(x,y,z))))+
    zeta3*(1-expit(true.coef[2]+sum(true.coef[3:5]*c(x,y,z))))
  var.cond.mod = sum(c((zeta1 - mean.cond.mod)^2,(zeta2-mean.cond.mod)^2,(zeta3-mean.cond.mod)^2)*p.cond(x,y))
  return(var.cond.mod)
}

mean.mod.var.zeta = hcubature(function(v){cond.var.mod.zeta(v[1],v[2],v[3])*dmvnorm(c(v[1],v[3]),mean=c(0,0), sigma=matrix(c(1,0.3,0.3,1),ncol=2))*0.5},
                              c(-100,-1,-100),c(100,1,100))$integral

re.LOR = mean.cond.var.zeta / var.zeta
re.LOR.mod = mean.mod.var.zeta / var.zeta