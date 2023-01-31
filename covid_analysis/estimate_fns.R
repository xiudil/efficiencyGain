#############################################################################
### estimating efficiency gain of covariate adjusted analysis: analytical ###
#############################################################################

# this file contains all functions to estimate RE for ordinal and survival outcomes for easy import
# for detailed documentation of these functions, please refer to our code for the simulation studies

### DIM
est_var_DIM_clean = function(y, w, u, pevent_mod, cdf_mod, link, cf = TRUE, nfold=10){
  n = dim(w)[1]
  
  ### unadjusted variance
  sigma_u = var(u(y))
  var_u = mean(((u(y)-mean(u(y)))^2-sigma_u)^2)
  
  ### fully adjusted variance
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
  
  ycat = dim(pevent_mod)[2]
  dlayer1 = rep(0,ycat-1+dim(w)[2])
  for (i in 1:(ycat-1)){
    dlayer1[i] = mean(2*(u(y)-mean_mod)*(u(i+1)-u(i))*cdf_mod[,i]*(1-cdf_mod[,i]))
  }
  dsum = matrix(0,n,dim(w)[2])
  for (j in 1:(ycat-1)){
    dsum = dsum + (u(j+1)-u(j))*as.matrix(w)*cdf_mod[,j]*(1-cdf_mod[,j])
  }
  dlayer1[ycat:(ycat-1+dim(w)[2])] = apply(2*(u(y)-mean_mod)*dsum,2,mean)
  dlayer2 = matrix(0,ycat-1+dim(w)[2],ycat-1+dim(w)[2])
  for (i in 1:(ycat-1)){
    dlayer2[i,i] = mean(cdf_mod[,i]*(1-cdf_mod[,i]))
    dlayer2[ycat:(ycat-1+dim(w)[2]),i] = apply(w*cdf_mod[,i]*(1-cdf_mod[,i]),2,mean)
    dlayer2[i,ycat:(ycat-1+dim(w)[2])] = dlayer2[ycat:(ycat-1+dim(w)[2]),i]
    dlayer2[ycat:(ycat-1+dim(w)[2]),ycat:(ycat-1+dim(w)[2])] = dlayer2[ycat:(ycat-1+dim(w)[2]),ycat:(ycat-1+dim(w)[2])] +  t(as.matrix(w)) %*% diag(cdf_mod[,i]*(1-cdf_mod[,i])) %*% as.matrix(w)/n #mean(w^2*cdf_mod[,i]*(1-cdf_mod[,i]))
  }
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
  names(variances) = c("RE_f","var_RE_f","RE_m","var_RE_m")
  return(variances)
}

### Mann-Whitney
est_var_MW_clean = function(y,w,pevent_mod,cdf_mod,link, cf=TRUE, nfold=10){
  ycat = dim(pevent_mod)[2]
  n = dim(w)[1]
  
  ### unadjusted variance
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
  dlayer1 = rep(0,ycat-1+dim(w)[2])
  for (i in 1:(ycat-1)){
    dlayer1[i] = mean(2*(hh-hmean_mod)*(h[i+1]-h[i])*cdf_mod[,i]*(1-cdf_mod[,i]))
  }
  dsum = 0
  for (j in 1:(ycat-1)){
    dsum = dsum + (h[j+1]-h[j])*as.matrix(w)*cdf_mod[,j]*(1-cdf_mod[,j])
  }
  dlayer1[ycat:(ycat-1+dim(w)[2])] = apply(2*(hh-hmean_mod)*dsum,2,mean)
  
  dlayer2 = matrix(0,ycat-1+dim(w)[2],ycat-1+dim(w)[2])
  for (i in 1:(ycat-1)){
    dlayer2[i,i] = mean(cdf_mod[,i] * (1-cdf_mod[,i]))
    dlayer2[ycat:(ycat-1+dim(w)[2]),i] = apply(as.matrix(w) * cdf_mod[,i] * (1-cdf_mod[,i]),2,mean)
    dlayer2[i,ycat:(ycat-1+dim(w)[2])] = dlayer2[ycat:(ycat-1+dim(w)[2]),i]
    dlayer2[ycat:(ycat-1+dim(w)[2]),ycat:(ycat-1+dim(w)[2])] = dlayer2[ycat:(ycat-1+dim(w)[2]),ycat:(ycat-1+dim(w)[2])] + 
      t(as.matrix(w)) %*% diag(cdf_mod[,i]*(1-cdf_mod[,i])) %*% as.matrix(w)/n
  }
  IFab = matrix(0,ycat-1+dim(w)[2],n)
  for (i in 1:(ycat-1)){
    IFab[i,] = (y <= i) - cdf_mod[,i]
    IFab[ycat:(ycat-1+dim(w)[2]),] = IFab[ycat:(ycat-1+dim(w)[2]),] + t(w*((y <= i) - cdf_mod[,i]))
  }
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
  IFp = matrix(0,ycat,n)
  for (i in 1:ycat){
    IFp[i,] = (y == i) - marginal[i] 
  }
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
  names(variances) = c("RE_f","var_RE_f","RE_m","var_RE_m")
  return(variances)
}

### Log odds ratio
est_var_LOR_clean = function(y,w,pevent_mod,cdf_mod,link, cf=TRUE, nfold=10){
  ycat = dim(pevent_mod)[2]
  n = dim(w)[1]
  
  ### unadjusted variance
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
  
  dtheta = rep(0,ycat-1)
  IFtheta = matrix(0,n,ycat-1)
  for (i in 1:(ycat-1)){
    dtheta[i] = mean((-2)/(ycat-1)*sums*(((y<=i)-cdf[i])/(cdf[i]*(1-cdf[i])))^2)
    IFtheta[,i] = (y<=i) - cdf[i]
  }
  var_u = mean((sums^2-sigma_u + IFtheta %*% dtheta)^2)
  
  ### fully adjusted
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
  names(variances) = c("RE_f","var_RE_f","RE_m","var_RE_m")
  return(variances)
}

calc_shc = function(dat,tt,method){
  k = length(tt)
  n = dim(dat)[1]
  S = matrix(0,n,k)
  h = S
  C = S
  H = S
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
    mod_super = survSuperLearner(time = dat$time, event = dat$event, X = data.frame("W"=dat$W), newX = data.frame("W"=dat$W), 
                                 new.times = tt, event.SL.library = event.SL.library, cens.SL.library = cens.SL.library, verbose = FALSE)
    C = mod_super$cens.SL.predict
    S = mod_super$event.SL.predict
    h = 1 - S / cbind(rep(1,n),S[,-k])
  } else if (method == "supercf"){
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

### RD and RR
# margg: marginal distribution of C, censoring time, in the future trial, user-specified
# g: conditional distribution of C|W in the future trial, user-specified
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

### RMST
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