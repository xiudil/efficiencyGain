############################################################
### helper functions for estimating RE (ordinal outcome) ###
############################################################


###########################################
### logit and expit functions
###########################################

logit = function(x){
  return(log(x/(1-x)))
}
expit = function(x){
  return(exp(x)/(1+exp(x)))
}

################################################
### functions for proportional odds model 
################################################

### find true coefficients in proportional odds model by solving estimating equations (EE)
###    coef = candidate coefficient value for the proportional odds model
###    pcov = probability distribution function of the univariate covariate
###    pout = conditional distribution of outcome given covariate, a matrix
###    output the values of the estimating equation corresponding to the candidate coefficient value (we then find coef such that EE=0)
# only used for setting = 1 to calculate true gain
prop_odd_EE = function(coef,pcov,pout){
  marginal = pcov %*% pout
  ycat = length(marginal)
  wcat = length(pcov)
  ee = rep(0,ycat)
  for (i in 1:(ycat-1)){
    ee[i] = sum(marginal[1:i]) - sum(expit(coef[i]+coef[ycat]*c(1:wcat))*pcov)
    if (i==1){
      ee[ycat] = ee[ycat] + sum(pout[,1]*c(1:wcat)*pcov) - sum(c(1:wcat)*expit(coef[i]+coef[ycat]*c(1:wcat))*pcov)
    } else{
      ee[ycat] = ee[ycat] + sum(rowSums(pout[,1:i])*c(1:wcat)*pcov) - sum(c(1:wcat)*expit(coef[i]+coef[ycat]*c(1:wcat))*pcov)
    }
  }
  return(ee)
}


### find true coefficients in proportional odds model by maximizing the objective function
### this function takes in a set of coefficients and outputs the corresponding objective
#prop_odd_obj = function(coef, ycat, cutoff){
#  obj = 0
#  for (i in 1:(ycat-1)){
#    thisint = hcubature(function(v){(sum(p.cond(v[1],v[2])[1:i])*log(expit(coef[i] + sum(coef[ycat:length(coef)]*v)))
#                                     +(1-sum(p.cond(v[1],v[2])[1:i]))*log(1-expit(coef[i] + sum(coef[ycat:length(coef)]*v))))*dmvnorm(c(v[1],v[3]), mean=c(0,0),sigma=matrix(c(1,0.3,0.3,1),ncol=2))*0.5},
#                        c(-cutoff,-1,-cutoff),c(cutoff,1,cutoff))$integral
#    #print(thisint)
#    obj = obj + thisint
#  }
#  return((-1)*obj)
#}

#########################
### the loss function of the proportional odds EE can be recast into one correspond,
#   by stacking the outcomes Y <= k
#   therefore, we fit proportional odds model by recasting into a logistic regression
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

# this function uses the above function to fit prop odds model to a dataset
# and outputs the estimated cdf and pmf conditional on covariate, n by K
#   neww is the covariate values to obtain predictions on, default to the covariate values used to fit the model
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

### generate data for setting I
gen_ordinal_out = function(i,covvec,ycat,pout){
  y = sample(ycat,size=1,prob=pout[covvec[i],])
  return(y)
}
gen_ordinal_data = function(n,pcov,pout){
  wcat = length(pcov)
  ycat = dim(pout)[2]
  
  w = sample(wcat,size = n, replace=TRUE, prob=pcov)
  y = sapply(1:n,gen_ordinal_out,covvec=w, ycat=ycat,pout=pout)
  
  dat = data.frame(w,y)
  colnames(dat) = c("W","Y")
  return(dat)
}
