####################
# COVID-19 Data    #
# 05/07/2022       #
####################
rm(list =ls())
library(table1)
library(MASS)
library(SuperLearner)
library(tidyverse)
library(nleqslv)
library(survival)
library(glmnet)
library(MASS)
library(survival)
library(survSuperLearner)


learners = create.Learner("SL.gam",tune=list(deg.gam = c(5,10,20)))
event.SL.library <- cens.SL.library <- c("survSL.km", "survSL.coxph", "survSL.gam")#, "survSL.rfsrc")


source("helper_fns.R")
source("covid_code_/estimate_fns.R")

######################
###  preprocessing ###
######################
load("mockdata.rda")
data$dm <- ifelse(data$type_1_dm=="Y" | data$type_2_dm=="Y", "Y","N" )
data$dm = as.factor(data$dm)


dim(data)
for (i in which(data$hospitalization_no>1)){
  print(data[i,"patient_died"] != data[i-1,"patient_died"])
}

# 0: censor, 1: discharge, 2: incubatation/ventillation, 3:death
data$y_ord <- ifelse(data$right_censor=="Y",0,NA)
data$y_ord <- ifelse(!is.na(data$minutes_to_intubation) |!is.na(data$minutes_to_ventilation),2,data$y_ord)
data$y_ord <- ifelse(data$patient_died=="Y",3,data$y_ord)
data$y_ord <- ifelse(is.na(data$y_ord),1,data$y_ord)
#subset(data,select = c("minutes_to_death","right_censor","minutes_to_intubation","minutes_to_ventilation","y_ord"))
data <- subset(data,y_ord!=0)
nrow(data)

count = as.data.frame(table(data$pid))
dupID2 = count$Var[count$Freq==2]
dupID3 = count$Var[count$Freq==3]
dupID4 = count$Var[count$Freq==4]

subset(look,select = c(pid,y_ord))
change2 = matrix(data$y_ord[which(data$pid %in% dupID2)],ncol=2 ,byrow=T )
rownames(change2) = dupID2
idchange2 = subset(change2, change2[,1]!=change2[,2])


change3 = matrix(data$y_ord[which(data$pid %in% dupID3)],ncol=3 ,byrow=T )
rownames(change3) = dupID3
idchange3 = subset(change3, change3[,3]!=change3[,2])


# take the more serious one and unique data
unique.data = data
unique.data$y_ord[unique.data$pid %in% rownames(idchange2)]  = 2
unique.data$y_ord[unique.data$pid %in% rownames(idchange3)]  = 2

final =unique.data[!duplicated(unique.data$pid), ]
data = final

# discretize age
data$age <- as.numeric(paste(data$age_at_admit))
data$age[is.na(data$age)] <- 90
data$age <- as.numeric(cut(data$age, seq(20,90,length.out=8),include.lowest = T))

# unfactor others:
data$Gender <- as.numeric(as.factor(data$Gender)) # 2--male
data$race <- as.numeric(as.factor(data$race)) # indian/alaskan (1) asian (2),black(3), native(4), unknown (5),white (6)
data$ethnicity <- as.numeric(as.factor(data$ethnicity)) # hispnaic(1), NH (2), unknown(3)
data$cardiovascular_disease <- as.numeric(as.factor(data$cardiovascular_disease))
data$hypertension <- as.numeric(as.factor(data$hypertension))
data$dm <- as.numeric(data$dm)
data$chronic_kidney_disease <- as.numeric(as.factor(data$chronic_kidney_disease))
data$cholesterol_meds <- as.numeric(as.factor(data$cholesterol_meds ))
data$hypertension_meds <- as.numeric(as.factor(data$hypertension_meds))
data$bmi.dis <- as.numeric(cut(data$bmi, seq(10,80,length.out=8),include.lowest = T))
data$event <- ifelse(data$patient_died=="Y",1,ifelse(data$patient_died=="N",0,NA))



# -------------ordinal----------------------------------------------------------
oneW <- function(which.W,u=I,link="logit",param){
  set.seed(1)
  if (which.W=="bmi"){
    dat <- subset(data,select = c(which.W,"y_ord"))
    dat <- dat[complete.cases(dat),]
  }else{
    dat <- subset(data,select = c(which.W,"y_ord"))
  }
  colnames(dat) <- c("W","Y")
  param_mod <-  ord_work_param_clean(dat$Y,data.frame(dat$W))
  pevent_mod <-  param_mod[[1]]
  cdf_mod <-  param_mod[[2]]
  if (param == "DIM"){
    result = est_var_DIM_clean(dat$Y,data.frame(dat$W),u,pevent_mod,cdf_mod,link,cf=F,nfold=10)
  } else if (param == "MW"){
    result = est_var_MW_clean(dat,pevent_mod,cdf_mod,link)
  } else if (param == "LOR"){
    result = est_var_LOR_clean(dat,pevent_mod,cdf_mod,link)
  }
  return(result)
}


#---------Results---------------------------------------------------------------
result_age <- sapply(c("DIM"), oneW,u=I,link="identity",which.W="age")
result_gender <- sapply(c("DIM"), oneW,u=I,link="identity",which.W="Gender")
result_race <- sapply(c("DIM"), oneW,u=I,link="identity",which.W="race")
result_card <- sapply(c("DIM"), oneW,u=I,link="identity",which.W="cardiovascular_disease")
result_hyper <- sapply(c("DIM"), oneW,u=I,link="identity",which.W="hypertension")
result_dm <- sapply(c("DIM"), oneW,u=I,link="identity",which.W="dm")
result_kidney <- sapply(c("DIM"), oneW,u=I,link="identity",which.W="chronic_kidney_disease")
result_chole_meds <- sapply(c("DIM"), oneW,u=I,link="identity",which.W="cholesterol_meds")
result_hyper_meds <- sapply(c("DIM"), oneW,u=I,link="identity",which.W="hypertension_meds")
result_bmi <- sapply(c("DIM"), oneW,u=I,link="identity",which.W="bmi")  # missing bmi data 



#---------Multiple W and use SuperLearner---------------------------------------
moreW <- function(which.W,u=I,link="logit",param){
  set.seed(1)
  dat= subset(data,select = c(which.W,"y_ord"))
  colnames(dat) = c(paste0("W",1:length(which.W)),"Y")
  param_mod <-  ord_work_param_clean(dat$Y,subset(data,select = which.W))
  pevent_mod <-  param_mod[[1]]
  cdf_mod <-  param_mod[[2]]
  if (param == "DIM"){
    result = est_var_DIM_clean(dat$Y,subset(data,select = which.W),u,pevent_mod,cdf_mod,link,cf=F,nfold=10)
  } else if (param == "MW"){
    result = est_var_MW_clean(dat,pevent_mod,cdf_mod,link)
  } else if (param == "LOR"){
    result = est_var_LOR_clean(dat,pevent_mod,cdf_mod,link)
  }
  return(result)
}

result_all <- sapply(c("DIM"), moreW,u=I,link="identity",which.W=c("age","Gender","cardiovascular_disease",
                                                                              "race","hypertension","dm","chronic_kidney_disease","cholesterol_meds","hypertension_meds"))
result_all 

################################################################################

#----------------------survival outcome-----------------------------------------
fit <- survfit(Surv(minutes_to_discharge_or_censor, event) ~ 1, data = data) 
plot(fit)

oneW.surv <-  function(which.W,t0,method){
  if (which.W=="bmi.dis"){
    data <- data[!is.na(subset(data,select=which.W)),]
  }else{
    data <- data
  }
  W = subset(data,select=which.W)
  summary(W)
  time = data$minutes_to_discharge_or_censor
  event = data$event
  tt = sort(unique(time[event ==1 & time<=t0]))
  
  dat = data.frame(cbind(time,event,W,1-event))
  colnames(dat) = c("time","event","W","censor")
  ### specify censoring distribution
  set.seed(1)
  g = matrix(rep(1-pexp((tt-1)*0.2,rate=0.00001),nrow(dat)),byrow=TRUE,ncol=length(tt))
  margg <- g[1,]
  
  
  estimates = calc_shc(dat,tt,method )
  S = estimates[[1]];h = estimates[[2]]; C = estimates[[3]];rm(estimates)
  result = est_gain_RD(dat,tt,S,h,C,margg,g)
  return(result)
}



# 70% survival at time 21000
result_age <- sapply(c("super"),oneW.surv,t0=21000,which.W="age")
result_gender <- sapply(c("super"),oneW.surv,t0=21000,which.W="Gender")
result_race <- sapply(c("super"),oneW.surv,t0=21000,which.W="race")
result_ethnicity <- sapply(c("super"),oneW.surv,t0=21000,which.W="ethnicity")
result_card <- sapply(c("super"),oneW.surv,t0=21000,which.W="cardiovascular_disease")
result_hyper <- sapply(c("super"),oneW.surv,t0=21000,which.W="hypertension")
result_dm <- sapply(c("super"),oneW.surv,t0=21000,which.W="dm")
result_kidney <- sapply(c("super"),oneW.surv,t0=21000,which.W="chronic_kidney_disease")
result_chole_meds <- sapply(c("super"),oneW.surv,t0=21000,which.W="cholesterol_meds")
result_hyper_meds <- sapply(c("super"),oneW.surv,t0=21000,which.W="hypertension_meds")
result_bmi <- sapply(c("super"),oneW.surv,t0=21000,which.W="bmi.dis") # discretize bmi


#more W
calc_shc_M = function(dat,tt,method){
  k = length(tt)
  n = dim(dat)[1]
  S = matrix(0,n,k)
  h = S
  C = S
  H = S
  if (method=="super"){
    mod_super = survSuperLearner(time = dat$time, event = dat$event, X = subset(dat,select  = which.W), newX = subset(dat,select =which.W), 
                                 new.times = tt, event.SL.library = event.SL.library, cens.SL.library = cens.SL.library, verbose = FALSE)
    C = mod_super$cens.SL.predict
    S = mod_super$event.SL.predict
    h = 1 - S / cbind(rep(1,n),S[,-k])
  }else if (method == "coxph"){
    # use elastic net to perform variable selection
    dat.complete = dat[complete.cases(dat),]
    dat.less = dat.complete[-which(dat.complete$time==0),]
    W.less = subset(dat.less,select = which.W)
    # remove time ==0 first (for glmnet to work)--only one entry
    cvfit = cv.glmnet(as.matrix(W.less),Surv(dat.less$time,dat.less$event),family = "cox")
    coef.min = coef(cvfit, s = cvfit$lambda.min)
    
    form = paste0("Surv(time,event)~", paste0(which.W[which(coef.min!=0)],collapse = "+"))
    form2 = paste0("+", paste0("I(",which.W[which(coef.min!=0)],"^2)",collapse = "+"))
    form3 = paste0("+", paste0("I(",which.W[which(coef.min!=0)],"^3)",collapse = "+"))
    form4 = paste0("+", paste0("I(",which.W[which(coef.min!=0)],"^4)",collapse = "+"))
    form5 = paste0("+", paste0("I(",which.W[which(coef.min!=0)],"^5)",collapse = "+"))
    form6 = paste0("+", paste0("I(",which.W[which(coef.min!=0)],"^6)",collapse = "+"))
    form7 = paste0("+", paste0("I(",which.W[which(coef.min!=0)],"^7)",collapse = "+"))
    
    mod_lm = coxph(as.formula(form),data = dat)
    mod_poly3 = coxph(as.formula(paste0(form,form2,form3)),data = dat)
    mod_poly5 = coxph(as.formula(paste0(form,form2,form3,form4,form5)),data = dat)
    mod_poly7 = coxph(as.formula(paste0(form,form2,form3,form4,form5,form6,form7)),data = dat)
    mods = list(mod_lm,mod_poly3,mod_poly5,mod_poly7)
    index = which.min(BIC(mod_lm,mod_poly3,mod_poly5,mod_poly7)$BIC)
    
    formC = paste0("Surv(time,censor)~", paste0(which.W[which(coef.min!=0)],collapse = "+"))
    mod_lm_c = coxph(as.formula(formC),data = dat)
    mod_poly3_c = coxph(as.formula(paste0(formC,form2,form3)),data = dat)
    mod_poly5_c = coxph(as.formula(paste0(formC,form2,form3,form4,form5)),data = dat)
    mod_poly7_c = coxph(as.formula(paste0(formC,form2,form3,form4,form5,form6,form7)),data = dat)
    mods_c = list(mod_lm_c,mod_poly3_c,mod_poly5_c,mod_poly7_c)
    index_c = which.min(BIC(mod_lm_c,mod_poly3_c,mod_poly5_c,mod_poly7_c)$BIC)
    
    for (i in 1:k){
      newdata = data.frame(cbind(rep(tt[i],n),rep(1,n),W[,which(coef.min!=0)],rep(1,n)))
      colnames(newdata)[c(1,2,ncol(newdata))] = c("time","event","censor")
      pred = predict(mods[[index]],type="expected",newdata = newdata)
      S[,i] = exp(-pred)
      H[,i] = pred
      pred_c = predict(mods_c[[index_c]],type="expected",newdata = newdata)
      C[,i] = exp(-pred_c)
      
    }
    C = cbind(rep(1,n),C[,-k])
    h = H - cbind(rep(0,n),H[,-k])
    
  } 
  return(list(S,h,C))
}
moreW.surv <-  function(which.W,t0,method){
  if(any(which.W=="bmi.dis")){
    data <- data[!is.na(data$bmi.dis),]
  } else {
    data <- data
  }
  W <<- subset(data,select=which.W)
  time = data$minutes_to_discharge_or_censor
  event = data$event
  tt = sort(unique(time[event ==1 & time<=t0]))
  
  dat = data.frame(cbind(time,event,W,1-event))
  colnames(dat)[ncol(dat)] = "censor"
  ### specify censoring distribution
  set.seed(1)
  g = matrix(rep(1-pexp((tt-1)*0.2,rate=0.00001),nrow(dat)),byrow=TRUE,ncol=length(tt))
  margg <- g[1,]
  
  estimates = calc_shc_M(dat,tt,method )
  S = estimates[[1]];h = estimates[[2]]; C = estimates[[3]];rm(estimates)
  result = est_gain_RD(dat,tt,S,h,C,margg,g)
  return(result)
}
which.W = colnames(data)[c(3,5,6,18,19,22,23,24,25,27)]
which.W[2] <- "bmi.dis"

result_all <- sapply(c("super"),moreW.surv ,t0=21000,which.W=which.W)


oneW.surv.RMST <-  function(which.W,t0,method){
  if (which.W=="bmi.dis"){
    data <- data[!is.na(subset(data,select=which.W)),]
  }else{
    data <- data
  }
  W = subset(data,select=which.W)
  time = data$minutes_to_discharge_or_censor
  event = data$event
  tt = sort(unique(time[event ==1 & time<=t0]))
  
  dat = data.frame(cbind(time,event,W,1-event))
  colnames(dat) = c("time","event","W","censor")
  set.seed(1)
  g = matrix(rep(1-pexp((tt-1)*0.2,rate=0.00001),nrow(dat)),byrow=TRUE,ncol=length(tt))
  margg <- g[1,]
  estimates = calc_shc(dat,tt,method )
  S = estimates[[1]];h = estimates[[2]]; C = estimates[[3]];rm(estimates)
  result = est_gain_RMST(dat,tt,S,h,C,margg,g,t0)
  return(result)
}

# 70% survival at time 21000
result_age <- sapply(c("super"),oneW.surv.RMST,t0=21000,which.W="age")
result_gender <- sapply(c("super"),oneW.surv.RMST,t0=21000,which.W="Gender")
result_race <- sapply(c("super"),oneW.surv.RMST,t0=21000,which.W="race")
result_ethnicity <- sapply(c("super"),oneW.surv.RMST,t0=21000,which.W="ethnicity")
result_card <- sapply(c("super"),oneW.surv.RMST,t0=21000,which.W="cardiovascular_disease")
result_hyper <- sapply(c("super"),oneW.surv.RMST,t0=21000,which.W="hypertension")
result_dm <- sapply(c("super"),oneW.surv.RMST,t0=21000,which.W="dm")
result_kidney <- sapply(c("super"),oneW.surv.RMST,t0=21000,which.W="chronic_kidney_disease")
result_chole_meds <- sapply(c("super"),oneW.surv.RMST,t0=21000,which.W="cholesterol_meds")
result_hyper_meds <- sapply(c("super"),oneW.surv.RMST,t0=21000,which.W="hypertension_meds")
result_bmi <- sapply(c("super"),oneW.surv.RMST,t0=21000,which.W="bmi.dis") 



moreW.surv.RMST <-  function(which.W,t0,method){
  if(any(which.W=="bmi.dis")){
    data <- data[!is.na(data$bmi.dis),]
  } else {
    data <- data
  }
  W <<- subset(data,select=which.W)
  
  time = data$minutes_to_discharge_or_censor
  event = data$event
  tt = sort(unique(time[event ==1 & time<=t0]))
  
  dat = data.frame(cbind(time,event,W,1-event))
  colnames(dat)[ncol(dat)] = "censor"
  set.seed(1)
  g = matrix(rep(1-pexp((tt-1)*0.2,rate=0.00001),nrow(dat)),byrow=TRUE,ncol=length(tt))
  margg <- g[1,]
  estimates = calc_shc_M(dat,tt,method)
  S = estimates[[1]];h = estimates[[2]]; C = estimates[[3]];rm(estimates)
  result = est_gain_RMST(dat,tt,S,h,C,margg,g,t0)
  return(result)
}

which.W = colnames(data)[c(3,5,6,18,19,22,23,24,25,27)]
which.W[2] <- "bmi.dis"
result_all <- sapply("super",moreW.surv.RMST ,t0=21000,which.W=which.W)


