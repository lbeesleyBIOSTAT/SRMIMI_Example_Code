

#################################################
#################################################
### Example Script for Binary MNAR Imputation ###
#################################################
#################################################

### Developed by Dr. Lauren J Beesley, PhD
### Contact: lbeesley@umich.edu
### Last Updated: 10/29/2021

### This script provides example code for imputing binary variables with MNAR missingness as described
### in Beesley, Bondarenko, et al. (2021) in Statistical Methods in Medical Research. 

library(dplyr)
library(mice)
library(MASS)
library(ggplot2)
mUnif = function(Y){return(runif(n=1, min = Y, max = MAX))}
mSample = function(Y){return(sample(x=c(0,1), size = 1, prob = c(1-Y,Y)))}
expit = function(x){exp(x)/(1+exp(x))}
logit = function(x){log(x/(1-x))}

#####################
### Simulate Data ###
#####################
Nobs = 2000
DAT = mvrnorm(n = Nobs, mu = c(0,0,0,0,0), Sigma =  rbind(c(1,0.4, 0.2, 0.3, 0.1),
                                                          c(0.4,1, 0.4, 0.1, 0.2),
                                                          c(0.2,0.4, 1, 0.2, 0.3),
                                                          c(0.3,0.1, 0.2, 1, 0.1),
                                                          c(0.1,0.2, 0.3, 0.1, 1)))
X1 = as.numeric(DAT[,1]>0)
X2 = as.numeric(DAT[,2]>0)
X3 = as.numeric(DAT[,3]>0)
X4 = DAT[,4]
X5 = DAT[,5]

phi = 1.5
rho = 1
prob1 = expit(0+ phi*X2 + phi*X3 + rho*X4 + rho*X5)
prob2 = expit(0+ phi*X1 + phi*X3 + rho*X4 + rho*X5)
prob3 = expit(0+ phi*X1 + phi*X2 + rho*X4 + rho*X5)
S1 = apply(matrix(prob1),1, mSample)
S2 = apply(matrix(prob2),1, mSample)
S3 = apply(matrix(prob3),1, mSample)

DAT_all = data.frame(X1 = X1,X2 = X2, X3 = X3, X4 = X4, X5 = X5, S1 = S1, S2 = S2, S3 = S3) #true data
DAT_cc = DAT_all[DAT_all$S1 == 1 & DAT_all$S2 == 1,] #complete case subjects only
DAT_miss = DAT_all 
DAT_miss[DAT_miss$S1==0,'X1'] = NA
DAT_miss[DAT_miss$S2==0,'X2'] = NA
DAT_miss[DAT_miss$S3==0,'X3'] = NA
DAT_miss$X1 = factor(ifelse(is.na(DAT_miss[,'X1']), NA,ifelse(DAT_miss[,'X1']==1,'Yes', 'No')), levels = c('No', 'Yes'))
DAT_miss$X2 = factor(ifelse(is.na(DAT_miss[,'X2']), NA,ifelse(DAT_miss[,'X2']==1,'Yes', 'No')), levels = c('No', 'Yes'))
DAT_miss$X3 = factor(ifelse(is.na(DAT_miss[,'X3']), NA,ifelse(DAT_miss[,'X3']==1,'Yes', 'No')), levels = c('No', 'Yes'))


########################
### Define Functions ###
########################

### This function implements logistic regression imputation with a fixed offset. 
### A vector of offset variable names in provided as an argument in mice called OFFSET_NAMES
mice.impute.logregOFFSET = function (y, ry, x, wy = NULL, ..., OFFSET_NAMES) {
  if (is.null(wy)) 
    wy <- !ry
  #print(head(x))
  if(sum(colnames(x) %in% OFFSET_NAMES)>=1){ ### method can only handle one offset at a time
    OFFSET = x[,colnames(x) %in% OFFSET_NAMES]
    OFFSET = as.matrix(OFFSET)[,1]
  }else{
    OFFSET = matrix(rep(0,length(y)))
  }
  x <- cbind(1,as.matrix(x[,!(colnames(x) %in% OFFSET_NAMES)]))
  expr <- expression(glm.fit(x = x[ry, , drop = FALSE], y = y[ry], 
                             family = binomial(link = logit), offset = OFFSET[ry]))
  fit <- eval(expr)
  fit.sum <- summary.glm(fit)
  beta <- coef(fit)
  rv <- t(chol(mice:::sym(fit.sum$cov.unscaled)))
  beta.star <- beta + rv %*% rnorm(ncol(rv))
  p <- 1/(1 + exp(-(x[wy, , drop = FALSE] %*% beta.star + OFFSET[wy])))
  vec <- as.numeric(runif(nrow(p)) <= p)
  if (is.factor(y)) {
    vec <- factor(vec, c(0, 1), levels(y))
  }
  return(vec)
}

### This function implements the propensity offset approach for BINARY variables
### S corresponds to the missingness model outcome
### Cov corresponds to the matrix of predictors for missingness model
### index gives the index of the binary covariate in Cov for which we are evaluating the propensity
propensity_function = function(S,Cov,index){
  #fit = glm(S~as.matrix(Cov), family = binomial())
  Cov = data.matrix(Cov)
  fit = glm(S~., data = data.frame(S,Cov),family = binomial())
  param = fit$coefficients
  Cov_1 = Cov_0 = Cov
  Cov_1[,index] = 1
  Cov_0[,index] = 0
  prob_1 = expit(cbind(1,data.matrix(Cov_1)) %*% matrix(param))
  prob_0 = expit(cbind(1,data.matrix(Cov_0)) %*% matrix(param))
  Z = log((1-prob_1)/(1-prob_0))
  return(param[index+1]*S+Z)
}

### This function implements logistic regression imputation where parameters are drawn via stratified bootstrap.
### This method is mainly useful for imputing rare binary variables.  
mice.impute.logreg.STRATboot = function (y, ry, x, wy = NULL, ...){
  if (is.null(wy))
    wy <- !ry
  xobs <- x[ry, , drop = FALSE]
  yobs <- y[ry]
  if(is.factor(y)){
    s1 <- sample(which(yobs==levels(y)[2]), sum(yobs==levels(y)[2]), replace = TRUE)
    s0 <- sample(which(yobs==levels(y)[1]), sum(yobs==levels(y)[1]), replace = TRUE)
  }else{
    s1 <- sample(which(yobs==1), sum(yobs==1), replace = TRUE)
    s0 <- sample(which(yobs==0), sum(yobs==0), replace = TRUE)
  }
  s = c(s1,s0)
  doty <- y
  doty[ry] <- yobs[s]
  dotx <- x
  dotx[ry, ] <- xobs[s, , drop = FALSE]
  x <- dotx
  y <- doty
  x <- cbind(1, as.matrix(x))
  expr <- expression(glm.fit(x = x[ry, , drop = FALSE], y = y[ry],
                             family = binomial(link = logit)))
  fit <- suppressWarnings(eval(expr))
  beta.star <- coef(fit)
  p <- 1/(1 + exp(-(x[wy, ] %*% beta.star)))
  vec <- (runif(nrow(p)) <= p)
  vec[vec] <- 1
  if (is.factor(y)) {
    vec <- factor(vec, c(0, 1), levels(y))
  }
  return(vec)
}


### Custom function for imputing using the exact conditional distribution
### This function can be used as a (crude) template for users to create their own custom mice methods
mice.impute.logregExact = function (y, ry, x, wy = NULL,...) {
  if (is.null(wy))
    wy <- !ry
  dat = data.frame(y=y,x) #x includes S values
  ### Covariate model draws
  dat_cc = dat[ry,]
  inds = sample(c(1:length(dat_cc[,1])), length(dat_cc[,1]), replace = TRUE)
  dat_boot = dat_cc[inds,]
  if(!('X1' %in% colnames(x))){
    fit = glm(y~X2+X3+X4+X5, family = binomial(), data = dat_boot)
    XB = as.matrix(cbind(1,x[,c('X2','X3','X4','X5')])) %*% as.matrix(fit$coefficients)
  }else if(!('X2' %in% colnames(x))){
    fit = glm(y~X1+X3+X4+X5, family = binomial(), data = dat_boot)
    XB = as.matrix(cbind(1,x[,c('X1','X3','X4','X5')])) %*% as.matrix(fit$coefficients)
  }else if(!('X3' %in% colnames(x))){
    fit = glm(y~X1+X2+X4+X5, family = binomial(), data = dat_boot)
    XB =  as.matrix(cbind(1,x[,c('X1','X2','X4','X5')])) %*% as.matrix(fit$coefficients)
  }
  COV_PROB = expit(XB)
  ### Missing Data Models
  inds = sample(c(1:length(dat[,1])), length(dat[,1]), replace = TRUE)
  dat_boot = dat[inds,]
  if(!('X1' %in% colnames(x))){
    fit2 = glm(S2~y+X3+X4+X5, data = dat_boot, family = binomial())
    fit3 = glm(S3~y+X2+X4+X5, data = dat_boot, family = binomial())
    param2 = as.numeric(coef(fit2))
    param3 = as.numeric(coef(fit3))
    PROP2 = expit(param2[1] + param2[2]*1+param2[3]*x[,'X3']+param2[4]*x[,'X4']+ param2[5]*x[,'X5'])
    PROP3 = expit(param3[1] + param3[2]*1+param3[3]*x[,'X2']+param3[4]*x[,'X4']+ param3[5]*x[,'X5'])
    CONST_1 = ifelse(dat$S2==1,PROP2,1-PROP2)*ifelse(dat$S3==1,PROP3,1-PROP3)
    PROP2 = expit(param2[1] + param2[2]*0+param2[3]*x[,'X3']+param2[4]*x[,'X4']+ param2[5]*x[,'X5'])
    PROP3 = expit(param3[1] + param3[2]*0+param3[3]*x[,'X2']+param3[4]*x[,'X4']+ param3[5]*x[,'X5'])
    CONST_0 = ifelse(dat$S2==1,PROP2,1-PROP2)*ifelse(dat$S3==1,PROP3,1-PROP3)
  }else if(!('X2' %in% colnames(x))){
    fit1 = glm(S1~y+X3+X4+X5, data = dat_boot, family = binomial())
    fit3 = glm(S3~X1+y+X4+X5, data = dat_boot, family = binomial())
    param1 = as.numeric(coef(fit1))
    param3 = as.numeric(coef(fit3))
    PROP1 = expit(param1[1] + param1[2]*1+param1[3]*x[,'X3']+param1[4]*x[,'X4']+ param1[5]*x[,'X5'])
    PROP3 = expit(param3[1] + param3[2]*x[,'X1']+param3[3]*1+param3[4]*x[,'X4']+ param3[5]*x[,'X5'])
    CONST_1 = ifelse(dat$S1==1,PROP1,1-PROP1)*ifelse(dat$S3==1,PROP3,1-PROP3)
    PROP1 = expit(param1[1] + param1[2]*0+param1[3]*x[,'X3']+param1[4]*x[,'X4']+ param1[5]*x[,'X5'])
    PROP3 = expit(param3[1] + param3[2]*x[,'X1']+param3[3]*0+param3[4]*x[,'X4']+ param3[5]*x[,'X5'])
    CONST_0 = ifelse(dat$S1==1,PROP1,1-PROP1)*ifelse(dat$S3==1,PROP3,1-PROP3)
  }else if(!('X3' %in% colnames(x))){
    fit1 = glm(S1~X2+y+X4+X5, data = dat_boot, family = binomial())
    fit2 = glm(S2~X1+y+X4+X5, data = dat_boot, family = binomial())
    param1 = as.numeric(coef(fit1))
    param2 = as.numeric(coef(fit2))
    PROP1 = expit(param1[1] + param1[2]*x[,'X2']+param1[3]*1+param1[4]*x[,'X4']+ param1[5]*x[,'X5'])
    PROP2 = expit(param2[1] + param2[2]*x[,'X1']+param2[3]*1+param2[4]*x[,'X4']+ param2[5]*x[,'X5'])
    CONST_1 = ifelse(dat$S1==1,PROP1,1-PROP1)*ifelse(dat$S2==1,PROP2,1-PROP2)
    PROP1 = expit(param1[1] + param1[2]*x[,'X2']+param1[3]*0+param1[4]*x[,'X4']+ param1[5]*x[,'X5'])
    PROP2 = expit(param2[1] + param2[2]*x[,'X1']+param2[3]*0+param2[4]*x[,'X4']+ param2[5]*x[,'X5'])
    CONST_0 = ifelse(dat$S1==1,PROP1,1-PROP1)*ifelse(dat$S2==1,PROP2,1-PROP2)
  }
  EXACT_PROB = (COV_PROB*CONST_1)/(COV_PROB*CONST_1 + (1-COV_PROB)*CONST_0)
  proposal = apply(matrix(EXACT_PROB[wy]),1, mSample)
  return(proposal)
}


##################
##################
### Imputation ###
##################
##################

METHOD_NAMES = c('No Missingness','Complete Case','SRMI','SRMI-MI','SRMI-Offset(Binary)',
                'SRMI-Interactions R','SRMI-Interactions X', 'SRMI-MI Stratified Bootstrap','SRMI-Exact' )
RESULTS = data.frame(METHODS = METHOD_NAMES)
RESULTS$Mean_X1 = NA
RESULTS$Var_X1 = NA

######################
### No Missingness ###
######################
Mean_X1 = mean(DAT_all$X1)
RESULTS[RESULTS$METHODS == 'No Missingness','Mean_X1'] = Mean_X1
RESULTS[RESULTS$METHODS == 'No Missingness','Var_X1'] = Mean_X1*(1-Mean_X1)/Nobs

#####################
### Complete Case ###
#####################
Mean_X1 = mean(DAT_cc$X1)
RESULTS[RESULTS$METHODS == 'Complete Case','Mean_X1'] = Mean_X1
RESULTS[RESULTS$METHODS == 'Complete Case','Var_X1'] = Mean_X1*(1-Mean_X1)/length(DAT_cc$X1)

###########
### MAR ###
###########
impute_MAR <- mice(DAT_miss, m=10, method=c("logreg", "logreg","logreg",NA,NA,NA,NA,NA), printFlag=F, maxit = 0)
pred <- impute_MAR$predictorMatrix
pred[,] = 0
pred["X1",c("X2", "X3", "X4", "X5")] <- rep(1,4)
pred["X2",c("X1", "X3", "X4", "X5")] <- rep(1,4)
pred["X3",c("X1", "X2", "X4", "X5")] <- rep(1,4)
impute_MAR <- mice(DAT_miss, m=10, predictorMatrix=pred, method=c("logreg", "logreg","logreg",NA,NA,NA,NA,NA), printFlag=F, maxit = 50)
Mean_X1 = mean(unlist(with(impute_MAR,mean(as.numeric(as.character(X1)=='Yes')))$analyses))
RESULTS[RESULTS$METHODS == 'SRMI','Mean_X1'] = Mean_X1
means_long = with(impute_MAR,mean(as.numeric(as.character(X1)=='Yes')))$analyses
vars_long = lapply(c(1:10),function(x){return(means_long[[x]]*(1-means_long[[x]])/Nobs)})
Var_X1 = mean(unlist(vars_long)) + (1+(1/10))*var(unlist(means_long))
RESULTS[RESULTS$METHODS == 'SRMI','Var_X1'] = Var_X1


##############################
### Missingness Indicators ###
##############################
pred <- impute_MAR$predictorMatrix
pred[,] = 0
pred["X1",c("X2", "X3", "X4", "X5", "S2", "S3")] <- rep(1,6)
pred["X2",c("X1", "X3", "X4", "X5", "S1", 'S3')] <- rep(1,6)
pred["X3",c("X1", "X2", "X4", "X5", "S1", 'S2')] <- rep(1,6)
impute_MNAR <- mice(DAT_miss, m=10, predictorMatrix=pred, method=c("logreg", "logreg","logreg",NA,NA,NA,NA,NA), printFlag=F, maxit = 50)
Mean_X1 = mean(unlist(with(impute_MNAR,mean(as.numeric(as.character(X1)=='Yes')))$analyses))
RESULTS[RESULTS$METHODS == 'SRMI-MI','Mean_X1'] = Mean_X1
means_long = with(impute_MNAR,mean(as.numeric(as.character(X1)=='Yes')))$analyses
vars_long = lapply(c(1:10),function(x){return(means_long[[x]]*(1-means_long[[x]])/Nobs)})
Var_X1 = mean(unlist(vars_long)) + (1+(1/10))*var(unlist(means_long))
RESULTS[RESULTS$METHODS == 'SRMI-MI','Var_X1'] = Var_X1

#################################
### Propensity Offsets Method ###
#################################
### Define offset functions via passive imputation
expr1_X2 <- expression(propensity_function(S=S1,Cov = data.frame(as.numeric(as.character(X2)=='Yes'),as.numeric(as.character(X3)=='Yes'),X4,X5), index = 1))
expr1_X3 <- expression(propensity_function(S=S1,Cov = data.frame(as.numeric(as.character(X2)=='Yes'),as.numeric(as.character(X3)=='Yes'),X4,X5), index = 2))
expr2_X1 <- expression(propensity_function(S=S2,Cov = data.frame(as.numeric(as.character(X1)=='Yes'),as.numeric(as.character(X3)=='Yes'),X4,X5), index = 1))
expr2_X3 <- expression(propensity_function(S=S2,Cov = data.frame(as.numeric(as.character(X1)=='Yes'),as.numeric(as.character(X3)=='Yes'),X4,X5), index = 2))
expr3_X1 <- expression(propensity_function(S=S3,Cov = data.frame(as.numeric(as.character(X1)=='Yes'),as.numeric(as.character(X2)=='Yes'),X4,X5), index = 1))
expr3_X2 <- expression(propensity_function(S=S3,Cov = data.frame(as.numeric(as.character(X1)=='Yes'),as.numeric(as.character(X2)=='Yes'),X4,X5), index = 2))
OFFSET_NAMES = c('offset_X1', 'offset_X2', 'offset_X3')
METHODS = c("logregOFFSET","logregOFFSET","logregOFFSET",NA,NA,NA,NA,NA,
            paste("~I(",expr2_X1,"+",expr3_X1,")"),
            paste("~I(",expr1_X2,"+",expr3_X2,")"),
            paste("~I(",expr1_X3,"+",expr2_X3,")"))
### Initialize missing values
DAT_miss_aug = data.frame(DAT_miss, 
    offset_X1 = propensity_function(S=DAT_miss$S2,Cov = data.frame(as.numeric(as.character(DAT_miss$X1)=='Yes'),as.numeric(as.character(DAT_miss$X3)=='Yes'),DAT_miss$X4,DAT_miss$X5), index = 1)+
                propensity_function(S=DAT_miss$S3,Cov = data.frame(as.numeric(as.character(DAT_miss$X1)=='Yes'),as.numeric(as.character(DAT_miss$X2)=='Yes'),DAT_miss$X4,DAT_miss$X5), index = 1),
    offset_X2 = propensity_function(S=DAT_miss$S1,Cov = data.frame(as.numeric(as.character(DAT_miss$X2)=='Yes'),as.numeric(as.character(DAT_miss$X3)=='Yes'),DAT_miss$X4,DAT_miss$X5), index = 1)+
                propensity_function(S=DAT_miss$S3,Cov = data.frame(as.numeric(as.character(DAT_miss$X1)=='Yes'),as.numeric(as.character(DAT_miss$X2)=='Yes'),DAT_miss$X4,DAT_miss$X5), index = 2),
    offset_X3 = propensity_function(S=DAT_miss$S1,Cov = data.frame(as.numeric(as.character(DAT_miss$X2)=='Yes'),as.numeric(as.character(DAT_miss$X3)=='Yes'),DAT_miss$X4,DAT_miss$X5), index = 2)+
                propensity_function(S=DAT_miss$S2,Cov = data.frame(as.numeric(as.character(DAT_miss$X1)=='Yes'),as.numeric(as.character(DAT_miss$X3)=='Yes'),DAT_miss$X4,DAT_miss$X5), index = 2)
)
### Perform imputation
impute_MNAR_PROP <- mice(DAT_miss_aug, m=10, method=METHODS, printFlag=F, maxit = 0, OFFSET_NAMES=OFFSET_NAMES)
pred <- impute_MNAR_PROP$predictorMatrix
pred[,] = 0
pred["X1",c("X2", "X3", "X4", "X5", "offset_X1")] <- rep(1,5)
pred["X2",c("X1", "X3", "X4", "X5", "offset_X2")] <- rep(1,5)
pred["X3",c("X1", "X2", "X4", "X5", "offset_X3")] <- rep(1,5)
impute_MNAR_PROP <- mice(DAT_miss_aug, m=10,method=METHODS, predictorMatrix = pred, printFlag=F, maxit = 50, OFFSET_NAMES=OFFSET_NAMES)
Mean_X1 = mean(unlist(with(impute_MNAR_PROP,mean(as.numeric(as.character(X1)=='Yes')))$analyses))
RESULTS[RESULTS$METHODS == 'SRMI-Offset(Binary)','Mean_X1'] = Mean_X1
means_long = with(impute_MNAR_PROP,mean(as.numeric(as.character(X1)=='Yes')))$analyses
vars_long = lapply(c(1:10),function(x){return(means_long[[x]]*(1-means_long[[x]])/Nobs)})
Var_X1 = mean(unlist(vars_long)) + (1+(1/10))*var(unlist(means_long))
RESULTS[RESULTS$METHODS == 'SRMI-Offset(Binary)','Var_X1'] = Var_X1





###########################
### Interactions with R ###
###########################
### Define passive imputation for interactions
METHODS = c("logreg","logreg","logreg",NA,NA,NA,NA,NA,
            paste("~I(as.numeric(as.character(X1) == 'Yes')*S1)"),
            paste("~I(as.numeric(as.character(X1) == 'Yes')*S2)"),
            paste("~I(as.numeric(as.character(X1) == 'Yes')*S3)"),
            paste("~I(as.numeric(as.character(X2) == 'Yes')*S1)"),
            paste("~I(as.numeric(as.character(X2) == 'Yes')*S2)"),
            paste("~I(as.numeric(as.character(X2) == 'Yes')*S3)"),
            paste("~I(as.numeric(as.character(X3) == 'Yes')*S1)"),
            paste("~I(as.numeric(as.character(X3) == 'Yes')*S2)"),
            paste("~I(as.numeric(as.character(X3) == 'Yes')*S3)"),
            paste("~I(X4*S1)"),paste("~I(X4*S2)"),paste("~I(X4*S3)"),
            paste("~I(X5*S1)"),paste("~I(X5*S2)"),paste("~I(X5*S3)"))
### Initialize interactions variables. First index is X variable, second index is missingness indicator.
DAT_miss_aug = data.frame(DAT_miss, 
                          I_11 = as.numeric( as.numeric(as.character(DAT_miss[,'X1'])=='Yes')*DAT_miss[,'S1']), 
                          I_12 = as.numeric( as.numeric(as.character(DAT_miss[,'X1'])=='Yes')*DAT_miss[,'S2']), 
                          I_13 = as.numeric( as.numeric(as.character(DAT_miss[,'X1'])=='Yes')*DAT_miss[,'S3']), 
                          I_21 = as.numeric( as.numeric(as.character(DAT_miss[,'X2'])=='Yes')*DAT_miss[,'S1']), 
                          I_22 = as.numeric( as.numeric(as.character(DAT_miss[,'X2'])=='Yes')*DAT_miss[,'S2']), 
                          I_23 = as.numeric( as.numeric(as.character(DAT_miss[,'X2'])=='Yes')*DAT_miss[,'S3']), 
                          I_31 = as.numeric( as.numeric(as.character(DAT_miss[,'X3'])=='Yes')*DAT_miss[,'S1']), 
                          I_32 = as.numeric( as.numeric(as.character(DAT_miss[,'X3'])=='Yes')*DAT_miss[,'S2']), 
                          I_33 = as.numeric( as.numeric(as.character(DAT_miss[,'X3'])=='Yes')*DAT_miss[,'S3']), 
                          I_41 = as.numeric( DAT_miss[,'X4']*DAT_miss[,'S1']),
                          I_42 = as.numeric( DAT_miss[,'X4']*DAT_miss[,'S2']), 
                          I_43 = as.numeric( DAT_miss[,'X4']*DAT_miss[,'S3']), 
                          I_51 = as.numeric( DAT_miss[,'X5']*DAT_miss[,'S1']), 
                          I_52 = as.numeric( DAT_miss[,'X5']*DAT_miss[,'S2']),
                          I_53 = as.numeric( DAT_miss[,'X5']*DAT_miss[,'S3']))
### Perform imputation
impute_MNAR_INTR <- mice(DAT_miss_aug, m=10, method=METHODS, printFlag=F, maxit = 0)
pred <- impute_MNAR_INTR$predictorMatrix
pred[,] = 0
pred["X1",c("X2", "X3", "X4", "X5",  "S2", "S3",   "I_32", "I_42", "I_52",   "I_23", "I_43", "I_53")] <- rep(1,12)
pred["X2",c("X1", "X3", "X4", "X5",  "S1", "S3",   "I_31", "I_41", "I_51",   "I_13", "I_43", "I_53")] <- rep(1,12)
pred["X3",c("X1", "X2", "X4", "X5",  "S1", "S2",   "I_21", "I_41", "I_51",   "I_12", "I_42", "I_52")] <- rep(1,12)
impute_MNAR_INTR <- mice(DAT_miss_aug, m=10,method=METHODS, predictorMatrix = pred, printFlag=F, maxit = 50)
Mean_X1 = mean(unlist(with(impute_MNAR_INTR,mean(as.numeric(as.character(X1)=='Yes')))$analyses))
RESULTS[RESULTS$METHODS == 'SRMI-Interactions R','Mean_X1'] = Mean_X1
means_long = with(impute_MNAR_INTR,mean(as.numeric(as.character(X1)=='Yes')))$analyses
vars_long = lapply(c(1:10),function(x){return(means_long[[x]]*(1-means_long[[x]])/Nobs)})
Var_X1 = mean(unlist(vars_long)) + (1+(1/10))*var(unlist(means_long))
RESULTS[RESULTS$METHODS == 'SRMI-Interactions R','Var_X1'] = Var_X1




###########################
### Interactions with X ###
###########################
### Define passive imputation for interactions
METHODS = c("logreg","logreg","logreg",NA,NA,NA,NA,NA,
            paste("~I(as.numeric(as.character(X1)=='Yes')*as.numeric(as.character(X1)=='Yes'))"),
            paste("~I(as.numeric(as.character(X1)=='Yes')*as.numeric(as.character(X2)=='Yes'))"),
            paste("~I(as.numeric(as.character(X1)=='Yes')*as.numeric(as.character(X3)=='Yes'))"),
            paste("~I(as.numeric(as.character(X1)=='Yes')*X4)"),
            paste("~I(as.numeric(as.character(X1)=='Yes')*X5)"),
            paste("~I(as.numeric(as.character(X2)=='Yes')*as.numeric(as.character(X2)=='Yes'))"),
            paste("~I(as.numeric(as.character(X2)=='Yes')*as.numeric(as.character(X3)=='Yes'))"),
            paste("~I(as.numeric(as.character(X2)=='Yes')*X4)"),
            paste("~I(as.numeric(as.character(X2)=='Yes')*X5)"),
            paste("~I(as.numeric(as.character(X3)=='Yes')*as.numeric(as.character(X3)=='Yes'))"),
            paste("~I(as.numeric(as.character(X3)=='Yes')*X4)"),
            paste("~I(as.numeric(as.character(X3)=='Yes')*X5)"),
            paste("~I(X4*X4)"),paste("~I(X4*X5)"),paste("~I(X5*X5)"))
### Initialize interactions variables. Indices correspond to variables in X
DAT_miss_aug = data.frame(DAT_miss, 
                          I_11 = as.numeric(as.numeric(as.character(DAT_miss[,'X1'])=='Yes')*as.numeric(as.character(DAT_miss[,'X1'])=='Yes')), 
                          I_21 = as.numeric(as.numeric(as.character(DAT_miss[,'X2'])=='Yes')*as.numeric(as.character(DAT_miss[,'X1'])=='Yes')), 
                          I_31 = as.numeric(as.numeric(as.character(DAT_miss[,'X3'])=='Yes')*as.numeric(as.character(DAT_miss[,'X1'])=='Yes')), 
                          I_41 = as.numeric(DAT_miss[,'X4']*as.numeric(as.character(DAT_miss[,'X1'])=='Yes')), 
                          I_51 = as.numeric(DAT_miss[,'X5']*as.numeric(as.character(DAT_miss[,'X1'])=='Yes')), 
                          I_22 = as.numeric(as.numeric(as.character(DAT_miss[,'X2'])=='Yes')*as.numeric(as.character(DAT_miss[,'X2'])=='Yes')), 
                          I_32 = as.numeric(as.numeric(as.character(DAT_miss[,'X3'])=='Yes')*as.numeric(as.character(DAT_miss[,'X2'])=='Yes')), 
                          I_42 = as.numeric(DAT_miss[,'X4']*as.numeric(as.character(DAT_miss[,'X2'])=='Yes')), 
                          I_52 = as.numeric(DAT_miss[,'X5']*as.numeric(as.character(DAT_miss[,'X2'])=='Yes')),
                          I_33 = as.numeric(as.numeric(as.character(DAT_miss[,'X3'])=='Yes')*as.numeric(as.character(DAT_miss[,'X3'])=='Yes')), 
                          I_43 = as.numeric(DAT_miss[,'X4']*as.numeric(as.character(DAT_miss[,'X3'])=='Yes')), 
                          I_53 = as.numeric(DAT_miss[,'X5']*as.numeric(as.character(DAT_miss[,'X3'])=='Yes')),
                          I_44 = as.numeric(DAT_miss[,'X4']*DAT_miss[,'X4']),
                          I_54 = as.numeric(DAT_miss[,'X5']*DAT_miss[,'X4']),
                          I_55 = as.numeric(DAT_miss[,'X5']*DAT_miss[,'X5']))
### Perform imputation
impute_MNAR_INTX <- mice(DAT_miss_aug, m=10, method=METHODS, printFlag=F, maxit = 1)
pred <- impute_MNAR_INTX$predictorMatrix
pred[,] = 0
pred["X1",c("X2", "X3", "X4", "X5",    "S2", "S3",   "I_32","I_42","I_52",    "I_43","I_53","I_44", "I_54", "I_55")] <- rep(1,14)
pred["X2",c("X1", "X3", "X4", "X5",    "S1", "S3",   "I_31","I_41","I_51",    "I_43","I_53","I_44", "I_54", "I_55")] <- rep(1,14)
pred["X3",c("X1", "X2", "X4", "X5",    "S1", "S2",   "I_21","I_41","I_51",    "I_42","I_52","I_44", "I_54", "I_55")] <- rep(1,14)
impute_MNAR_INTX <- mice(DAT_miss_aug, m=10,method=METHODS, predictorMatrix = pred, printFlag=F, maxit = 50)
Mean_X1 = mean(unlist(with(impute_MNAR_INTX,mean(as.numeric(as.character(X1)=='Yes')))$analyses))
RESULTS[RESULTS$METHODS == 'SRMI-Interactions X','Mean_X1'] = Mean_X1
means_long = with(impute_MNAR_INTX,mean(as.numeric(as.character(X1)=='Yes')))$analyses
vars_long = lapply(c(1:10),function(x){return(means_long[[x]]*(1-means_long[[x]])/Nobs)})
Var_X1 = mean(unlist(vars_long)) + (1+(1/10))*var(unlist(means_long))
RESULTS[RESULTS$METHODS == 'SRMI-Interactions X','Var_X1'] = Var_X1


#####################################################
### Missingness Indicators + Stratified Bootstrap ###
#####################################################
pred <- impute_MAR$predictorMatrix
pred[,] = 0
pred["X1",c("X2", "X3", "X4", "X5", "S2", "S3")] <- rep(1,6)
pred["X2",c("X1", "X3", "X4", "X5", "S1", 'S3')] <- rep(1,6)
pred["X3",c("X1", "X2", "X4", "X5", "S1", 'S2')] <- rep(1,6)
impute_STRAT <- mice(DAT_miss, m=10, predictorMatrix=pred, method=c("logreg.STRATboot", "logreg.STRATboot","logreg.STRATboot",NA,NA,NA,NA,NA), printFlag=F, maxit = 50)
Mean_X1 = mean(unlist(with(impute_STRAT,mean(as.numeric(as.character(X1)=='Yes')))$analyses))
RESULTS[RESULTS$METHODS == 'SRMI-MI Stratified Bootstrap','Mean_X1'] = Mean_X1
means_long = with(impute_STRAT,mean(as.numeric(as.character(X1)=='Yes')))$analyses
vars_long = lapply(c(1:10),function(x){return(means_long[[x]]*(1-means_long[[x]])/Nobs)})
Var_X1 = mean(unlist(vars_long)) + (1+(1/10))*var(unlist(means_long))
RESULTS[RESULTS$METHODS == 'SRMI-MI Stratified Bootstrap','Var_X1'] = Var_X1


#####################################
### Exact Imputation Distribution ###
#####################################
pred <- impute_MAR$predictorMatrix
pred[,] = 0
pred["X1",c("X2", "X3", "X4", "X5", "S2", "S3")] <- rep(1,6)
pred["X2",c("X1", "X3", "X4", "X5", "S1", 'S3')] <- rep(1,6)
pred["X3",c("X1", "X2", "X4", "X5", "S1", 'S2')] <- rep(1,6)
DAT_miss_temp = DAT_miss
DAT_miss_temp$X1 = as.numeric(as.character(DAT_miss_temp$X1)=='Yes')
DAT_miss_temp$X2 = as.numeric(as.character(DAT_miss_temp$X2)=='Yes')
DAT_miss_temp$X3 = as.numeric(as.character(DAT_miss_temp$X3)=='Yes')
impute_EXACT <- mice(DAT_miss_temp, m=10, predictorMatrix=pred, method=c("logregExact", "logregExact","logregExact",NA,NA,NA,NA,NA), printFlag=F, maxit = 50)
Mean_X1 = mean(unlist(with(impute_EXACT,mean(X1))$analyses))
RESULTS[RESULTS$METHODS == 'SRMI-Exact','Mean_X1'] = Mean_X1
means_long = with(impute_EXACT,mean(X1))$analyses
vars_long = lapply(c(1:10),function(x){return(means_long[[x]]*(1-means_long[[x]])/Nobs)})
Var_X1 = mean(unlist(vars_long)) + (1+(1/10))*var(unlist(means_long))
RESULTS[RESULTS$METHODS == 'SRMI-Exact','Var_X1'] = Var_X1


#######################
### Compare Results ###
#######################
### Note: these are results for a single simulated dataset. The true mean of random variable X1 is 0.50 (dashed black line). 
### The mean of the simulated data sample without missingness is shown by the gray dashed line.
ggplot(RESULTS, aes(x = Mean_X1, y = METHODS, xmin = Mean_X1-1.96*sqrt(Var_X1), xmax = Mean_X1+1.96*sqrt(Var_X1))) +
  geom_pointrange(aes(fill = METHODS), position = ggstance::position_dodgev(height=1), shape = 21, size = 1.5) +
  geom_vline(xintercept=0.50,linetype=2)+
  scale_y_discrete(limits = rev(METHOD_NAMES))+
  xlab("Mean of X1") + ylab('Method')+labs(title='Estimated Mean of X1 and 95% Confidence Interval')+
  theme_classic() +theme(legend.position ='none')+
  geom_vline(aes(xintercept = Mean_X1), color = 'gray', linetype = 2, data = RESULTS[RESULTS$METHODS =='No Missingness', ])




