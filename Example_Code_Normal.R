

#################################################
#################################################
### Example Script for Normal MNAR Imputation ###
#################################################
#################################################

### Developed by Dr. Lauren J Beesley, PhD
### Contact: lbeesley@umich.edu
### Last Updated: 10/29/2021

### This script provides example code for imputing normally-distributed variables with MNAR missingness as described
### in Beesley, Bondarenko, et al. (2021) in Statistical Methods in Medical Research. 

library(dplyr)
library(mice)
library(MASS)
library(ggplot2)
library(miceadds)
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
X1 = DAT[,1]
X2 = DAT[,2]
X3 = DAT[,3]
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

########################
### Define Functions ###
########################

### This function implements linear regression imputation with a fixed offset. 
### A vector of offset variable names in provided as an argument in mice called OFFSET_NAMES
mice.impute.normOFFSET = function (y, ry, x, wy = NULL, ..., OFFSET_NAMES) {
  if (is.null(wy)) 
    wy <- !ry
  #print(colnames(x))
  if(sum(colnames(x) %in% OFFSET_NAMES)>=1){ ### method can only handle one offset at a time
    OFFSET = x[,colnames(x) %in% OFFSET_NAMES]
    OFFSET = as.matrix(OFFSET)[,1]
  }else{
    OFFSET = matrix(rep(0,length(y)))
  }
  x <- as.matrix(x[,!(colnames(x) %in% OFFSET_NAMES)])
  fit = glm(y~x, family = 'gaussian', offset = OFFSET)
  sigma.star <- sqrt(sum((residuals(fit))^2)/rchisq(1, Nobs-dim(x)[2]-1))
  beta.star <- MASS::mvrnorm(n=1, mu = as.numeric(coef(fit)), Sigma = (sigma.star^2)*summary(fit)$cov.unscaled)
  vec = OFFSET[wy] + cbind(1,x[wy, ]) %*% matrix(beta.star) + rnorm(sum(wy)) * sigma.star
  return(vec)
}

### This function implements the propensity offset approach for BINARY variables
### S corresponds to the missingness model outcome
### Cov corresponds to the matrix of predictors for missingness model
### index gives the index of the binary covariate in Cov for which we are evaluating the propensity
propensity_function = function(S,Cov,index = 1){
  Cov = data.matrix(Cov)
  fit = glm(S~., data = data.frame(S,Cov),family = binomial())
  predictions = predict(fit, newdata = data.frame(S,Cov),type = 'response')
  Z=coef(fit)[index+1]*(S-predictions)
  return(Z)
}

### Custom function for imputing using the exact conditional distribution
### This function can be used as a template for users to create their own custom mice methods
mice.impute.normExact = function (y, ry, x, wy = NULL,...) {
  if (is.null(wy))
    wy <- !ry
  dat = data.frame(y=y,x) #x includes S values
  ### Covariate model draws
  dat_cc = dat[ry,]
  inds = sample(c(1:length(dat_cc[,1])), length(dat_cc[,1]), replace = TRUE)
  dat_boot = dat_cc[inds,]
  if(!('X1' %in% colnames(x))){
    fit = glm(y~X2+X3+X4+X5, family = gaussian(), data = dat_boot)
    XB = as.matrix(cbind(1,x[,c('X2','X3','X4','X5')])) %*% as.matrix(fit$coefficients)
  }else if(!('X2' %in% colnames(x))){
    fit = glm(y~X1+X3+X4+X5, family = gaussian(), data = dat_boot)
    XB = as.matrix(cbind(1,x[,c('X1','X3','X4','X5')])) %*% as.matrix(fit$coefficients)
  }else if(!('X3' %in% colnames(x))){
    fit = glm(y~X1+X2+X4+X5, family = gaussian(), data = dat_boot)
    XB =  as.matrix(cbind(1,x[,c('X1','X2','X4','X5')])) %*% as.matrix(fit$coefficients)
  }
  ### Missing Data Models
  MISSING = as.numeric(wy)
  CONST = 1
  iterlim = 1
  proposal = y
  inds = sample(c(1:length(dat[,1])), length(dat[,1]), replace = TRUE)
  dat_boot = dat[inds,]
  if(!('X1' %in% colnames(x))){
    fit2 = glm(S2~y+X3+X4+X5, data = dat_boot, family = binomial())
    fit3 = glm(S3~y+X2+X4+X5, data = dat_boot, family = binomial())
    param2 = as.numeric(coef(fit2))
    param3 = as.numeric(coef(fit3))
    while(sum(MISSING) != 0 & iterlim != 50){
      draw = rnorm(n=sum(MISSING), mean = XB[MISSING == 1], sd = sqrt(summary(fit)$dispersion))
      U = runif(n=sum(MISSING), min=0,max=1)
      proposal[MISSING == 1] = draw
      PROP2 = expit(param2[1] + param2[2]*proposal+param2[3]*x[,'X3']+param2[4]*x[,'X4']+ param2[5]*x[,'X5'])
      PROP3 = expit(param3[1] + param3[2]*proposal+param3[3]*x[,'X2']+param3[4]*x[,'X4']+ param3[5]*x[,'X5'])
      CONST = ifelse(dat$S2==1,PROP2,1-PROP2)*ifelse(dat$S3==1,PROP3,1-PROP3)
      MISSING[MISSING == 1] = as.numeric(U>CONST[MISSING==1])
      iterlim = iterlim + 1
    }
  }else if(!('X2' %in% colnames(x))){
    fit1 = glm(S1~y+X3+X4+X5, data = dat_boot, family = binomial())
    fit3 = glm(S3~X1+y+X4+X5, data = dat_boot, family = binomial())
    param1 = as.numeric(coef(fit1))
    param3 = as.numeric(coef(fit3))
    while(sum(MISSING) != 0 & iterlim != 50){
      draw = rnorm(n=sum(MISSING), mean = XB[MISSING == 1], sd = sqrt(summary(fit)$dispersion))
      U = runif(n=sum(MISSING), min=0,max=1)
      proposal[MISSING == 1] = draw
      PROP1 = expit(param1[1] + param1[2]*proposal+param1[3]*x[,'X3']+param1[4]*x[,'X4']+ param1[5]*x[,'X5'])
      PROP3 = expit(param3[1] + param3[2]*x[,'X1']+param3[3]*proposal+param3[4]*x[,'X4']+ param3[5]*x[,'X5'])
      CONST = ifelse(dat$S1==1,PROP1,1-PROP1)*ifelse(dat$S3==1,PROP3,1-PROP3)
      MISSING[MISSING == 1] = as.numeric(U>CONST[MISSING==1])
      iterlim = iterlim + 1
    }
  }else if(!('X3' %in% colnames(x))){
    fit1 = glm(S1~X2+y+X4+X5, data = dat_boot, family = binomial())
    fit2 = glm(S2~X1+y+X4+X5, data = dat_boot, family = binomial())
    param1 = as.numeric(coef(fit1))
    param2 = as.numeric(coef(fit2))
    while(sum(MISSING) != 0 & iterlim != 50){
      draw = rnorm(n=sum(MISSING), mean = XB[MISSING == 1], sd = sqrt(summary(fit)$dispersion))
      U = runif(n=sum(MISSING), min=0,max=1)
      proposal[MISSING == 1] = draw
      PROP1 = expit(param1[1] + param1[2]*x[,'X2']+param1[3]*proposal+param1[4]*x[,'X4']+ param1[5]*x[,'X5'])
      PROP2 = expit(param2[1] + param2[2]*x[,'X1']+param2[3]*proposal+param2[4]*x[,'X4']+ param2[5]*x[,'X5'])
      CONST = ifelse(dat$S1==1,PROP1,1-PROP1)*ifelse(dat$S2==1,PROP2,1-PROP2)
      MISSING[MISSING == 1] = as.numeric(U>CONST[MISSING==1])
      iterlim = iterlim + 1
    }
  }
  ### Propose new imputations
  return(proposal[wy])
}

##################
##################
### Imputation ###
##################
##################

METHOD_NAMES = c('No Missingness','Complete Case','SRMI','SRMI-MI','SRMI-Offset(Normal)',
                 'SRMI-Interactions R','SRMI-TriCube', 'SRMI-Exact' )
RESULTS = data.frame(METHODS = METHOD_NAMES)
RESULTS$Mean_X1 = NA
RESULTS$Var_X1 = NA


######################
### No Missingness ###
######################
RESULTS[RESULTS$METHODS == 'No Missingness','Mean_X1'] = mean(DAT_all$X1)
RESULTS[RESULTS$METHODS == 'No Missingness','Var_X1'] = var(DAT_all$X1)/Nobs

#####################
### Complete Case ###
#####################
RESULTS[RESULTS$METHODS == 'Complete Case','Mean_X1'] = mean(DAT_cc$X1)
RESULTS[RESULTS$METHODS == 'Complete Case','Var_X1'] = var(DAT_cc$X1)/length(DAT_cc$X1)

###########
### MAR ###
###########
impute_MAR <- mice(DAT_miss, m=10, method=c("norm", "norm","norm",NA,NA,NA,NA,NA), printFlag=F, maxit = 0)
pred <- impute_MAR$predictorMatrix
pred[,] = 0
pred["X1",c("X2", "X3", "X4", "X5")] <- rep(1,4)
pred["X2",c("X1", "X3", "X4", "X5")] <- rep(1,4)
pred["X3",c("X1", "X2", "X4", "X5")] <- rep(1,4)
impute_MAR <- mice(DAT_miss, m=10, predictorMatrix=pred, method=c("norm", "norm","norm",NA,NA,NA,NA,NA), printFlag=F, maxit = 50)
Mean_X1 = mean(unlist(with(impute_MAR,mean(X1))$analyses))
RESULTS[RESULTS$METHODS == 'SRMI','Mean_X1'] = Mean_X1
means_long = with(impute_MAR,mean(as.numeric(as.character(X1)=='Yes')))$analyses
vars_long = with(impute_MAR,var(X1)/Nobs)$analyses
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
impute_MNAR <- mice(DAT_miss, m=10, predictorMatrix=pred, method=c("norm", "norm","norm",NA,NA,NA,NA,NA), printFlag=F, maxit = 50)
Mean_X1 = mean(unlist(with(impute_MNAR,mean(X1))$analyses))
RESULTS[RESULTS$METHODS == 'SRMI-MI','Mean_X1'] = Mean_X1
means_long = with(impute_MNAR,mean(as.numeric(as.character(X1)=='Yes')))$analyses
vars_long = with(impute_MNAR,var(X1)/Nobs)$analyses
Var_X1 = mean(unlist(vars_long)) + (1+(1/10))*var(unlist(means_long))
RESULTS[RESULTS$METHODS == 'SRMI-MI','Var_X1'] = Var_X1

#################################
### Propensity Offsets Method ###
#################################
### Define offset functions via passive imputation
sigma2_cc1 = summary(glm(X1~X2+X3+X4+X5, family = gaussian(), data = DAT_cc))$dispersion
sigma2_cc2 = summary(glm(X2~X1+X3+X4+X5, family = gaussian(), data = DAT_cc))$dispersion
sigma2_cc3 = summary(glm(X3~X1+X2+X4+X5, family = gaussian(), data = DAT_cc))$dispersion
expr1_X2 <- expression(sigma2_cc2*propensity_function(S=S1,Cov = data.frame(X2,X3,X4,X5), index = 1))
expr1_X3 <- expression(sigma2_cc3*propensity_function(S=S1,Cov = data.frame(X2,X3,X4,X5), index = 2))
expr2_X1 <- expression(sigma2_cc1*propensity_function(S=S2,Cov = data.frame(X1,X3,X4,X5), index = 1))
expr2_X3 <- expression(sigma2_cc3*propensity_function(S=S2,Cov = data.frame(X1,X3,X4,X5), index = 2))
expr3_X1 <- expression(sigma2_cc1*propensity_function(S=S3,Cov = data.frame(X1,X2,X4,X5), index = 1))
expr3_X2 <- expression(sigma2_cc2*propensity_function(S=S3,Cov = data.frame(X1,X2,X4,X5), index = 2))
OFFSET_NAMES = c('offset_X1', 'offset_X2', 'offset_X3')
METHODS = c("normOFFSET","normOFFSET","normOFFSET",NA,NA,NA,NA,NA,
            paste("~I(",expr2_X1,"+",expr3_X1,")"),
            paste("~I(",expr1_X2,"+",expr3_X2,")"),
            paste("~I(",expr1_X3,"+",expr2_X3,")"))
### Initialize missing values
DAT_miss_aug = data.frame(DAT_miss, 
    offset_X1 = sigma2_cc1*propensity_function(S=DAT_miss$S2,Cov = data.frame(DAT_miss$X1,DAT_miss$X3,DAT_miss$X4,DAT_miss$X5), index = 1)+
                sigma2_cc1*propensity_function(S=DAT_miss$S3,Cov = data.frame(DAT_miss$X1,DAT_miss$X2,DAT_miss$X4,DAT_miss$X5), index = 1),
    offset_X2 = sigma2_cc2*propensity_function(S=DAT_miss$S1,Cov = data.frame(DAT_miss$X2,DAT_miss$X3,DAT_miss$X4,DAT_miss$X5), index = 1)+
                sigma2_cc2*propensity_function(S=DAT_miss$S3,Cov = data.frame(DAT_miss$X1,DAT_miss$X2,DAT_miss$X4,DAT_miss$X5), index = 2),
    offset_X3 = sigma2_cc3*propensity_function(S=DAT_miss$S1,Cov = data.frame(DAT_miss$X2,DAT_miss$X3,DAT_miss$X4,DAT_miss$X5), index = 2)+
                sigma2_cc3*propensity_function(S=DAT_miss$S2,Cov = data.frame(DAT_miss$X1,DAT_miss$X3,DAT_miss$X4,DAT_miss$X5), index = 2)
)
### Perform imputation
impute_MNAR_PROP <- mice(DAT_miss_aug, m=10, method=METHODS, printFlag=F, maxit = 0, OFFSET_NAMES=OFFSET_NAMES)
pred <- impute_MNAR_PROP$predictorMatrix
pred[,] = 0
pred["X1",c("X2", "X3", "X4", "X5", "offset_X1")] <- rep(1,5)
pred["X2",c("X1", "X3", "X4", "X5", "offset_X2")] <- rep(1,5)
pred["X3",c("X1", "X2", "X4", "X5", "offset_X3")] <- rep(1,5)
impute_MNAR_PROP <- mice(DAT_miss_aug, m=10,method=METHODS, predictorMatrix = pred, printFlag=F, maxit = 50, OFFSET_NAMES=OFFSET_NAMES)
Mean_X1 = mean(unlist(with(impute_MNAR_PROP,mean(X1))$analyses))
RESULTS[RESULTS$METHODS == 'SRMI-Offset(Normal)','Mean_X1'] = Mean_X1
means_long = with(impute_MNAR_PROP,mean(as.numeric(as.character(X1)=='Yes')))$analyses
vars_long = with(impute_MNAR_PROP,var(X1)/Nobs)$analyses
Var_X1 = mean(unlist(vars_long)) + (1+(1/10))*var(unlist(means_long))
RESULTS[RESULTS$METHODS == 'SRMI-Offset(Normal)','Var_X1'] = Var_X1





###########################
### Interactions with R ###
###########################
### Define passive imputation for interactions
METHODS = c("norm","norm","norm",NA,NA,NA,NA,NA,
            paste("~I(X1*S1)"),paste("~I(X1*S2)"),paste("~I(X1*S3)"),
            paste("~I(X2*S1)"),paste("~I(X2*S2)"),paste("~I(X2*S3)"),
            paste("~I(X3*S1)"),paste("~I(X3*S2)"),paste("~I(X3*S3)"),
            paste("~I(X4*S1)"),paste("~I(X4*S2)"),paste("~I(X4*S3)"),
            paste("~I(X5*S1)"),paste("~I(X5*S2)"),paste("~I(X5*S3)"))
### Initialize interactions variables. First index is X variable, second index is missingness indicator.
DAT_miss_aug = data.frame(DAT_miss, 
                          I_11 = as.numeric( DAT_miss[,'X1']*DAT_miss[,'S1']), 
                          I_12 = as.numeric( DAT_miss[,'X1']*DAT_miss[,'S2']), 
                          I_13 = as.numeric( DAT_miss[,'X1']*DAT_miss[,'S3']), 
                          I_21 = as.numeric( DAT_miss[,'X2']*DAT_miss[,'S1']), 
                          I_22 = as.numeric( DAT_miss[,'X2']*DAT_miss[,'S2']), 
                          I_23 = as.numeric( DAT_miss[,'X2']*DAT_miss[,'S3']), 
                          I_31 = as.numeric( DAT_miss[,'X3']*DAT_miss[,'S1']), 
                          I_32 = as.numeric( DAT_miss[,'X3']*DAT_miss[,'S2']), 
                          I_33 = as.numeric( DAT_miss[,'X3']*DAT_miss[,'S3']), 
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
Mean_X1 = mean(unlist(with(impute_MNAR_INTR,mean(X1))$analyses))
RESULTS[RESULTS$METHODS == 'SRMI-Interactions R','Mean_X1'] = Mean_X1
means_long = with(impute_MNAR_INTR,mean(as.numeric(as.character(X1)=='Yes')))$analyses
vars_long = with(impute_MNAR_INTR,var(X1)/Nobs)$analyses
Var_X1 = mean(unlist(vars_long)) + (1+(1/10))*var(unlist(means_long))
RESULTS[RESULTS$METHODS == 'SRMI-Interactions R','Var_X1'] = Var_X1


#####################
### Cubic Splines ###
#####################
pred <- impute_MAR$predictorMatrix
pred[,] = 0
pred["X1",c("X2", "X3", "X4", "X5", "S2", "S3")] <- rep(1,6)
pred["X2",c("X1", "X3", "X4", "X5", "S1", 'S3')] <- rep(1,6)
pred["X3",c("X1", "X2", "X4", "X5", "S1", 'S2')] <- rep(1,6)
DAT_miss_factor = DAT_miss
DAT_miss_factor$S1 = factor(DAT_miss_factor$S1)
DAT_miss_factor$S2 = factor(DAT_miss_factor$S2)
DAT_miss_factor$S3 = factor(DAT_miss_factor$S3)
impute_MNAR_SPLINE <- mice(DAT_miss_factor, m=10, predictorMatrix=pred, method=c("tricube.pmm", "tricube.pmm","tricube.pmm",NA,NA,NA,NA,NA), printFlag=F, remove.constant = FALSE, remove.collinear = FALSE)
Mean_X1 = mean(unlist(with(impute_MNAR_SPLINE,mean(X1))$analyses))
RESULTS[RESULTS$METHODS == 'SRMI-TriCube','Mean_X1'] = Mean_X1
means_long = with(impute_MNAR_SPLINE,mean(X1))$analyses
vars_long = with(impute_MNAR_SPLINE,var(X1)/Nobs)$analyses
Var_X1 = mean(unlist(vars_long)) + (1+(1/10))*var(unlist(means_long))
RESULTS[RESULTS$METHODS == 'SRMI-TriCube','Var_X1'] = Var_X1

#####################################
### Exact Imputation Distribution ###
#####################################
pred <- impute_MAR$predictorMatrix
pred[,] = 0
pred["X1",c("X2", "X3", "X4", "X5", "S2", "S3")] <- rep(1,6)
pred["X2",c("X1", "X3", "X4", "X5", "S1", 'S3')] <- rep(1,6)
pred["X3",c("X1", "X2", "X4", "X5", "S1", 'S2')] <- rep(1,6)
impute_EXACT <- mice(DAT_miss, m=10, predictorMatrix=pred, method=c("normExact", "normExact","normExact",NA,NA,NA,NA,NA), printFlag=F, maxit = 50)
Mean_X1 = mean(unlist(with(impute_EXACT,mean(X1))$analyses))
RESULTS[RESULTS$METHODS == 'SRMI-Exact','Mean_X1'] = Mean_X1
means_long = with(impute_MNAR_SPLINE,mean(X1))$analyses
vars_long = with(impute_MNAR_SPLINE,var(X1)/Nobs)$analyses
Var_X1 = mean(unlist(vars_long)) + (1+(1/10))*var(unlist(means_long))
RESULTS[RESULTS$METHODS == 'SRMI-Exact','Var_X1'] = Var_X1



#######################
### Compare Results ###
#######################
### Note: these are results for a single simulated dataset. The true mean of random variable X1 is 0 (dashed black line). 
### The mean of the simulated data sample without missingness is shown by the gray dashed line.
ggplot(RESULTS, aes(x = Mean_X1, y = METHODS, xmin = Mean_X1-1.96*sqrt(Var_X1), xmax = Mean_X1+1.96*sqrt(Var_X1))) +
  geom_pointrange(aes(fill = METHODS), position = ggstance::position_dodgev(height=1), shape = 21, size = 1.5) +
  geom_vline(xintercept=0,linetype=2)+
  scale_y_discrete(limits = rev(METHOD_NAMES))+
  xlab("Mean of X1") + ylab('Method')+labs(title='Estimated Mean of X1 and 95% Confidence Interval')+
  theme_classic() +theme(legend.position ='none')+
  geom_vline(aes(xintercept = Mean_X1), color = 'gray', linetype = 2, data = RESULTS[RESULTS$METHODS =='No Missingness', ])




