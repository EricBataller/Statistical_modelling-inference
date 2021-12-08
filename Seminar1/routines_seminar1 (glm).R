#########################################################################################
##
## R ROUTINES FOR SEMNINAR 1 - STATISTICAL MODELLING AND INFERENCE.
##
## AUTHOR: DAVID ROSSELL, PAUL ROGNON, UNIVERSITAT POMPEU FABRA
##
#########################################################################################

# INDEX
# 1. SETTING PENALIZATION PARAMETER VIA BIC
# 2. CROSS-VALIDATION
# 3. LASSO POST-SELECTION INFERENCE
# 4. ROUTINES TO SIMULATE DATA

#########################################################################################
## 1. SETTING PENALIZATION PARAMETER VIA BIC
#########################################################################################

lasso.bic <- function(y,x,extended=FALSE) {
  #Select model in LASSO path with best BIC (using LASSO regression estimates)
  #Input
  # - y: vector with response variable
  # - x: design matrix
  # - extended: whether to use EBIC (Chen and Chen 2008) instead of BIC
  #
  #Output: list with the following elements
  # - coef: LASSO-estimated regression coefficient with lambda set via BIC
  # - ypred: predicted y
  # - lambda.opt: optimal value of lambda
  # - lambda: data.frame with bic and number of selected variables for each value of lambda
  require(glmnet)
  fit <- glmnet(x=x,y=y,family='gaussian',alpha=1)
  pred <- cbind(1,x) %*% rbind(fit$a0,fit$beta)
  n <- length(y)
  p <- colSums(fit$beta!=0) + 1
  if (!extended){
    bic <- n * log(colSums((y-pred)^2)/length(y)) + n*(log(2*pi)+1) + log(n)*p 
  } else {
    bic <- n * log(colSums((y-pred)^2)/length(y)) + n*(log(2*pi)+1) + log(n)*p + 2*log(choose(ncol(x),p))
  }
  sel <- which.min(bic)
  beta <- c(fit$a0[sel],fit$beta[,sel]); names(beta)[1]= 'Intercept'
  ypred <- pred[,sel]
  ans <- list(coef=beta,ypred=ypred,lambda.opt=fit$lambda[sel],lambda=data.frame(lambda=fit$lambda,bic=bic,nvars=p))
  return(ans)
}

lasso.bic.glm <- function(y,x,extended=FALSE, family='gaussian') {
  #Select model in LASSO path with best BIC (using LASSO regression estimates)
  #Input
  # - y: vector with response variable
  # - x: design matrix
  # - extended: whether to use EBIC (Chen and Chen 2008) instead of BIC
  # - family: exponential family distribution to fit model with (for logit use binomial)
  #Output: list with the following elements
  # - coef: LASSO-estimated regression coefficient with lambda set via BIC
  # - ypred: predicted y
  # - lambda.opt: optimal value of lambda
  # - lambda: data.frame with bic and number of selected variables for each value of lambda
  require(glmnet)
  fit <- glmnet(x=x,y=y,family=family,alpha=1)
  pred <- cbind(1,x) %*% rbind(fit$a0,fit$beta)
  n <- length(y)
  p <- colSums(fit$beta!=0) + 1
  if (!extended){
# Deviance - gives -2*Log likelihood
    bic <- deviance(fit) + log(n)*p 
  } else {
    bic <- deviance(fit) + log(n)*p + 2*log(choose(ncol(x),p))
  }
  sel <- which.min(bic)
  beta <- c(fit$a0[sel],fit$beta[,sel]); names(beta)[1]= 'Intercept'
  ypred <- pred[,sel]
  ans <- list(coef=beta,ypred=ypred,lambda.opt=fit$lambda[sel],lambda=data.frame(lambda=fit$lambda,bic=bic,nvars=p))
  return(ans)
}

#########################################################################################
## 2. CROSS-VALIDATION
#########################################################################################

kfoldCV.mle <- function(y,x,K=10,seed) {
## Perform K-fold cross-validation for least-squares regression estimate
## Input
## - y: response
## - x: data.frame with predictors, intercept should not be present
## - K: number of folds in K-fold cross-validation
## - seed: random number generator seed (optional)
## Output
## - pred: cross-validated predictions for y
## - ssr: residual sum of squares, sum((y-pred)^2)
  if (!missing(seed)) set.seed(seed)
  subset <- rep(1:K,ceiling(nrow(x)/K))[1:nrow(x)]
  subset <- sample(subset,size=nrow(x),replace=FALSE)
  pred <- double(nrow(x))
  if (ncol(x)>0) {
    for (k in 1:K) {
      sel <- subset==k
      fit <- lm(y[!sel] ~ ., data=x[!sel,,drop=FALSE])
      pred[sel] <- predict(fit, newdata=x[sel,,drop=FALSE])
    }
  } else {
    for (k in 1:K) {
      sel <- subset==k
      pred[sel] <- mean(y[!sel],na.rm=TRUE)
    }
  }
  return(list(pred=pred,ssr=sum((pred-y)^2,na.rm=TRUE)))
}

kfoldCV.glm <- function(y,x,K=10,seed, family = "gaussian", predict_type = "link") {
  ## Perform K-fold cross-validation for least-squares regression estimate
  ## Input
  ## - y: response
  ## - x: data.frame with predictors, intercept should not be present
  ## - K: number of folds in K-fold cross-validation
  ## - family: generalises to exponential family regressions
  ## - predict_type: type of prediction. For logit use 'response' to get probabilities.
  ## - seed: random number generator seed (optional)
  ## Output
  ## - pred: cross-validated predictions for y
  ## - ssr: residual sum of squares, sum((y-pred)^2)
  if (!missing(seed)) set.seed(seed)
  subset <- rep(1:K,ceiling(nrow(x)/K))[1:nrow(x)]
  subset <- sample(subset,size=nrow(x),replace=FALSE)
  pred <- double(nrow(x))
  if (ncol(x)>0) {
    for (k in 1:K) {
      sel <- subset==k
      fit <- glm(y[!sel] ~ ., data=x[!sel,,drop=FALSE], family = family)
      pred[sel] <- predict.glm(fit, newdata=x[sel,,drop=FALSE], type = predict_type)
    }
  } else {
    for (k in 1:K) {
      sel <- subset==k
      pred[sel] <- mean(y[!sel],na.rm=TRUE)
    }
  }
  return(list(pred=pred,ssr=sum((pred-y)^2,na.rm=TRUE)))
}


kfoldCV.lasso <- function(y,x,K=10,seed,criterion='cv') {
## Perform K-fold cross-validation for LASSO regression estimate (lambda set either via cross-val or BIC or EBIC)
## Input
## - y: response
## - x: data.frame with predictors, intercept should not be present
## - K: number of folds in K-fold cross-validation
## - seed: random number generator seed (optional)
## - criterion: the criterion to select the penalization parameter, either cross-val or BIC or EBIC
## Output
## - pred: cross-validated predictions for y
## - ssr: residual sum of squares, sum((y-pred)^2)
  require(glmnet)
  if (!missing(seed)) set.seed(seed)
  subset <- rep(1:K,ceiling(nrow(x)/K))[1:nrow(x)]
  subset <- sample(subset,size=nrow(x),replace=FALSE)
  pred <- double(nrow(x))
  cat("Starting cross-validation")
  if (ncol(x)>0) {  #if there are some covariates
    for (k in 1:K) {
        sel <- subset==k
        if (criterion=='cv') {
            fit <- cv.glmnet(x=x[!sel,,drop=FALSE], y=y[!sel], alpha = 1, nfolds=10)
            b= as.vector(coef(fit,s='lambda.min'))
            pred[sel] <- b[1] + x[sel,,drop=FALSE] %*% as.matrix(b[-1])
        } else if (criterion=='bic'){
            fit <- lasso.bic(y=y[!sel],x=x[!sel,,drop=FALSE])
            pred[sel] <- fit$coef[1] + x[sel,,drop=FALSE] %*% matrix(fit$coef[-1],ncol=1)
        } else if (criterion=='ebic'){
            fit <- lasso.bic(y=y[!sel],x=x[!sel,,drop=FALSE],extended = TRUE)
            pred[sel] <- fit$coef[1] + x[sel,,drop=FALSE] %*% matrix(fit$coef[-1],ncol=1)
        } else if (criterion=='bayes'){
          fit.bayesreg <- modelSelection(x=x[!sel,,drop=FALSE], y=y[!sel], priorCoef=zellnerprior(taustd=1),
                                         priorDelta=modelbbprior(1,1))
          intercept <- coef(fit.bayesreg)[1,1]
          pred[sel] <- predict(fit.bayesreg, newdata= x[sel,,drop=FALSE])[,1] + intercept
        } else { stop("method.lambda not implemented") }
        cat(".")
    }
  } else { #if there are no covariates, just use the intercept
    for (k in 1:K) {
      sel <- subset==k
      pred[sel] <- mean(y[!sel],na.rm=TRUE)
    }
  }
  cat("\n")
  return(list(pred=pred,ssr=sum((pred-y)^2,na.rm=TRUE)))
}

kfoldCV.lasso.glm <- function(y,x,K=10,seed,criterion='cv', family = 'gaussian', predict_type = "link") {
  ## Perform K-fold cross-validation for LASSO regression estimate (lambda set either via cross-val or BIC or EBIC)
  ## Input
  ## - y: response
  ## - x: data.frame with predictors, intercept should not be present
  ## - K: number of folds in K-fold cross-validation
  ## - seed: random number generator seed (optional)
  ## - family: generalises to exponential family regressions
  ## - predict_type: type of prediction. For logit use 'response' to get probabilities (see ?predict.glmnet for options).
  ## - criterion: the criterion to select the penalization parameter, either cross-val or BIC or EBIC
  ## Output
  ## - pred: cross-validated predictions for y
  ## - ssr: residual sum of squares, sum((y-pred)^2)
  require(glmnet)
  if (!missing(seed)) set.seed(seed)
  subset <- rep(1:K,ceiling(nrow(x)/K))[1:nrow(x)]
  subset <- sample(subset,size=nrow(x),replace=FALSE)
  pred <- double(nrow(x))
  cat("Starting cross-validation")
  if (ncol(x)>0) {  #if there are some covariates
    for (k in 1:K) {
      sel <- subset==k
      if (criterion=='cv') {
        fit <- cv.glmnet(x=x[!sel,,drop=FALSE], y=y[!sel], family = family, nfolds=10)
        pred[sel] <- predict(fit, newx= x[sel,,drop=FALSE], type = predict_type) 
      } else if (criterion=='bic'){
        fit <- lasso.bic.glm(y=y[!sel],x=x[!sel,,drop=FALSE], family = family, extended = FALSE)
        model <- glmnet(x=x[!sel,,drop=FALSE],y=y[!sel],lambda = fit$lambda.opt,family=family)
        pred[sel] <- predict(model, newx= x[sel,,drop=FALSE], type = predict_type)
      } else if (criterion=='ebic'){
        fit <- lasso.bic.glm(y=y[!sel],x=x[!sel,,drop=FALSE], family = family, extended = TRUE)
        model <- glmnet(x=x[!sel,,drop=FALSE],y=y[!sel],lambda = fit$lambda.opt,family=family)
        pred[sel] <- predict(model, newx= x[sel,,drop=FALSE], type = predict_type)
      } else { stop("method.lambda not implemented") }
      cat(".")
    }
  } else { #if there are no covariates, just use the intercept
    for (k in 1:K) {
      sel <- subset==k
      pred[sel] <- mean(y[!sel],na.rm=TRUE)
    }
  }
  cat("\n")
  return(list(pred=pred,ssr=sum((pred-y)^2,na.rm=TRUE)))
}

kfoldCV.ridge <- function(y,x,K=10,seed,criterion='cv') {
  ## Perform K-fold cross-validation for Ridge regression estimate (lambda set via cross-val)
  ## Input
  ## - y: response
  ## - x: data.frame with predictors, intercept should not be present
  ## - K: number of folds in K-fold cross-validation
  ## - seed: random number generator seed (optional)
  ## - criterion: the criterion to select the penalization parameter, (only cross-val for now)
  ## Output
  ## - pred: cross-validated predictions for y
  ## - ssr: residual sum of squares, sum((y-pred)^2)
  require(glmnet)
  if (!missing(seed)) set.seed(seed)
  subset <- rep(1:K,ceiling(nrow(x)/K))[1:nrow(x)]
  subset <- sample(subset,size=nrow(x),replace=FALSE)
  pred <- double(nrow(x))
  cat("Starting cross-validation")
  if (ncol(x)>0) {  #if there are some covariates
    for (k in 1:K) {
      sel <- subset==k
      if (criterion=='cv') {
        fit <- cv.glmnet(x=x[!sel,,drop=FALSE], y=y[!sel], alpha = 0, nfolds=10)
        b= as.vector(coef(fit,s='lambda.min'))
        pred[sel] <- b[1] + x[sel,,drop=FALSE] %*% as.matrix(b[-1])
      } else { stop("method.lambda not implemented") }
      cat(".")
    }
  } else { #if there are no covariates, just use the intercept
    for (k in 1:K) {
      sel <- subset==k
      pred[sel] <- mean(y[!sel],na.rm=TRUE)
    }
  }
  cat("\n")
  return(list(pred=pred,ssr=sum((pred-y)^2,na.rm=TRUE)))
}


#########################################################################################
## 3. LASSO POST-SELECTION INFERENCE
#########################################################################################

lassopost  <- function(y,x, method.lambda='bic', verbose = FALSE) {
  #Run LASSO + post-selection inference on selected coefficients. lambda set via 10-fold cross-validation
  #Input
  # - y: response variable
  # - x: design matrix
  #
  # Output: matrix with parameter estimates, 95% CIs and P-values for the variables selected by LASSO (lambda set via cross-validation), using the post-selection inference method of Lee et al (2016)
  require(glmnet)
  require(selectiveInference)
  x= x[,apply(x,2,'sd') > 0] #remove intercept
  xstd= scale(x); ystd= scale(y) #important!
  if (is.null(colnames(x))) {
    colnames(xstd)= paste('X',1:ncol(xstd),sep='')
  } else {
    colnames(xstd)= colnames(x)
  }
  #Estimate residual variance
  sigmahat= estimateSigma(x=xstd, y=ystd, standardize=FALSE)$sigmahat
  #Set lambda
  if (method.lambda == 'cv') {
    cvfit= cv.glmnet(x=xstd, y=ystd, nfolds=10)
    lambda= cvfit$lambda.min
  } else if (method.lambda == 'bic') {
    lambda= lasso.bic(y=ystd, x=xstd)$lambda.opt
  } else { stop("method.lambda not implemented") }
  #LASSO estimate for selected lambda            
  gfit= glmnet(x=xstd, y=ystd, standardize=FALSE, lambda=lambda)
  b= coef(gfit)[-1]
  #Post-selection inference
  lcv= fixedLassoInf(x=xstd,y=ystd,beta=b,lambda=lambda*nrow(x),sigma=sigmahat,verbose=verbose)
  sel.lasso= (b!=0)
  ci.lassopost= matrix(0,nrow=length(b),ncol=4)
  colnames(ci.lassopost)= c('estimate','ci.low','ci.up','pvalue'); rownames(ci.lassopost)= colnames(xstd)
  ci.lassopost[b != 0,]= cbind(lcv$coef0, lcv$ci, lcv$pv)
  ci.lassopost[b == 0,4]= NA
  return(ci.lassopost)
}


#########################################################################################
## 4. ROUTINES TO SIMULATE DATA
#########################################################################################


simdata.cs <- function(seed,n,theta,phi=1,rho) {
    #Simulate n observations from a linear model y= X %*% theta + e, where e ~ N(0,phi)
    # and the rows in X ~ Normal(0,Sigma) and Sigma is compound symmetric (Sigma[i,i]=1, Sigma[i,j]=rho)
    #Input
    # - seed: random number seed
    # - n: sample size
    # - theta: true regression coefficients
    # - phi: residual variance
    # - rho: true pairwise correlation between variables
    #Output: list with the following elements
    # - y: response
    # - x: n * length(theta) matrix with predictors
    require(mvtnorm)
    S <- diag(length(theta))
    S[upper.tri(S)] <- S[lower.tri(S)] <- rho
    set.seed(seed)
    x <- rmvnorm(n=n,sigma=S)
    y <- x %*% theta + rnorm(n,sd=sqrt(phi))
    return(list(y=y,x=x))
}



