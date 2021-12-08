

PATH= 'C:/Users/ericb/Desktop/Data Science/BSE_Data_Science_master/Statistical Modeling and Inference/Part 1 - Regression, Penalized likelihood, Bayesian regression, Non-parametric regression/StatInference_seminar1'
library(hdm)
library(glmnet)
library(ggplot2)
library(tidyverse)
library(HDCI)
library(gridExtra)
source(file.path(PATH,"routines_seminar1.R"))
source(file.path(PATH,"routines_seminar1 (glm).R"))

spam <- read.table(file.path(PATH, 'data/spam.data'), quote="\"", comment.char="")
spam.names <- c(read.table("data/spambase.names", sep = ":", skip = 33, nrows = 53, as.is = TRUE)[,1],"char_freq_#",
                read.table("data/spambase.names", sep = ":", skip = 87, nrows = 3, as.is = TRUE)[,1], "spam.01")

names(spam) <- spam.names
spam <- spam[sample(nrow(spam)),]

y= spam$spam.01
x= model.matrix(spam.01~.,data=spam) # y~. tells R to use all variables in the dataframe as predictors y(x1,x2,x3,...) except spam.01 that is our target 
dim(x)

#Logistic regression
####################################################################################################
fit.mle= glm(y ~ x[,-1], family = "binomial") #1st column in x is the intercept, already added by lm
b.mle= coef(fit.mle)
summary(fit.mle)
cor(fit.mle$fitted.values,y)^2 #In Sample Efron's Pseudo-R^2=0.7726014

cv.mle= kfoldCV.glm(y=y,x=data.frame(x[,-1]),K=10,seed=1, family = "binomial", predict_type = "response")
r2.mle= cor(cv.mle$pred, y)^2
r2.mle #10-fold cross-validated Efron's Pseudo-R^2 = 0.7526937 

logLik_fit = sum(log(cv.mle$pred**(y)*(1-cv.mle$pred)**(1-y))) #Log likelihood of the fitted model
logLik_null = sum(log((sum(y)/length(y))**(y)*(1-sum(y)/length(y))**(1-y))) #LogLik of Null model where we don't take into account any covariate --> p = spam/emails (independently of any x) 
r2_McFadden.mle = 1-logLik_fit/logLik_null # McFadden Pseudo-r2 = 0.62863
####################################################################################################

#Lasso
####################################################################################################
fit.lasso = cv.glmnet(x=x[,-1], y=y, nfolds=10, family = "binomial")
fit.lasso #Optimal lambda = 0.000335 --> Nonzero = 54 --> Index 69 of the lambdas tried
length(fit.lasso$lambda) #91 lambdas tried
fit.lasso$lambda.min #Optimal lambda = 0.000335
b.lasso= as.vector(coef(fit.lasso, s='lambda.min'))
rownames(coef(fit.lasso, s='lambda.min'))
round(b.lasso, 3)

rownames(coef(fit.lasso, s = 'lambda.min'))[coef(fit.lasso, s = 'lambda.min')[,1]!= 0] # returns names of nonzero coefs
coef(fit.lasso, s = 'lambda.min')[coef(fit.lasso, s = 'lambda.min')[,1]!= 0] # returns nonzero coefs (Takes also the intercept --> 55)

lasso_cv.mle = kfoldCV.lasso.glm(y=y,x=x[,-1], K=10, seed=1, family = "binomial", predict_type = "response")
lasso_r2.mle= cor(lasso_cv.mle$pred, y)^2
lasso_r2.mle #10-fold cross-validated Efron's Pseudo-R^2 = 0.7351376 

logLik_fit = sum(log(lasso_cv.mle$pred**(y)*(1-lasso_cv.mle$pred)**(1-y))) #Log likelihood of the fitted model
logLik_null = sum(log((sum(y)/length(y))**(y)*(1-sum(y)/length(y))**(1-y))) #LogLik of Null model where we don't take into account any covariate --> p = spam/emails (independently of any x) 
lasso_r2_sMcFadden.mle = 1-logLik_fit/logLik_null # McFadden Pseudo-r2 = 0.6417641
####################################################################################################

#Lasso Bic
####################################################################################################
fit.lassobic = lasso.bic.glm(x=x[,-1], y=y, family = "binomial", extended = FALSE)
fit.lassobic #Optimal lambda = 0.000335 --> Nonzero = 55 --> Index 69 of the lambdas tried
length(fit.lassobic$lambda[['lambda']]) #91 lambdas tried
fit.lassobic$lambda.opt #Optimal lambda = 0.0001591292
b.lassobic= as.vector(fit.lassobic[["coef"]])
rownames(coef(fit.lasso, s='lambda.min'))
round(b.lasso, 3)

lassobic_cv.mle = kfoldCV.lasso.glm(y=y,x=x[,-1], K=10, seed=1, family = "binomial", predict_type = "response", criterion="bic")
lassobic_r2.mle= cor(lassobic_cv.mle$pred, y)^2
lassobic_r2.mle #10-fold cross-validated Efron's Pseudo-R^2 = 0.7521989 

logLik_fit = sum(log(lassobic_cv.mle$pred**(y)*(1-lassobic_cv.mle$pred)**(1-y))) #Log likelihood of the fitted model
logLik_null = sum(log((sum(y)/length(y))**(y)*(1-sum(y)/length(y))**(1-y))) #LogLik of Null model where we don't take into account any covariate --> p = spam/emails (independently of any x) 
lassobic_r2_McFadden.mle = 1-logLik_fit/logLik_null # McFadden Pseudo-r2 = 0.643267
####################################################################################################




