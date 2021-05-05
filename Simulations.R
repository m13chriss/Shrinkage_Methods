######################################################################
#########
install.packages("tidyverse")
library(tidyverse)
library(ggplot2)



# Generate Data 
generate = function(sample_size){
  eps = rnorm(n = sample_size, mean = 0, sd = 1)
  beta.vec = c(0.9)
  X = data.frame(X2 = rnorm(n = sample_size, mean = 1, sd = 0.65)) 
  Y = as.matrix(X) %*% beta.vec + eps
  return(data.frame(Y, X))}



residuals = function(sample_size, slope){
  # Generate the sample data
  data = generate(sample_size = sample_size)
  # Calculating sum of squared residuals divided by sample size
  sr = mean((data$Y - as.matrix(data$X) %*% slope)**2)
  return(sr)
}

  
visualization_ols = function(sample_size){
  # Deciding a list of slopes
  slopes = seq(-0.5, 1, by = 0.01)
  # Calculate squared residuals for each one of the slopes
  squared_residuals = c()
  for(s in slopes){
    sr = residuals(sample_size = sample_size, slope = s)
    squared_residuals = append(squared_residuals, sr)
  }
  
  
  
  # The curve
  a = data.frame(slopes, squared_residuals)   
  ggplot(a, aes(x=slopes, y=squared_residuals)) + 
    geom_point(pch = 1) + 
    geom_smooth() +
    labs(x="Slope Values", y= "Sum of Squared Residuals", 
         title = expression(paste("OLS - Sum of squared Residuals  ", lambda == 0))) 
}


visualization_ols(10000)


residuals_lasso = function(sample_size, slope, lambda){
  # Generate the sample data
  data = generate(sample_size = sample_size)
  # Calculating sum of squared residuals + ridge penalty
  lasso = mean((data$Y - as.matrix(data$X) %*% slope)**2) + lambda * abs(slope)         
  return(lasso)
}

residuals_ridge = function(sample_size, slope, lambda){
  # Generate the sample data
  data = generate(sample_size = sample_size)
  # Calculating sum of squared residuals + ridge penalty
  ridge = mean((data$Y - as.matrix(data$X) %*% slope)**2) + lambda * (slope**2)         
  return(ridge)
}

srs = function(sample_size, f){
  # Calculates sum of sqrt residuals for each lambda
  lambda = c(0, 1, 2, 4)
  # Deciding a list of slopes
  slopes = seq(-0.5, 1, by = 0.01)
  # Calculate squared residuals for each one of the slopes
  squared_residuals1 = c()
  squared_residuals2 = c()
  squared_residuals3 = c()
  squared_residuals4 = c()
  for(s in slopes){
    sr1 = f(sample_size, slope = s, lambda = lambda[1])
    sr2 = f(sample_size, slope = s, lambda = lambda[2])
    sr3 = f(sample_size, slope = s, lambda = lambda[3])
    sr4 = f(sample_size, slope = s, lambda = lambda[4])
    squared_residuals1 = append(squared_residuals1, sr1)
    squared_residuals2 = append(squared_residuals2, sr2)
    squared_residuals3 = append(squared_residuals3, sr3)
    squared_residuals4 = append(squared_residuals4, sr4)
  }
  # Gathers sqrt residuals rowwise in a data frame to be used in plot
  vis.data.ridge1 = data.frame(slopes, sr = squared_residuals1)
  vis.data.ridge2 = data.frame(slopes, sr = squared_residuals2)
  vis.data.ridge3 = data.frame(slopes, sr = squared_residuals3)
  vis.data.ridge4 = data.frame(slopes, sr = squared_residuals4)
  data.set = bind_rows(vis.data.ridge1, vis.data.ridge2, vis.data.ridge3, vis.data.ridge4)
  data.set$g = c(rep("lambda = 0", length(slopes)), rep("lambda = 1", length(slopes)),
                 rep("lambda = 2", length(slopes)), rep("lambda = 4", length(slopes)))
  colnames(data.set) <- c("Slopes", "sr", "lambdas")
  
  return(data.set)
}

vis.data.lasso = srs(100000, f = residuals_lasso)

# LASSO 
ggplot(vis.data.lasso, aes(x=Slopes, y=sr, color = lambdas)) + 
  geom_line(size = 1.7) + 
  labs(x="Slope Values", y="Sum of Squared Residuals + lambda * |slope|", 
       title = "LASSO - Sum of squared Residuals\n+ l1 penalty") +
  labs(color = "Lambda Values") +
  scale_color_manual(labels = parse(text = c("lambda == 0", "lambda == 1", "lambda == 2", "lambda == 4")),
                     values = c("violet", "green", "red", "blue")) +
  theme(legend.justification=c(1,0), legend.position=c(1, 0.80))


# Ridge
vis.data.ridge = srs(100000, f = residuals_ridge)
ggplot(vis.data.ridge, aes(x=Slopes, y=sr, color = lambdas)) + 
  geom_line(size = 1.7) + 
  labs(x="Slope Values", y="Sum of Squared Residuals + lambda * slope^2", 
       title = "RIDGE - Sum of squared Residuals\n+ l2 penalty") +
  labs(color = "Lambda Values") +
  scale_color_manual(labels = parse(text = c("lambda == 0", "lambda == 1", "lambda == 2", "lambda == 4")),
                     values = c("violet", "green", "red", "blue")) +
  theme(legend.justification=c(1,0), legend.position=c(1, 0.80))


##############################################################################
##############################################################################
  


library(MASS)  # Package needed to generate correlated precictors
library(glmnet)  # Package to fit ridge/lasso/elastic net models
library(ggplot2)
install.packages("ggforce")
library(ggforce)



########################################
##A) Small Signal and lot noise#########
########################################




# Generate data
set.seed(19)  # Set seed for reproducibility
n <- 500  # Number of observations
p <- 1500  # Number of predictors included in model
real_p <- 15  # Number of true predictors
x <- matrix(rnorm(n*p), nrow=n, ncol=p)
y <- apply(x[,1:real_p], 1, sum)*3/sqrt(real_p) + rnorm(n)



# Split data into train (2/3) and test (1/3) sets
train_rows <- sample(1:n, .66*n)
x.train <- x[train_rows, ]
x.test <- x[-train_rows, ]

y.train <- y[train_rows]
y.test <- y[-train_rows]
lambdas_to_try <- 10^seq(-50, 5, length.out = 100)



# Fit models 
# (For plots on left): 
fit.lasso1 <- glmnet(x.train, y.train, family="gaussian", alpha=1)
fit.ridge1 <- glmnet(x.train, y.train, family="gaussian", alpha=0)
fit.lasso <- cv.glmnet(x.train, y.train, family="gaussian", alpha=1)
fit.ridge <- cv.glmnet(x.train, y.train, family="gaussian", alpha=0, lambda = lambdas_to_try)
yhat10 <- predict(fit.lasso, s=fit.lasso$lambda.1se, newx=x.test)
mse10 <- mean((y.test - yhat10)^2)
mse10
yhat0 <- predict(fit.ridge, s=fit.ridge$lambda.1se, newx=x.test)
mse0 <- mean((y.test - yhat0)^2)
mse0
# Plot solution paths:
par(mfrow = c(2,2))
plot(fit.lasso1, xvar="lambda")
plot(fit.lasso)

plot(fit.ridge1, xvar="lambda")
plot(fit.ridge)



#########################################################################
#########For ELNET selection of Alpha and Lambda (Small Signal)##########
#########################################################################




library(ensr)
ggforce::facet_zoom

ensr_obj <- ensr(y = y.train, x = x.train, standardize = FALSE)
ensr_obj
summary.1=summary(object = ensr_obj)
summary.1[cvm == min(cvm)]
plot(ensr_obj)
plot(ensr_obj) +
  theme_minimal() +
  facet_zoom(x = 0.50 < alpha & alpha < 0.75, y = 1e-03 < lambda & lambda < 1e+00)

scale_index





#################################################################################
#################################################################################




##########################################
###B) Big Signal and lot of noise#########
##########################################



# Generate data
set.seed(2003)
n <- 500   # Number of observations
p <- 1500     # Number of predictors included in model
real_p <- 700  # Number of true predictors
x <- matrix(rnorm(n*p), nrow=n, ncol=p)
y <- apply(x[,1:real_p], 1, sum)*3/sqrt(real_p) + rnorm(n)



# Split data into train and test sets
train_rows <- sample(1:n, .66*n)
x.train <- x[train_rows, ]
x.test <- x[-train_rows, ]

y.train <- y[train_rows]
y.test <- y[-train_rows]
lambdas_to_try <- 10^seq(-50, 5, length.out = 100)
lambdas_to_try2 <- 10^seq(-5, 10, length.out = 100)




# Fit models:
fit.lasso <- glmnet(x.train, y.train, family="gaussian", alpha=1)
fit.ridge <- glmnet(x.train, y.train, family="gaussian", alpha=0)


fit.lasso.cv <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=1, 
                          family="gaussian")
fit.ridge.cv <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=0,
                          family="gaussian", lambda = lambdas_to_try)
yhat10 <- predict(fit.lasso.cv, s=fit.lasso.cv$lambda.1se, newx=x.test)
mse10 <- mean((y.test - yhat10)^2)
mse10
yhat0 <- predict(fit.ridge.cv, s=fit.ridge.cv$lambda.1se, newx=x.test)
mse0 <- mean((y.test - yhat0)^2)
mse0




# Plot solution paths:
par(mfrow=c(2,2))
plot(fit.lasso, xvar="lambda")
plot(fit.lasso.cv)

plot(fit.ridge, xvar="lambda")
plot(fit.ridge.cv)


#########################################################################
#########For ELNET selection of Alpha and Lambda (Big Signal)############
#########################################################################


library(ensr)


ensr_obj <- ensr(y = y.train, x = x.train, standardize = FALSE)
ensr_obj
summary.2=summary(object = ensr_obj)
summary.2[cvm == min(cvm)]
plot(ensr_obj)
plot(ensr_obj) +
  theme_minimal() +
  facet_zoom(x = -0.1 < alpha & alpha < 0.10, y = 1e+00 < lambda & lambda < 1e+02)



##################################################################################
##################################################################################


##########################################################
########C) Correlated X Variables#########################
##########################################################


# Generate data
set.seed(19873)
n <- 100    # Number of observations
p <- 50     # Number of predictors included in model
CovMatrix <- outer(1:p, 1:p, function(x,y) {.7^abs(x-y)})
x <- mvrnorm(n, rep(0,p), CovMatrix)
y <- 10 * apply(x[, 1:2], 1, sum) + 
  5 * apply(x[, 3:4], 1, sum) +
  apply(x[, 5:14], 1, sum) +
  rnorm(n)



# Split data into train and test sets
train_rows <- sample(1:n, .66*n)
x.train <- x[train_rows, ]
x.test <- x[-train_rows, ]

y.train <- y[train_rows]
y.test <- y[-train_rows]



# Fit models:
fit.lasso <- glmnet(x.train, y.train, family="gaussian", alpha=1)
fit.ridge <- glmnet(x.train, y.train, family="gaussian", alpha=0)

fit.lasso.cv <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=1, 
                          family="gaussian")
fit.ridge.cv <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=0,
                          family="gaussian")
yhat10 <- predict(fit.lasso.cv, s=fit.lasso.cv$lambda.1se, newx=x.test)
mse10 <- mean((y.test - yhat10)^2)
mse10
yhat0 <- predict(fit.ridge.cv, s=fit.ridge.cv$lambda.1se, newx=x.test)
mse0 <- mean((y.test - yhat0)^2)
mse0
# Plot solution paths:
par(mfrow=c(2,2))
plot(fit.lasso, xvar="lambda")
plot(fit.lasso.cv)

plot(fit.ridge, xvar="lambda")
plot(fit.ridge.cv)


#####################################################################################
#########For ELNET selection of Alpha and Lambda (Correlated X variables)############
#####################################################################################


library(ensr)


ensr_obj <- ensr(y = y.train, x = x.train, standardize = FALSE)
ensr_obj
summary.1=summary(object = ensr_obj)
summary.1[cvm == min(cvm)]
plot(ensr_obj)
plot(ensr_obj) +
  theme_minimal() +
  facet_zoom(x = 0.6 < alpha & alpha < 1, y = 1e-02 < lambda & lambda < 1e+02)



################################################################################
################################################################################



library(glmnet)
library(ggplot2)




############################################################
######## Loop for varying signal in LASSO/Ridge against MSE# 
############################################################
mse_ridge = c()
mse_lasso = c()

for (i in 1:298) {
  set.seed(19)
  n <- 500    # Number of observations
  p <- 1500     # Number of predictors included in model
  real_p <- 10+ i*5  # Number of true predictors
  x <- matrix(rnorm(n*p), nrow=n, ncol=p)
  y <- apply(x[,1:real_p], 1, sum)*3/sqrt(real_p) + rnorm(n)
  
  # Split data into train and test sets
  train_rows <- sample(1:n, .66*n)
  x.train <- x[train_rows, ]
  x.test <- x[-train_rows, ]
  
  y.train <- y[train_rows]
  y.test <- y[-train_rows]
  
  
  fit.lasso.cv <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=1, 
                            family="gaussian")
  fit.ridge.cv <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=0,
                            family="gaussian")
  yhat_lasso <- predict(fit.lasso.cv, s=fit.lasso.cv$lambda.1se, newx=x.test)
  yhat_ridge <- predict(fit.ridge.cv, s=fit.ridge.cv$lambda.1se, newx=x.test)
  
  mse_lasso_add <- mean((y.test - yhat_lasso)^2)
  mse_ridge_add <- mean((y.test - yhat_ridge)^2)
  
  mse_lasso <- c(mse_lasso, mse_lasso_add)
  mse_ridge <- c(mse_ridge, mse_ridge_add)
  
}

mse_lasso
mse_ridge
#############################################################################
#############################################################################
mses<-data.frame(cbind(mse_lasso, mse_ridge, c(seq(15, 1500, by = 5))))
ggplot(data=df2, aes(x=Real_p, y=mse, group=supp)) +
  geom_line(aes(color=supp))+
  geom_point(aes(color=supp))



#################################################
## Loop for varying signal in ELNET against MSE##
#################################################

mse_elnet = c()

for (i in 1:298) {
  set.seed(19)
  n <- 500    # Number of observations
  p <- 1500     # Number of predictors included in model
  real_p <- 10+ i*5  # Number of true predictors
  x <- matrix(rnorm(n*p), nrow=n, ncol=p)
  y <- apply(x[,1:real_p], 1, sum)*3/sqrt(real_p) + rnorm(n)
  
  # Split data into train and test sets
  train_rows <- sample(1:n, .66*n)
  x.train <- x[train_rows, ]
  x.test <- x[-train_rows, ]
  
  y.train <- y[train_rows]
  y.test <- y[-train_rows]
  
  mse_vec=c()
  for (k in 1:99) {
    fit.elnet = cv.glmnet(x.train, y.train, type.measure="mse", 
                          alpha=k/100,family="gaussian")
    yhat <- predict(fit.elnet, s=fit.elnet$lambda.1se, newx=x.test)
    mse <- mean((y.test - yhat)^2)
    mse_vec = c(mse_vec, mse)
    
  }
  
  mse_elnet=c(mse_elnet, min(mse_vec))
}
mse_elnet
############################################################################
############################################################################







#############################################################
###Plotting the varying results for LASSO, Ridge and ELNET###
#############################################################




mses<-data.frame(cbind(mse_lasso, mse_ridge, mse_elnet, c(seq(15, 1500, by = 5))))


df2 <- data.frame(supp=rep(c("mse_lasso", "mse_ridge", "mse_elnet"), each=298),
                  mse=rep(c(mse_lasso, mse_ridge, mse_elnet)),
                  Real_p=c(seq(15, 1500, by = 5)))
head(df2)
library(ggplot2)
ggplot(data=df2, aes(x=Real_p, y=mse, group=supp)) +
  geom_line(aes(color=supp))+
  geom_point(aes(color=supp))
