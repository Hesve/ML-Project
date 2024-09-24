## ----echo=FALSE, results=FALSE-------------------------------------------------------------------------------------
knitr::opts_chunk$set(
 comment = '', fig.width = 6, fig.height = 6, eval=TRUE, echo=TRUE, warning=FALSE, message=FALSE
)


## ------------------------------------------------------------------------------------------------------------------
#Libraries
 library(tidyverse)
 library(xtable)
 library(uuml)
 library(gridExtra)
 library(glmnet)
 library(parallel)



## ------------------------------------------------------------------------------------------------------------------

data("binary")
binary$gre_sd <- (binary$gre - mean(binary$gre)) /sd(binary$gre)
binary$gpa_sd <- (binary$gpa - mean(binary$gpa)) /sd(binary$gpa)
X <- model.matrix(admit ~ gre_sd + gpa_sd, binary)
y <- binary$admit


## ------------------------------------------------------------------------------------------------------------------
#' gradient function
#' @param y  vector of observations outcomes
#' @param X Data Matrix
#' @param theta Parameters
#' @return The gradient

# Gradient function
ll_grad <- function(y, X, theta) {
  n <- length(y)
  P <- 1 / (1 + exp(-X %*% theta))
  grad <- t(X) %*% (y-P) / n
  return(grad)
}


round(t(ll_grad(y, X, theta = c(0, 0, 0))), 4)
round(t(ll_grad(y, X, theta = c(-1, 0.5, 0.5))), 4)


## ------------------------------------------------------------------------------------------------------------------
theta <- c(0, 0, 0) 
ll(y, X, theta)
glm(y ~X  -1 , family = binomial(link = "logit"))$coefficients %>% 
  print()



## ----ordinary gradient descent-------------------------------------------------------------------------------------
#' Batch gradient descent
#' @param  learn_rate the step size
#' @param epochs number of iterations

batch_gsd <- function(y, X, theta = rep(0, ncol(X)), learn_rate, epochs) {
    results <- matrix(0.0, ncol = ncol(X) + 2L, nrow = epochs)
  colnames(results) <- c("epoch", "nll", colnames(X))
  for (epoch in 1:epochs) {
    gradient <- -ll_grad(y, X, theta) #negative gradient
    theta <- theta - (learn_rate * gradient)
    
    results[epoch, "epoch"] <- epoch
    results[epoch, "nll"] <- ll(y, X, theta)
    results[epoch, -(1:2)] <- theta
  }
  attributes(results)$learn_rate <- learn_rate
  attributes(results)$epochs <- epochs
  return(results)
}



## ----Stochastic gradient descent-----------------------------------------------------------------------------------
stochastic_gsd <- function(y, X, theta = rep(0, ncol(X)), 
                                   learn_rate, epochs, seed=1337 ){
  set.seed(seed)
  results <- matrix(0.0, ncol = ncol(X) + 2L, nrow = epochs)
  colnames(results) <- c("epoch", "nll", colnames(X))
    for (epoch in 1:epochs) {
      index_order <- sample(length(y)) 
      #shuffle the data by indexing in a random order
      
      for (i in index_order ){
    gradient <- -ll_grad(y[i], X[i,, drop=FALSE], theta) 
    #dont drop "redundant" information to keep it as a matrix
   theta <- theta - (learn_rate * gradient)
      }
    results[epoch, "epoch"] <- epoch
    results[epoch, "nll"] <- ll(y, X, theta)
    results[epoch, -(1:2)] <- theta
      
    }
  attributes(results)$learn_rate <- learn_rate
  attributes(results)$epochs <- epochs
  return(results)
}



## ----minibatch gradient descent------------------------------------------------------------------------------------
minibatch_gsd <- function(y, X, theta = rep(0, ncol(X)), 
                          sample_size, learn_rate, epochs, 
                          batch_size, seed=1337){
  set.seed(seed)
  n <- length(y)
  num_batches <- ceiling(n/batch_size)

results <- matrix(0.0, ncol = ncol(X) + 2L, nrow = epochs)
colnames(results) <- c("epoch", "nll", colnames(X))

# Run algorithm
for(epoch in 1:epochs){
 index_order <- sample(length(y)) 
 #shuffle the data by indexing in a random order
 
### Put the algorithm code here ### 
  for (i in 1:num_batches){

     batch_index <- which(findInterval(index_order,
                                       seq(1, n, by = batch_size)) == i)
     #which x and y  to subset for each iteration by 
     #dividing indexes to different batches
    batch_X <- X[batch_index,]
    batch_y <- y[batch_index]
    
    gradient <- -ll_grad(y=batch_y, X = batch_X, theta=theta ) / batch_size
   theta <- theta - (learn_rate * gradient)
  }
 # Store epoch, nll and output results
results[epoch, "epoch"] <- epoch
results[epoch, "nll"] <- ll(y, X, theta)
results[epoch, -(1:2)] <- theta
}

attributes(results)$learn_rate <- learn_rate
attributes(results)$epochs <- epochs
return(results)
}



## ----echo = FALSE--------------------------------------------------------------------------------------------------
#' Function to plot negative log likelihood and the results for one chosen 
#' parameter
#' @param  variable_name name of the parameter to plot results for
#' @return  plots for the negative log-likelihood as well as the parameter 
#' of interest
#'
plot_GD <- function(results, variable_name ){
  variable_name #this is just to instantiate the object within the local
 # environment due to some weird bug later with the plot_GD_multi function 
  n_epochs <- attributes(results)$epochs
  learn_rate <- attributes(results)$learn_rate
  result_df <- data.frame(results)
  
# Check for Adam implementation attributes (added later).
attributes_to_check <- c("B1", "B2")

if ("B1" %in% names(attributes(results)) && "B2" %in% names(attributes(results))) {
    B1 <- attributes(results)$B1
    B2 <- attributes(results)$B2
    title_text <- bquote(eta == .(learn_rate) ~
                           ", " ~ beta[1] == .(B1) ~ ", " ~ beta[2] == .(B2))
    
  } else {
    title_text <- bquote(eta == .(learn_rate))
  }

y_var <- ensym(variable_name) # ensym + injection with !! to inject the variable
#rather than the string provided by variable_name

   #plot negative Log-likelihood
   ll_plot <- ggplot(result_df, aes(x=epoch, y=nll)) + geom_line() +
     ylab("Negative Log-Likelihood") +
     ggtitle(title_text) + 
     theme(axis.text = element_text(size = 12))
  variable_plot <-   ggplot(result_df, aes(x = epoch, 
                                           y = !!y_var)) +
    geom_line() + 
    geom_hline(yintercept = 0.31, col="blue",linewidth=1, linetype = "dashed") +
    ggtitle(title_text) + 
    theme(axis.text = element_text(size = 12))
  
   plots <- grid.arrange(ll_plot, variable_plot, ncol=2)
  return(plots)
}



## ----echo = FALSE--------------------------------------------------------------------------------------------------
#' wrapper functions to make the previous functions more
#'  generalized to multiple values of learn_rates
#' @param learn_rates a vector of different learning rates to try. 
#' All other arguments are the same

batch_gsd_multi <- function(y, X, learn_rates,
                            thetas = rep(0, ncol(X)), n_epochs){
  res <- mapply(FUN = batch_gsd, learn_rate=learn_rates,
                MoreArgs = list(y=y, X=X, 
                                theta = thetas, epochs = n_epochs), 
                SIMPLIFY = FALSE)
return(res)
}

stochastic_gsd_multi <- function(y, X, learn_rates, 
                                 thetas = rep(0, ncol(X)), n_epochs, seed=1337){
  res <- mapply(FUN = stochastic_gsd, learn_rate=learn_rates,
                MoreArgs = list(y=y, X=X, 
                                theta = thetas, epochs = n_epochs,
                                seed=seed), 
                SIMPLIFY = FALSE)
  return(res)
}

minibatch_gsd_multi <- function(y, X, learn_rates, 
                                 thetas = rep(0, ncol(X)), n_epochs,
                                seed=1337, batch_size) {
  
  res <- mapply(FUN = minibatch_gsd, learn_rate=learn_rates,
                MoreArgs = list(y=y, X=X, 
                                theta = thetas, epochs = n_epochs, 
                                seed=seed, batch_size = batch_size), 
                SIMPLIFY = FALSE)
  return(res)
  }
  

#' function to generalize the plot_GD function to multiple plots in the same grid

plot_GD_multi <-function(results_list, variable_name){
  extracted_plots <- lapply(results_list, FUN= plot_GD,
                            variable_name = variable_name)
  
  merged_plots <- do.call(grid.arrange, c(extracted_plots, 
                                          nrow=length(extracted_plots)))
  #using grid.arrange on all extracted plots
  return(merged_plots)
}




## ----task_1.3_main, echo=FALSE-------------------------------------------------------------------------------------

Task_1.3_main <- function(y, X, thetas = rep(0, ncol(X)) , learn_rates, 
                        n_epochs, seed, batch_size, variable_name){
  
batch_results <- batch_gsd_multi(y = y, X=X,
                                 learn_rates = learn_rates[[1]],
                                 n_epochs = 500, thetas=thetas)
print("Ordinary GD done")

stochastic_results <- stochastic_gsd_multi(y = y, X=X, 
                                           learn_rates = learn_rates[[2]],
                                           n_epochs = n_epochs, seed=seed)
print("Stochastic GD done")

minibatch_results <-  minibatch_gsd_multi(y = y, X=X, 
                                          learn_rates = learn_rates[[3]],
                                          n_epochs = n_epochs,
                                          batch_size=batch_size,
                                          seed=seed)
print("Minibatch GD done")
batch_plots <- plot_GD_multi(batch_results, 
                             variable_name = variable_name)
stochastic_plots <- plot_GD_multi(stochastic_results, 
                                  variable_name = variable_name)
minibatch_plots <- plot_GD_multi(minibatch_results, 
                                 variable_name = variable_name)
list_results <- list(batch_GD =batch_plots, 
                     stochastic_GD = stochastic_plots,
                     minibatch_GD = minibatch_plots)
return(list_results)
}



## ----fig.show="hide", results="hide"-------------------------------------------------------------------------------

plots <- Task_1.3_main(y=y,X=X, learn_rates = list(c(0.03,1,50),
                                                 c(0.0001, 0.001, 1),
                                                 c(0.05,1,15)),
                     n_epochs=500, seed=1337,
                     batch_size = 25, 
                     variable_name = "gre_sd")


## ----batch_results, fig.cap="Results for ordinary gradient descent", fig.width=8, fig.height=8---------------------
plot(plots$batch_GD)


## ----stochastic_results, fig.cap = "Results for stochastic gradient descent", fig.width=8, fig.height=8------------
plot(plots$stochastic_GD)


## ----minibatch_results, fig.cap ="Results for minibatch gradient descent", fig.width=8, fig.height=8---------------
plot(plots$minibatch_GD)


## ------------------------------------------------------------------------------------------------------------------
fit <- lm( y ~ ., data=prob2_train)


## ----eval=FALSE----------------------------------------------------------------------------------------------------
## Residuals:
## ALL 200 residuals are 0: no residual degrees of freedom!
## 
## Residual standard error: NaN on 0 degrees of freedom
## Multiple R-squared:      1,	Adjusted R-squared:    NaN
## F-statistic:   NaN on 199 and 0 DF,  p-value: NA
## 


## ------------------------------------------------------------------------------------------------------------------
X <- as.matrix(prob2_train[, -241])
y <- as.matrix(prob2_train[, "y"])
X_test <- as.matrix(prob2_test[, -241]) 
y_test <- as.matrix(prob2_test[, "y"])
lasso_fit <- glmnet(x=X, y=y, family = "gaussian", lambda=1)

lasso_coeffs <- coef(lasso_fit)


## ------------------------------------------------------------------------------------------------------------------
head(lasso_coeffs)
tail(lasso_coeffs)


## ------------------------------------------------------------------------------------------------------------------
used_coefs <- lasso_coeffs[lasso_coeffs@i+1,] #@ due to S4 class
print(used_coefs)
length(used_coefs)


## ----eval=FALSE----------------------------------------------------------------------------------------------------
## findInterval(cut(sample(1:n,n),breaks=nfolds),1:n)


## ------------------------------------------------------------------------------------------------------------------
#' function to calculate MSE and RMSE
errors_func <- function(y_obs, pred){
  assertthat::assert_that(length(y_obs) == length(pred), msg="lengths differ")
  n <- length(y_obs)
  Sum_square <- sum( (y_obs - pred)^2)
  MSE <- Sum_square / n
  RMSE <- sqrt(MSE)
  res <- list("MSE"= MSE, "RMSE" = RMSE)
  return(res)
}



## ------------------------------------------------------------------------------------------------------------------
cross_validation <- function(X, y, n_folds, lambda, alpha=1, seed=1337){
  assertthat::are_equal(nrow(X),length(y), msg="X and y dimension mismatch")
  checkmate::assert_int(n_folds)
  set.seed(seed)
  
  n <- length(y)
  errors_matrix <- matrix(0, nrow= n_folds, ncol=2)
  colnames(errors_matrix) <- c("MSE", "RMSE")

  Group= findInterval(cut(sample(1:n,n),breaks=n_folds),1:n) 
  #making n_folds equal groups.
  for (i in 1:n_folds){
  index <- which(Group==i) 
  #index for which observations belongs to the current fold group
  
  X_test <- X[index,]
  X_train <- X[-index,]
  y_test <- y[index]
  y_train <- y[-index]
  model <-  glmnet(x=X_train, y=y_train, family = "gaussian", 
                     lambda=lambda, alpha=alpha )
  predictions <- predict(model, s=lambda, newx=X_test, type="response")
  
  fold_errors <- errors_func(y_obs = y_test, pred= predictions )
  errors_matrix[i,1] <- fold_errors$MSE
  errors_matrix[i,2] <- fold_errors$RMSE
  
  
  }
  means <- apply(errors_matrix, MARGIN=2, FUN=mean)
  return(means)
}



## ------------------------------------------------------------------------------------------------------------------
results <- cross_validation(X= X,y= y, n_folds=10, lambda=1)
print(results)




## ------------------------------------------------------------------------------------------------------------------
results2 <- cross_validation(X= X,y= y, n_folds=5, lambda=1)
print(results2)


## ------------------------------------------------------------------------------------------------------------------
#' @param export_vec vector of objects or functions to export to each cluster

cross_validation_parallel <- function(X, y, n_folds, lambdas, alpha=1, seed=1337, export_vec){
  n_cores <- parallel::detectCores()- 1 #saving 1 core for OS 

  cl <- makeCluster(n_cores)
  clusterExport(cl, export_vec)

  results <- t(parSapply(cl, lambdas, function(lambda) {
    cross_validation(X, y, n_folds, lambda, alpha, seed)
  }))
  
  
  row_names <- t(parSapply(cl, lambdas, function(lambda){paste("lambda = ", lambda)}))
  rownames(results) <- row_names
  stopCluster(cl) #putting cores to sleep
  return(results)
}


## ------------------------------------------------------------------------------------------------------------------
lambda_values <- seq(0,1, by=0.01)
lambda_cv <- cross_validation_parallel(X=X, y=y, n_folds=10, lambdas=lambda_values, 
                          export_vec = c("X", "y", "cross_validation", 
                                         "glmnet", "errors_func") )



## ----fig.cap="RMSE for different lambda values"--------------------------------------------------------------------
plot(y=lambda_cv[,2], x=lambda_values, ylab = "RMSE")


## ----eval=FALSE----------------------------------------------------------------------------------------------------
## #creating a sequence from 0 to 0.2 in 100 000 pieces.
## new_lambdas <- seq(0,0.2, length.out = 100000)
## new_lambda_cv <- cross_validation_parallel(X=X, y=y, n_folds=10, lambdas=new_lambdas,
##                           export_vec = c("X", "y", "cross_validation",
##                                          "glmnet", "errors_func") )
## 
## #find the minimum
## which(new_lambda_cv == min(new_lambda_cv[,2]), arr.ind = TRUE)
## multi_res[15,]


## ----eval=FALSE----------------------------------------------------------------------------------------------------
## which(new_lambda_cv == min(new_lambda_cv[,2]), arr.ind = TRUE)


## ----eval=FALSE----------------------------------------------------------------------------------------------------
##                                row col
## lambda =  0.0717627176271763 35882   2


## ----eval = FALSE--------------------------------------------------------------------------------------------------
## new_lambda_cv[35882,]
## 
##  MSE      RMSE
## 0.5612816 0.7456794


## ------------------------------------------------------------------------------------------------------------------
best_lambda_model <- glmnet(x = X, y= y, family = "gaussian", 
                            alpha=1, lambda= 0.0717627176271763)
predictions <- predict(best_lambda_model, newx = X_test, type="response")

errors_func(y=y_test, pred= predictions)

#just to compare that it gives the same results
uuml::rmse(y_test, predictions)



## ------------------------------------------------------------------------------------------------------------------
Adam <- function(y, X, eta = 0.001, B1 = 0.9, B2 = 0.999,
                 epsilon = 10^-8, theta = rep(0, ncol(X)), epochs, seed=1337){
  assertthat::assert_that(length(y) == nrow(X), msg="y and X length differ")
  assertthat::assert_that(( B1 >= 0 & B1 < 1) & (B2 >=0 & B1 <1 ), 
                          msg = "invalid Beta value")
  set.seed(seed)
  
  results <- matrix(0.0, ncol = ncol(X) + 2L, nrow = epochs)
  colnames(results) <- c("epoch", "nll", colnames(X))
  
  m <- 0 #initalise 1st moment vector
  v <- 0 #initalize 2nd moment vector
  
  for (epoch in 1:epochs) { 
      index_order <- sample(length(y))
      X <- X[index_order,]
      y <- y[index_order]
      
      #shuffle the data by indexing in a random order
    gradient <- -ll_grad(y, X, theta) # neg.gradients w.r.t. stochastic objective at timestep t
     
    m <- B1*m + (1-B1)*gradient # Update biased first moment estimate)
    v <- B2*v + (1-B2)*gradient^2 # Update biased second raw moment estimate)
    
    m_hat <- m/(1-B1^epoch) # Compute bias-corrected first moment estimate)
    
    v_hat <- v/(1-B2^epoch) #Compute bias-corrected second raw moment estimate)
    
    theta <- theta - ( eta *m_hat / ( sqrt(v_hat) + epsilon) ) 
    # Update parameters)
    
    results[epoch, "epoch"] <- epoch
    results[epoch, "nll"] <- ll(y, X, theta)
    results[epoch, -(1:2)] <- theta
  }
  attributes(results)$B1 <- B1
  attributes(results)$B2 <- B2
  attributes(results)$learn_rate<-eta
  attributes(results)$epochs <- epochs
  return(results)
}

Adam_multi <- function(y, X, etas, B1s, B2s, epsilon, theta, epochs, seed){
  result_list <- mapply(FUN = Adam, etas, B1s, B2s, 
         MoreArgs = list(y = y, X = X, epsilon = epsilon,
                         theta = theta, epochs = epochs, seed = seed), 
         SIMPLIFY = FALSE)
  return(result_list)
}


task_3_main <- function(y, X, etas, B1s, B2s, 
                        theta = rep(0, ncol(X)), epochs = 500, 
                        epsilon= 10^-8, seed=1337, variable_name){
  
 res <-  Adam_multi(y=y, X=X, etas = etas, B1s = B1s, B2s = B2s,theta = theta,
             epsilon = epsilon, epochs=epochs, seed=seed)
 plots <- plot_GD_multi(res, variable_name)
 return(list(res, plots))
}




## ----task_3_mains, results='hide', fig.show='hide'-----------------------------------------------------------------
etas_vec = c(0.001, 0.1, 1)
Beta_ones <- c(0.9, 0.85, 0.5)
Beta_twos <- c(0.999, 0.95, 0.8)
X <- model.matrix(admit ~ gre_sd + gpa_sd, binary)
y <- binary$admit

adam_results <- task_3_main(y=y, X=X, etas = etas_vec,
            B1s = Beta_ones, B2s = Beta_twos, variable_name = "gre_sd" )


## ----adam_results, fig.cap = "Results for Adam implementation", fig.width=8, fig.height=8--------------------------
plot(adam_results[[2]])


