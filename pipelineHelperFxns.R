####################################
## Helper Functions





## Set vague priors of mean=0, scale=2.5
##  for development of original model
getPriOriginal<-function(n ){
  loc <- rep(0, n)
  sca <- rep(2.5, n)
  return(list(location=loc, scale=sca))
}


getPriUpdate<-function(m_b, lambdaF){
  n <- nrow(m_b) - 3
  loc <- m_b[1:n, "mean"]
  sca <- ((m_b[1:n, "sd"])^2)/lambdaF
  sca <- sqrt(sca)
  return(list(location=loc, scale=sca))
}

## Set up the data in the format we need it for STAN
## We need to work with censored and uncensored data separately
getStanData<-function(dat, co){
  X <- dat[,co]
  time <- dat$eventtime
  event <- dat$status
  N_censored <- length(which(dat$status==0))
  X_censored <- X[event==0,]
  X_uncensored <- X[event!=0,]
  
  ## Sort the priors as well 
  prior_mean_betas <- priors$location
  prior_sd_betas <- priors$scale
  
  N <- nrow(X)
  NC <- ncol(X_censored)
  stan_data <- list(N=N, 
                    N_uncensored=N-N_censored, 
                    N_censored = N_censored,
                    NC=NC,
                    X_censored = X_censored, 
                    X_uncensored=X_uncensored,
                    times_censored=time[event==0],
                    times_uncensored=time[event!=0],
                    prior_mean_mu = prior_int$location,
                    prior_sd_mu = prior_int$scale,
                    prior_mean_betas = prior_mean_betas,
                    prior_sd_betas = prior_sd_betas)
  return(stan_data)
}

getBayesPreds_R<-function(postDraws, dat, co){
  x_test <- dat[,co]
  x_test$eventtime <- dat$eventtime
  stanout <- apply(x_test, MARGIN=1, get_stanpreds, postDraws=postDraws, co=co)
  return(stanout)
}

get_stanpreds<-function(x_test, postDraws, co){
  out <- vector() #the lp, the expected, the surv prob
  x <- as.matrix(x_test[1:length(co)])
  lp <- postDraws$betas %*% (x)
  out[1] <- median(lp)
  
  ## for Weibull baseline hazard:
  H0t <- (x_test[length(x_test)]^postDraws$alpha) * exp(postDraws$intercept)
  
  expec <- as.numeric(exp(lp)) * as.numeric(H0t) 
  out[2] <- median(expec)
  
  lambdaPred <- lp+as.matrix(postDraws$intercept)
  # generate Weibull survival times
  weib_params <- cbind(lambdaPred, postDraws$alpha)
  survTs <- apply(weib_params, MARGIN=1, rweib_KT, n=1)
  out[3] <- length(which(survTs>predHor))/length(survTs)
  return(out)
}

## Function returns true if the model corresponding to met2
##   is superior to the model with met1 on at least one metric
isSuperior <- function(met1, met2){
  sup2 <- FALSE
  for(i in 1:length(met2)){
    if(metricDirxn[i] == "Higher"){
      sup2 <- ifelse(met2[i] > met1[i], TRUE, sup2)
    }else{
      sup2 <- ifelse(met2[i] < met1[i], TRUE, sup2)
    }
  }
  return(sup2)
}

## Function returns true if the model corresponding to met2
##   is non-inferior to the model with met1 on ALL metrics
isNonInferior <- function(met1, met2){
  nonI <- TRUE
  for(i in 1:length(met1)){
    if(nonInfMargin[i][2,1] == 0){ ##then percentage
      marg <- nonInfMargin[i][1,1]*met1[i]
    }else{ ## then value
      marg <- nonInfMargin[i][1,1]
    }
    if(metricDirxn[i] == "Higher"){
      nonI <- ifelse(met2[i] < met1[i] - marg , FALSE, nonI)
    }else{
      nonI <- ifelse(met2[i] > met1[i] + marg, FALSE, nonI)
    }
  }
  return(nonI)
}

