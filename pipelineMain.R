
#############################################
##  Code to accompany the manuscript, ""Implementation of a dynamic model 
##  updating pipeline provides a systematic process for maintaining 
##  performance of prediction models"
##
##  Here we demonstrate the use of a proactive pipeline implemented
##   with Bayesian updating.
##
##  Helper functions are found in the file "pipeline HelperFxns.R"
##
##  August 2024
##
#############################################

## ----
## Load necessary libraries

library(survival)
library(dplyr)
library(rstan)

source("pipelineHelperFxns.R")

## ----
## Prepare for Bayesian updating

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
##Load STAN model file 
sm <- stan_model(file="pipeline_sampleSTAN.stan")


#############################################
## Step 0. Define some settings for the 
##         dynamic updating pipeline
##
##############################################

## We use split sample validation here for simplicity
##    but cross-validation or bootstrapping may be more
##    appropriate to control for optimism
## Define the %age of data that should be training vs holdout
TRAIN_PCT <- 0.75

## Prediction horizon is two years
predHor <- 730 ##(in days)

## We will update using Bayesian dynamic updating
##  with forgetting factors of 0.9 and 0.8
LAMBDA_Fa <- 0.9
LAMBDA_Fb <- 0.8

## For STAN:
ITER <- 5000  
CHAINS <- 4

## Formulas to be used in survival models
## First, the left hand side 
formLHS <- "Surv(eventtime,status)~"
## Second, the right-hand for Kaplan-Meier
formKM <- "1"  


## Specify non-inferiority margins (amount) and margin type 
##   margin type can be Percent==0, Value==1, or a threshold==2 that doesn't depend on performance in current data
##                        metric, amount, type
nonInfMargin <- data.frame("m1"= c(0.02, 0), ##Non-inferior if within 2%
                           "m2"= c(0.005, 1)) ##Non-inferior if within 0.005

## For the metrics above, is Higher or Lower better?                           
metricDirxn <- c("Higher", "Lower" )                         


#############################################
## Step 1. Develop the "original" model
## 
##    We will use the pbcseq data (survival library)
##    This is meant only to illustrate the pipeline,
##    not to be a proper analysis of the PBC data.
##
##    Note also that this is a small dataset and the sample
##    size is insufficient for proper validation
##    
##    Take the day = 0 observations to be the development dataset, P_0.
##    We will use the predictors trt and log(bili)
##    and treat both transplant and death as an event.
##    The outcome is 2-year survival, so we take 2-years of
##    follow-up from the data.
##
##    We will fit a Bayesian model with vague priors for
##    the original model, f_0
##
#############################################

## Get the data
data(pbc)
pbcseq$log_bili <- log(pbcseq$bili)
pbcseq$status <- ifelse(pbcseq$status>0, 1, 0)
pbcseq$eventtime <- pbcseq$futime
keep <- c("id", "eventtime", "status", "day", "trt", "log_bili" )
dat0 <- pbcseq[pbcseq$day==0, keep]
## censor those with futime > 2 years (i.e. pretend we only have 2 years of follow-up)
dat0$status <- if_else(dat0$eventtime > 730, 0, dat0$status)
dat0$futime <- if_else(dat0$eventtime > 730, 730, dat0$eventtime)
updateCovs <- c("trt", "log_bili")

## Split into train / holdout datasets
set.seed(357) 
size<-floor(TRAIN_PCT*nrow(dat0))
trainRows <- sample(seq(1,size), size=size, replace=FALSE)
qDatTrain <- dat0[trainRows,]
qDatHold <- dat0[-trainRows,]

## Original model developed using Bayesian methods
priors <- getPriOriginal(n=length(updateCovs))
prior_int <- list()
prior_int$location <- 0
prior_int$scale <- 5
## Format data for STAN
stanData <- getStanData( dat=qDatTrain, co=updateCovs )
fitStan <- sampling(sm, data=stanData, chains=CHAINS, iter=ITER)
mSumm <- as.data.frame(summary(fitStan)$summary)
Bdraws <- rstan::extract(fitStan)



#############################################
## Step 2. New data is collected 
##
##   We will use the pbcseq data that was collected over
##   the next 1 year to represent an update dataset
##
#############################################

## Get next period's data and keep only the most recent record for each id
dat1 <- pbcseq[pbcseq$day>1 & pbcseq$day<366 & pbcseq$eventtime>365,keep]
dat1 <- as.data.frame( dat1 %>% group_by(id) %>%
        slice_tail(n = 1) %>%  
        ungroup())

## censor those with eventtime > 3 years (i.e. pretend we only have 2 years of follow-up)
dat1$status <- if_else(dat1$eventtime > 1096, 0, dat1$status)
dat1$eventtime <- if_else(dat1$eventtime > 1096, 1096, dat1$eventtime)

updateCovs <- c("trt", "log_bili") 
## NOTE: If a new predictor became available, update the above to include
##    it and define a prior for the new variable

## split into train / holdout datasets
size<-floor(TRAIN_PCT*nrow(dat1))
trainRows <- sample(seq(1,size), size=size, replace=FALSE)
pDatTrain <- dat1[trainRows,]
pDatHold <- dat1[-trainRows,]


#############################################
## Step 3. Update the current model using the new data
## 
##  We illustrate with Bayesian updating
##     with two forgetting factors. 
##  Multiple update candidates can be generated with 
##     multiple forgetting factors and the best performer chosen
##
#############################################

## Update candidate 1
## Obtain priors based on the previous period's selected model
priors <- getPriUpdate( m_b = mSumm, lambdaF=LAMBDA_Fa)
## Format STAN data
stanData <- getStanData( dat=pDatTrain, co=updateCovs )
fitStan_a <- sampling(sm, data=stanData, chains=CHAINS, iter=ITER)
mSumm_1a <-  as.data.frame(summary(fitStan_a)$summary)
Bdraws_1a <- rstan::extract(fitStan_a)

## Update candidate 2
## Obtain priors based on the previous period's selected model
priors <- getPriUpdate( m_b = mSumm, lambdaF=LAMBDA_Fb)
## Format STAN data
stanData <- getStanData( dat=pDatTrain, co=updateCovs )
fitStan_b <- sampling(sm, data=stanData, chains=CHAINS, iter=ITER)
mSumm_1b <-  as.data.frame(summary(fitStan_b)$summary)
Bdraws_1b <- rstan::extract(fitStan_b)


#############################################
## Step 4. Evaluate the update candidate(s) and the previous
##         period's model in the new data
##         
#############################################

## Get a vector of performance metrics for the updated model, mSumm_1
##     and a vector of the same performance metrics for the previous
##     period's model
## Choose which performance metrics to implement
## Also choose how to compare them -- via hypothesis testing framework
##     or by comparing point estimates

## The implementation of this will vary by case

## Suppose we've calculated two metrics (m1, m2) for each model and the results are:
metrics_0 <- c(0.80, 0.1)
metrics_1a <- c(0.82, 0.12)
metrics_1b <- c(0.815, 0.1)

## For this simple example, we compare the point estimates but confidence
##  intervals could be constructed on the performance metrics (e.g. by bootstrapping)
##  and a hypothesis testing framework used to compare them

## Is the updated model superior on at least one criteria?
sup_1a <- isSuperior(metrics_0, metrics_1a)
sup_1b <- isSuperior(metrics_0, metrics_1b)

nonI_1a <-isNonInferior(metrics_0, metrics_1a)
nonI_1b <-isNonInferior(metrics_0, metrics_1b)

## In this example, the model with metrics_1a is superior to
##   last period's model on performance metric m1 BUT, based on
##   the non-inferiority margins defined, it is not non-inferior 
##   on all metrics.  Therefore, this model is NOT an acceptable update.

## The model with metrics_1b is superior to last period's model on
##   performance metrics m1 and it is non-inferior on all metrics.
##   This update is acceptable and will become the model for the next period.

## The procedure is repeated as new data is acquired.



