library(rstan)
library(bayesplot)
library(ggplot2)
library(coda)
library(circular)
library(here)
library(expm) 
library(msm)
library(fGarch)
library(splines2)
library(splines)
library(gridExtra)
library(ggpubr)


# This script was used to simulate all CAV data sets used
# The code is slightly adapted, but mostly corresponds to the script of the same 
# name from the original paper (Williams et al. (2020))

msm_res <- matrix(0, 22, 100)
models_msm  <- list()

for(n in 1:50){
  set.seed(n)
  Population_Study <- TRUE
  Dir <- "./New_work/Data/" # /NoBiasData/ and Population_Study <- FALSE when no delayed enrollment
  
  # Set the sample size.  Note that the cav data set has 622 subjects.
  N <- 700 # 1500 and n = 1 for extension results
  # Choose the discretization of time.
  dt <- 1/365
  
  
  # The years and iyears columns were both centered at round(mean(years),0) = 4 in the cav data set.  These are the true parameter values for the uncentered data ( intercept - coef*mean ).
  betaMat <- matrix(c(-2.54,  0.11, -0.56,
                      -2.94, -0.24,  0.15,
                      -1.10, -0.15, -0.03,
                      -3.92,  0.23,  0.21,
                      -2.12,  0.08,  1.17), nrow=5, byrow=TRUE)
  
  errorMat <- matrix(c(.99, .01,   0, 0,
                       .24, .70, .06, 0,
                       0, .11, .89, 0,
                       0,   0,   0, 1), nrow=4, byrow=TRUE)
  
  initProbs <- c( .95, .04, .01, 0)
  
  
  
  #----------------------------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------------------
  # Collect information about the real data set.
  #----------------------------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------------------
  
  
  
  ptnum <- unique(cav$PTNUM)
  N_cav <- length(ptnum)
  interObsTime <- NULL
  propMale <- 0
  propDeaths <- 0
  NumObs <- rep(0,N_cav)
  for(i in 1:N_cav){   
    subject <- cav[cav$PTNUM==ptnum[i],,drop=FALSE]
    
    # The number of observations for each subject.
    NumObs[i] <- nrow(subject)
    
    # The times between observations.
    if(!(4 %in% subject$state)){  interObsTime <- c( interObsTime, round( diff(subject$years), 3))  }
    
    # Determine whether the subject is male.
    propMale <- propMale + as.integer(subject$sex[1]==1)
    
    # Determine whether the subject's death was observed.
    propDeaths <- propDeaths + as.integer(4 %in% subject$state)
  }
  propMale <- propMale / N_cav
  propDeaths <- propDeaths / N_cav
  
  
  
  #----------------------------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------------------
  # Fill in the states using the transition rate matrix and error matrix estimated on the real data set.
  #----------------------------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------------------
  
  
  
  Q <- function(iyears,sex,betaMat){
    
    q1  = exp( c(1,iyears,sex) %*% betaMat[1,] )  # Transition from state 1 to state 2.
    q2  = exp( c(1,iyears,sex) %*% betaMat[2,] )  # Transition from state 2 to state 3.
    q3  = exp( c(1,iyears,sex) %*% betaMat[3,] )  # Transition from state 1 to death.
    q4  = exp( c(1,iyears,sex) %*% betaMat[4,] )  # Transition from state 2 to death.
    q5  = exp( c(1,iyears,sex) %*% betaMat[5,] )  # Transition from state 3 to death.
    
    qmat = matrix(c( 0,q1, 0,q2,
                     0, 0,q3,q4,
                     0, 0, 0,q5,
                     0, 0, 0, 0),nrow=4,byrow=TRUE) 
    diag(qmat) = -rowSums(qmat)
    
    return(qmat)
  }
  
  
  rawData <- NULL
  propDeaths_sim <- 0
  NumObs_sim <- NULL
  for(i in 1:N){
    
    # Sample the gender, as proportional to the cav data set.
    sex <- as.integer(runif(1,0,1) < propMale)
    
    # Sample for an initial state.
    trueState <- sample(1:4, size=1, prob=initProbs)
    
    # Sample the remaining states until death.
    years <- 0
    time1 <- 0
    s <- trueState
    while(s < 4){
      
      # Infinitesimal transition rates.
      qmat <- Q(floor(time1),sex,betaMat)
      
      # Possible next states.
      moveToStates <- which(qmat[s,] > 0)
      
      # Sample the wait times before transition to each of the next possible states.
      waitTimes <- rexp( n=length(moveToStates), rate= qmat[s,moveToStates])
      
      # If any of the wait times are smaller than dt, then transition to the state with the minimum wait time.
      min_waitTime <- min(waitTimes)
      if(min_waitTime < dt){  s <- moveToStates[ which(waitTimes == min_waitTime) ]  }
      
      time1 <- time1 + dt
      
      years <- c( years, time1)
      trueState <- c( trueState, s)	
    }
    timeOfDeath <- tail(years,1)
    
    # Sample inter-observation times from the cav data set.  Maximum of 20 years in study.
    visitTimes <- NULL
    # If this a population study, then first observations may occur after zero years of enrollment.
    if(Population_Study == 'TRUE'){  
      
      # first observation occurs after zero years, with probability .75.
      time2 <- (runif(1,0,1) < .75) * rnorm( n=1, mean=5, sd=1)  
      
    } else{  time2 <- 0  }
    
    while(time2 < min( 20, timeOfDeath)){
      
      visitTimes <- c( visitTimes, time2)
      time2 <- time2 + sample( interObsTime, size=1)
    }
    
    # If first visit time occurs after death, then subject is NOT entered into the study.
    if( !is.null(visitTimes) ){
      
      # If death occured before the study ended, then record the time of death.
      if( timeOfDeath < 20 ){  visitTimes <- c( visitTimes, timeOfDeath) }
      
      
      n_i <- length(visitTimes)
      state <- NULL
      for(k in 1:n_i){  state <- c( state, tail( trueState[ years <= visitTimes[k] ], 1))  }
      
      ptnum <- rep(i,n_i)
      years <- visitTimes
      rawData <- rbind( rawData, data.frame(ptnum,years,sex,state) )
      
      if(4 %in% state){  propDeaths_sim <- propDeaths_sim + 1  }
      NumObs_sim <- c( NumObs_sim, n_i)
    }
    #print(i)
  }
  
  colnames(rawData) <- c('ptnum','years','sex','state')
  N <- length(unique(rawData$ptnum))
  propDeaths_sim <- propDeaths_sim / N
  
  ## Use only to simulate F-V response variable
  # rawData$obs <- 0 
  # for(i in 1:nrow(rawData)){
  #   if(rawData$state[i] == 1){
  #     rawData$obs[i] <- rgamma(1, 6, 0.95)
  #   } else{
  #     if(rawData$state[i] == 2){
  #       rawData$obs[i] <- rgamma(1, 48, 4)
  #     } else{
  #       if(rawData$state[i] == 3){
  #         rawData$obs[i] <- rsnorm(1, 18, 3, 0.7)
  #       }
  #     }
  #   }
  # }
  # states <- rawData$state
  # # Add noise to the states.
  for(i in 1:nrow(rawData)){	rawData$state[i] <- sample(1:4, size=1, prob=errorMat[rawData$state[i],])  }

    
  #----------------------------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------------------
  # Add censored rows.
  #----------------------------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------------------
  
  
  
  obstype_hmm <- rep(1,nrow(rawData))
  obstype_hmm[rawData$state == 4] = 2
  obstype_hmm[!duplicated(rawData$ptnum)] = 3
  
  iyears <- floor(rawData$years)
  
  obstrue <- rep(0,nrow(rawData))
  
  hold <- cbind(rawData,obstrue,obstype_hmm,iyears)
  hold <- hold[,c('ptnum','years','iyears','sex','state','obstrue','obstype_hmm', "obs")]
  
  
  
  tempRow <- rep(0,ncol(hold))
  names(tempRow) <- c('ptnum','years','iyears','sex','state','obstrue','obstype_hmm', "obs")
  
  num <- 1
  cavData <- NULL
  for(i in unique(rawData$ptnum)){
    
    current <- NULL
    subject <- hold[hold$ptnum==i,,drop=FALSE]
    
    #------------------------------------
    #censoredAges <- unique( c( min(subject$age), ceiling(min(subject$age)):max(subject$age)) )
    censoredAges <- 0:max(subject$years)
    for(t in censoredAges ){
      
      # If 't' corresponds to an observed age, then the next row will include the observed clinical visit data.
      if(t %in% subject$years){	
        current <- rbind( current, subject[subject$iyears==floor(t),]) 
      } else{
        
        # Create a CENSORED row for each subject at each INTEGER year of years.
        tempRow['ptnum'] <- i
        tempRow['years'] <- t 
        tempRow['iyears'] <- t
        tempRow['sex'] <- subject$sex[1]
        tempRow['state'] <- 99
        tempRow['obstrue'] <- 1  
        tempRow['obstype_hmm'] <- 0
        
        current <- rbind( current, tempRow)
        
        # If 't' corresponds to an observed INTEGER years, then the subject was observed some time during this years.  According, the next row will include the observed clinical visit data.  Recall that integer years is simply the floor(years).
        if(t %in% subject$iyears){ current <- rbind( current, subject[subject$iyears==t,]) }
      }
      
    }
    #------------------------------------
    
    cavData <- rbind( cavData, current)
    #print(num)
    num <- num+1
  }
  colnames(cavData) <- c('ptnum','years','iyears','sex','state','obstrue','obstype_hmm', "obs")
  
  cavData$c_iyears <- cavData$iyears - floor(mean(cavData$years))
  
  save(cavData,file=paste0(Dir,'cavData', n, '.rda'))

}


