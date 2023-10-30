# This is the routine, as received from Jonathan
library(MASS)
library(mvtnorm)

		 
mcmcRoutine <- function(par, fn, qmatrix, qcoef, rcoef, prior_means, prior_stddev, steps, burnin){
	
	# I'm assuming that 'par' inputs the initial coefficient values.
	
	
	
	# Return the 'qmatrix' indices corresponding to the estimable transition rates.  In order to update the coefficient
	# vector for each transition rate separately, I need to know what indices of 'par' correspond to each transition.
	# That task is done here --------------------------------------------------------------------------------------------
	TransIndices = which(qmatrix > 0, arr.ind=TRUE)
	nTrans = dim(TransIndices)[1] 
	par_TransIndices = list()
	for(l in 1:nTrans){
		# Which rows of 'qcoef' correspond to a given transition rate?  These correspond to the indices of 'par'.
		par_TransIndices[[l]] = which(qcoef[,'state1'] == TransIndices[l,'row'] & qcoef[,'state2']==TransIndices[l,'col']) 
	}
	# Now account for the indices of 'par' corresponding to 'rcoef'.
	par_TransIndices[[nTrans+1]] = (dim(qcoef)[1]+1):( dim(qcoef)[1] + dim(rcoef)[1] )
	#--------------------------------------------------------------------------------------------------------------------



	betaList = list()  # This list stores the MCMC parameter vector chain for each transition.  
	                   # Each transition has it's own matrix whose rows represent the sample chain.
	for(l in 1:(nTrans+1)){ betaList[[l]] = matrix(0,steps,length(par_TransIndices[[l]])) }

	Sigma = list()  # This list stores the proposal covariance matrices for each transition.
	for(l in 1:(nTrans+1)){ Sigma[[l]] = diag( rep(1,length(par_TransIndices[[l]])) ) }

	tau = rep(0.05,nTrans+1) # This vector is used to tune the proposals during the burnin period.
		
	prev_logPosteriorDens = -Inf  # Initialize.
	proposal = par # Initialize.

	### The algorithm...
	accept = rep(0,nTrans+1)
	for(ttt in 2:steps){

		### Update each of the beta vectors, in turn.
		for(l in 1:(nTrans+1)){ 
			
			# The indices of 'par' corresponding to transition 'l'.
			tempInd = par_TransIndices[[l]]

			# Propose a beta vector as in random walk MH.
			proposal[tempInd] = par[tempInd] + rmvnorm(1,rep(0,length(tempInd)),tau[l]*Sigma[[l]])  

			# Compute the log of the posterior density of the proposed parameter vector.
			logPosteriorDens =  fn(proposal) + dmvnorm(proposal, prior_means, diag(prior_stddev^2), log=TRUE)

			# Accept/reject proposed parameters.
			logMH_Ratio = logPosteriorDens - prev_logPosteriorDens

			if(!is.finite(logMH_Ratio)){ print('Warning: infinite Metropolis-Hastings ratio occurred.') }
			if( log(runif(1,0,1)) < logMH_Ratio ){ # If TRUE, then accept the proposed parameter.
				
				# The updates if the proposal is accepted.
				prev_logPosteriorDens = logPosteriorDens
				par = proposal
				betaList[[l]][ttt,] = proposal[tempInd]
				accept[l] = accept[l] + 1 
		
			} else{ 
			
				betaList[[l]][ttt,] = par[tempInd] 
				proposal = par			
			
			}
		
			# During the burnin period, update the proposal covariance in each step to capture 
			# the relationships within the parameters vectors for each transition.  This helps
			# with mixing.
			if(30 < ttt && ttt < burnin){ Sigma[[l]] = cov(betaList[[l]]) }
			# Tune the proposal covariance every 50 steps for each transition to achieve
			# reasonable acceptance ratios.
			if(ttt %% 50 == 0 && ttt < burnin){ 
				if(accept[l]/ttt < 0.5){ 
					    tau[l] = max(tau[l] - 0.005,0.0001) # Make smaller proposals if acceptance ratio is too low.
				} else{ tau[l] = tau[l] + 0.005 } # Make larger proposals if acceptance ratio is too high.
			}

			#cat('transNum',l,'\n')
		}
	
		print(ttt)
	}
	
	accept = accept / steps
	print(accept)
	
	list(accept, TransIndices, betaList)
}
