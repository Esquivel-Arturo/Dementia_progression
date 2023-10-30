# A rewrite of the routine from Jonathan.  Changes are
#   sensible line widths for the code 
#   use the same argument names as mcmc0
#   cmap instead of qmat, qcoef, 
#   return a named list
# His key idea of clustered updates remains

mcmc2 <- function(par, fn, iter, burnin=1000, prior.std, 
                  prop.std = .3, cmap, prior.means, reset=100, print=0) {
    # par contains the initial coefficient values, which are also the
    #  prior means if prior.means is missing.
    npar <- length(par)
    if (missing(prior.means)) prior.means <- par
    if (length(prior.means) != npar) stop("wrong length for prior.means")

    if (missing(cmap)) stop("cmap argument is required")
    if (max(cmap) != npar) 
        stop("cmap does not agree with the number of parameters")

    if (missing(prior.std)) stop("prior.std argument is required")
    if (length(prior.std) == 1) prior.std <- rep(prior.std, npar)
    if (length(prior.std) != npar) stop("wrong length for prior.std")

    if (print > 0) {
        pfile <- file("mcmctrace", "w")
        cat("burn in\n", file=pfile)
    }

    # Group the parameters into sets.  Each column of cmap is a linear
    #  predictor.  The default grouping of parameters is by linear predictor.
    #  Some parameters may be in more than one linear predictor, in  which
    #  case that set of linear predictors is collapsed.  
    pset <- 1:ncol(cmap)   # the default sets, one per linear predictor
    for (i in 1:npar) {
        j <- col(cmap)[cmap==i]   # all columns that have parameter i
        k <- which(pset %in% j)   
        pset[k] <- min(pset[k])   # all those columns become a set
    }
    pset <- match(pset, sort(unique(pset)))  #remove holes in the sequence
    pgroup <- integer(npar)  # which group does each parameter belong to
    pgroup[cmap] <- (pset[col(cmap)])[cmap > 0]  

    # Each parameter is now marked as being in group "pgroup"
    #  There will be a separate covariance for each group
    ngroup <- max(pgroup)
    ng2    <- c(table(pgroup)) #number per parameter group
    sigma <- tapply(1:npar, pgroup, 
                    function(i) diag(prior.std[i], nrow=length(i)))
    pindx <- split(1:npar, factor(pgroup))

    tau <- rep(prop.std, ngroup)  #this is used to tune the proposals
    beta <- matrix(0., iter, npar)
    
    zeros <- lapply(ng2, function(i) rep(0.0, i))
    pvar  <- diag(prior.std^2)

    # burn in
    bpar <- par
    prior.log <- fn(par) + dmvnorm(par, prior.means, pvar, log=TRUE)
    accept <- double(ngroup)
    for (i in 1:burnin) {
        logunif <- log(runif(ngroup, 0, 1))
        for (j in 1:ngroup) {
            current <- bpar
            index <- pindx[[j]]
            current[index] <- current[index] + rmvnorm(1, zeros[[j]],
                                                 tau[j]* sigma[[j]])
            logPosterior = fn(current)+ dmvnorm(beta, prior.means, pvar, log=TRUE)
            # Accept/reject proposed parameters.
            logratio = logPosterior - prior.log
            if(!is.finite(logratio))
                warning("infinite Metropolis-Hastings ratio occurred") 
            if (logunif[j] < logratio) { 
                # accept the proposed parameter
                prior.log <- logPosterior
                bpar <- current
                accept[j] <- accept[j] +1/reset
            }
        }
        if (i%%print ==0) 
            cat("i=", i, "par=", format(current), "\n", file=pfile)

        if (i%%reset ==0) {
            # Tune the proposal variance
            k <- seq(round(i/2), i)
            for (j in 1:ngroup) 
                sigma[[j]] <- cov(beta[k, pindx[j]])
            
            # tune up the proposal std to get a better acceptance rate
            #   model: acceptance rate approx = exp(-k*sd)
            #   solve for k* newsd = 0.9, which gives 40% acceptance
            # If the current sd is just right we will tend to oscillate:
            #   a phat < .4 decreases the sd too much, etc.  So damp it
            #   a little.  This also stops explosion if phat=0 or 1
            #   
            phat <- (.4 + 2* accept)/3
            k <- -log(phat)/tau
            tau <- .92/k 
            accept[] <- 0
            if (print > 0) {
                temp <- format(rbind(accept, k, prop.std), digits=2)
                cat(" accept  ", temp[1,], "\n", 
                    "k       ", temp[2,], "\n",
                    "prop.std", temp[3,], "\n", file=pfile)
            }
        }
    }

        
    # Now for the iterations in earnest
    # No more updates of the proposals
    if (print>0) cat("iteration\n", file=pfile)
    beta[1,] <- bpar
    loglik <- matrix(0, iter, ngroup)
    accept <- double(ngroup)
    for (i in 1:iter) {
        logunif <- log(runif(ngroup, 0, 1))
        for (j in 1:ngroup) {
            current <- beta[i,]
            index <- pindx[[j]]
            current[index] <- current[index] + rmvnorm(1, zeros[[j]],
                                                 tau[j]* sigma[[j]])
            logPost = fn(current)+ dmvnorm(beta, prior.means, pvar, log=TRUE)
            # Accept/reject proposed parameters.
            logratio = logPost - prior.log
            if(!is.finite(logratio))
                warning("infinite Metropolis-Hastings ratio occurred") 
            if (logunif[j] < logratio) { 
                # accept the proposed parameter
                prior.log <- logPost
                beta[i,] <- current
                accept[j] <- accept[j] +1
            }
            loglik[i,j] <- prior.log
        }
        
        if (i%%print ==0) 
            cat("i=", i, "par=", format(current), "\n", file=pfile)
    }
    
    list(par= beta, loglik=loglik, accept=accept/iter, groups = pgroup,
         corlist = sigma)
}
