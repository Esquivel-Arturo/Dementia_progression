# Automatically generated from the hmmcode
hmmscore <- function(par, fn, gr, iter=30,
                     eps=1e-5, scale=4, shrink=1, debug=FALSE,
                     method=1, ...) {
    if (scale <2) stop("invalid value for scale")
    if (iter < 1) stop("iteration count must be 1 or more")

    ncoef <- length(par)
    
    # Some tracing information: the coefficients, loglik, and LM scale
    #  as we go along.  The second 2 are a matrix.
    cmat <- matrix(0, iter+1, ncoef)
    cmat[1,] <- par
    logmat <- matrix(0, iter+1, 4,
               dimnames=list(0:iter, c("loglik", "LM scale", "lambda", "step")))
    lm <- 4  # starting point
    
    # initial step
    fit <- gr(par)
    logmat[1,] <- c(fit$loglik, 0, 0,0)
    if (!is.list(fit) || is.null(fit$S))
        stop("hmmscore needs to call hmmboth")
    npar <- length(par)

    #  A simple solve() will blow up if the S +dmat is too close to
    #   singular; use a tempered svd that ignores the near zeros
    #  Solve using a generalized inverse -- the algorithm is essentially that
    #   of the ginv function in MASS
    pfun <- function(coef, deriv, S, dmat, shrink) {
        smat <- svd(S + dmat, nv=0)  #don't need V, since it is symmetric
        dpos <- (smat$d > max(smat$d[1]*eps, 0))
        dd <- ifelse(dpos, 1/smat$d, 0)
        # all the parentheses save a tiny bit of time 
        if (all(dpos)) x <- drop(smat$u %*% (dd*(t(smat$u) %*% deriv)))
        else if (!any(dpos)) stop("zero hessian in update") # impossible I think
        else x <-drop(smat$u[,dpos] %*%(dd[dpos] * (t(smat$u[,dpos, drop=FALSE]) %*% deriv)))
 
        coef + shrink * x
    }

    # The usual LM would start at diag(S *.001), we are more conservative
    lm <- 1/4
    log0 <- fit$loglik  # the initial loglik
    for (i in 1:iter) {
        switch(method,
            {S <- fit$S; dmat <- diag(diag(S)*lm)},
            {S <- fit$S2; dmat <- diag(diag(S)*lm)},
            {S <- fit$S2; dmat <- diag(fit$deriv^2 *lm)} )
        newpar <- pfun(par, fit$deriv, S, dmat, shrink)
        cmat[i+1,] <- newpar
        newfit <- gr(newpar)
        
        if (newfit$loglik <= 2*log0) {
            # fn() returns 2*intial -1 when the computation fails due
            #  to overflow/underflow; usually horrible parameters
            logmat[i+1,] <- c(newfit$loglik, lm, NA, 3) #bombed
            lm <- lm * scale * scale
        }
        else {
            # Use a quadratic approx to find a (possibly) improved guess
            # Equation 9.7.11 of Numerical Recipes in C, 1992
            gder <- sum(fit$deriv * (newpar -par))  # directional deriv
            lambda <- 0.5*gder/(fit$loglik + gder - newfit$loglik)
            logmat[i+1,] <- c(newfit$loglik, lm, lambda, 1) #standard step
 
            if (abs((newfit$loglik - fit$loglik)/fit$loglik) < eps) {
                trdata <- list(coef=cmat[1:(i+1),], loglik=logmat[1:(i+1),])
                if (newfit$loglik > fit$loglik) # very last step was better
                    return(c(newfit, list(coef=newpar, iter=i, converged=TRUE,
                                          trace=trdata)))  #all done
                else {
                    logmat[i+1, 4] <- 0  # no step
                    return(c(fit, list(coef=par, iter=i, converged=TRUE,
                                       trace=trdata)))
                }
            }
            
            if (newfit$loglik < fit$loglik) {
                lm <- lm*scale   # increase the penalty
                logmat[i+1, 4] <- 0  # no progress
                if (debug > 1) browser()
                
                if (gder >0) { # backtrack
                    try2 <- (1-lambda)*par + lambda*newpar
                    fit2 <- gr(try2)  # we expect this to usually work
                    if (fit2$loglik > newfit$loglik) { # shrinkage worked
                        logmat[i+1,4] <- fit2$loglik
                        fit <- fit2
                        par <- try2
                    }
                }    
            }
            else { #successful step
                fit <- newfit
                par <- newpar
                lm <- lm/scale
            }
            if (debug) cat("iter", i, "loglik", logmat[i+1,], "\n")
        }
    }
    trdata <- list(coef=cmat, loglik=logmat)
    c(fit, list(coef=par, iter= iter, converged = FALSE, trace=trdata))
}        
