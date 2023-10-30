predict.hmm <- function(object, newdata, type=c("link", "rate", "prevalence")) {
    if (!inherits(object, "hmm")) stop("only valid for hmm objects")
    type <- match.arg(type)
    if (type=="prevalence") Terms <- terms(object)
    else Terms <- delete.response(terms(object))

    if (missing(newdata) || is.null(newdata)) {
        X <- object$x
        if (is.null(X)) stop("program not yet finished")
    }
    else {
        mf <- model.frame(Terms, newdata, na.action=na.pass,
                          xlev=object$xlevels)
        X <- model.matrix(Terms, mf, contrasts.arg=object$contrasts)
    }

    eta <- X %*% coef(object, type="matrix")
    if (type == "link") eta
    else if (type=="rate") exp(eta)
    else {
        # For a prevalence curve the data is restricted to be a single
        #  set of sequential times, along with its set of time dependent
        #  covariates.
        if (!missing(newdata) || is.null(object$y)) Y <- model.response(mf)
        if (is.Surv(Y)) {
            if (attr(Y, "type") != "counting")
                stop("Y must have start and stop times")
            n <- nrow(mf)
            tstart <- Y[,1]
            tstop  <- Y[,2]
            if (any(tstart[-1] != tstop[-n]) || any(diff(tstop) <=0))
                stop("(start, stop) times must form an ordered sequence with no breaks")    
            ytime <- Y[,2]-Y[,1]
        }
        else ytime = diff(Y[,1])

        # The initial prevalence is assumed to hold at the smallest Y time
        nstate <- object$nstate
        prev <- matrix(0., nrow= length(ytime), ncol= nstate)
        prev[1,] <- object$istate
        rmat <- matrix(0.0, nstate, nstate)
        qmap <- which(object$qmatrix != 0)
        for (i in 1:n) { 
            rmat[qmap] <- eta[i,]
            diag(rmat) <- diag(rmat) - rowSums(rmat)
            prev[i+1,] <- prev[i,] %*% expm((tstart[i] - tstop[i]) * rmat)
        }
        dimnames(prev) <- list(c(tstart[1], tstop), dimnames(object$qmatrix)[[1]])
        prev
    }
}

            
            
      
