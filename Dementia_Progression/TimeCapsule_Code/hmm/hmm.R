# Automatically generated from the hmmcode
hmm <- function(formula, data, subset, weights, na.action= na.pass, 
                id, istate, otype, qmatrix, qcoef, rfun, rcoef, 
                pfun= hmminit, pcoef, entry, death,
                mfun=hmmtest, mpar= list(), 
                mc.cores= getOption("mc.cores", 2L),
                icoef, scale=c(TRUE, FALSE), penalty, constraint,
                debug=0, makefork=FALSE) {
    Call <- match.call()
    time0 <- proc.time()
    indx <- match(c("formula", "data", "subset", "weights",
                    "id", "istate", "otype"), names(Call), nomatch=0)
    if (indx[1] ==0) stop("a formula argument is required")
    if (indx[5] ==0) stop("an id argument is required")
    temp <- Call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")
    temp$na.action <- na.pass   #deal with NA later
    mf <- eval(temp, parent.frame())

    if (nrow(mf) ==0) stop("data has 0 rows")
    Terms <- terms(mf)
    attr(Terms, "intercept") <- 1  # ignore any "-1" in formula
    termnames <- c("(Intercept)", attr(Terms, 'term.labels'))

    id <- model.extract(mf, "id")
    if (length(id)==0) stop("an id variable is required")
    otype <- model.extract(mf, "otype")
    if (is.null(otype)) otype <- rep(1L, nrow(mf))
    if (any(otype==2) & missing(death))
        stop("death argument is missing and the data has deaths")
    first <- which(!duplicated(id))

    # missing values are done in a standard way if Y is a survival object
    # the do.call() construct allows na.action to be a function or a string
    if (missing(na.action)) na.action <- options("na.action")
    Y <- model.response(mf)
    X <- model.matrix(Terms, mf)
    xassign <- attr(X, "assign")
    if (inherits(Y, "Surv")) {
        stop("the survival portion off hmm is not longer supported")
        temp <- do.call(na.action, cbind(Y, X))
        if (nrow(temp) != nrow(Y)) { # some were tossed
            na.action <- attr(temp, "na.action")
            Y <- Y[-temp, ,drop=FALSE]
            X <- X[-temp, ,drop=FALSE]
            id <- id[-temp]
        }
        else na.action <- NULL
        if (attr(Y, "type") != "mcounting")
            stop("survival must be a multi-state counting process")
        ylevels <- (attributes(Y))$states
        
        yobs  <- Y[,3, drop=FALSE]
        otype <- otype * pmax(yobs, 1)
        itime <- Y[first,1]     # do I ever use this?
        ytime <- Y[,2] - Y[,1]  # length of each interval
        censor <- (yobs==0)
        ny <- 1
        
        istate <- model.extract(mf, "istate")
        if (missing(istate)) 
            stop("for a survival type endpoint, istate is required")
        else istate <- matrix(istate)
    }
    else {
        # the data set will have a row for each visit. 
        if (!is.matrix(Y) || ncol(Y) < 2) 
            stop("response must have at least two columns")
        ny <- ncol(Y) -1
        if (ny != length(rfun))
            stop("the number of mapping functions must be equal to the number of responses")
        last  <- !duplicated(id, fromLast=TRUE)
                                      
        if (is.character(Y)) {   # the caller used cbind(), fix this first
            ymiss <- is.na(Y)
            ylevels <- lapply(2:ncol(Y), function (i) {
                temp <- suppressWarnings(as.numeric(Y[,i]))
                if (any(is.na(temp) & !ymiss[,i])) 
                    levels(as.factor(Y[,i])) else NULL
                })
            tempy <- Y
            Y <- matrix(0, nrow(Y), ncol(Y))
            Y[,1] <- suppressWarnings(as.numeric(tempy[,1]))
            for (i in 2:ncol(Y)) {
                if (is.null(ylevels[[i-1]])) 
                    Y[,i] <- suppressWarnings(as.numeric(tempy[,i]))
                else Y[,i] <- as.numeric(as.factor(tempy[,i]))
            }
        }
        else ylevels <- attr(Y, "ylevels") # if they used hbind
                                      
        not.needed <- (last & otype %in% c(0, 2))
        X[not.needed,] <- 0   #turn any missings into zeros, and keep the row
        
        toss <- which(is.na(Y[,1]) | apply(is.na(X), 1, any))
        if (any(toss)) {
            na.action <- toss
            class(na.action) <- "omit"
            Y <- Y[-toss,, drop=FALSE]
            X    <- X[-toss,, drop=FALSE]
            id   <- id[-toss]
            otype <- otype[-toss]
            first <- which(!duplicated(id))
            last <- !duplicated(id, fromLast=TRUE)
        }
        else na.action <- NULL
        ytime <- c(diff(Y[,1]), 0)  # time to next obs, ignored for last
        yobs  <- Y[,-1, drop=FALSE] # the observed Y values
    }
    ymiss <- is.na(yobs)  #some individual elements may be missing
    weights <- model.weights(mf)
    if (length(weights) >0) stop("weights are not supported")
    # standardize the X matrix
    xmean  <- rep(0.0, ncol(X))
    xscale <- rep(1.0, ncol(X))
    keep <- which(!(otype==2 | (last & otype==0)))
    if (is.logical(scale)) {
        if (length(scale)==1) scale <- c(scale, scale)
        if (scale[1]) { # recenter the data
            xmean <- colMeans(X[keep,,drop=FALSE])
            for (i in 1:ncol(X)) if (all(X[keep,i] == X[keep[1],i])) xmean[i] <- 0
            }
        if (scale[2]) {
            xscale <- apply(X[keep,,drop=FALSE], 2, sd)
            xscale <- ifelse(xscale < 1e-8, 1, xscale)  #avoid division by 0
            }
        scale <- list(center=xmean, scale=xscale)
    }
    else if (is.list(scale)) {
        if (length(scale)!=2 || !all(sapply(scale, is.numeric)))
            stop("scale must be logical, or a list containing the center and scale")
        if (any(sapply(scale, length) != ncol(X)))
            stop("wrong length for numeric components of scale")
        xmean <- scale[[1]]
        xscale <- scale[[2]]
        }
                       
    if (any(xmean!=0 | xscale !=1))
        X <- scale(X, center=xmean, scale=xscale)
    temp <- diff(match(id, id))
    if (any(temp < 0)) {
        indx <- 1 + min(which(temp<0))
        stop("all the rows for each subject must be contiguous in the data set",
             paste("(row", indx, ")"))
    }

    if (any(ytime[!last] <=0)) {
        temp <- seq(along=ytime)[!last]
        indx <- min(temp[ytime[!last] <=0])
        stop("the rows for each subject must be in time order, with no duplicate times", paste("(row", indx, ")"))
    }
    if (missing(qmatrix)) stop("the qmatrix argument is required")
    if (!is.matrix(qmatrix) || (nrow(qmatrix) != ncol(qmatrix))) 
        stop("qmatrix must be a square matrix")
    nstate <- nrow(qmatrix)
    temp <- dimnames(qmatrix)
    if (length(temp[[1]]) > 0) {
        statenames <- temp[[1]]
        if (length(temp[[2]]) > 0 && temp[[2]] != temp[[1]])
            stop("row and column names for qmatrix must be identical, if present")
    }
    else if (length(temp[[2]]) >0) statenames <- temp[[2]]
    else statenames <- 1:nstate

    diag(qmatrix) <- 0
    if (any(qmatrix < 0)) stop("qmatrix elements must be >=0")
    qmap <- which(qmatrix != 0)   
    ntransitions <- length(qmap)
    if (all(qmatrix[row(qmatrix) > col(qmatrix)] == 0)) uppertri <- TRUE
    else {
        uppertri <- FALSE
        # the line below will eventually go away
        stop("transition matrix must be upper triangular")
        }
    qcoef2 <- data.frame(state1 = row(qmatrix)[qmatrix>0],
                         state2 = col(qmatrix)[qmatrix>0],
                         term   = 0,    #intercept
                         coef  = 1:ntransitions + 100,
                         init   = log(qmatrix[qmatrix>0]),
                         lp     = 1:ntransitions)
    qtest <- function(x, allowed, e1, e2, label="qcoef") {
         if (is.numeric(x)) {
            if (any(x != floor(x)) || any(x < 1))
                stop(paste(label, "numeric", e2, "must be integers greater than 0"))
            if (any(x > length(allowed))) 
                stop(label, ": numeric", e1, "that is > number of", e2)
            x
        }
        else {
            temp <- match(x, allowed, nomatch=0)
            if (any(temp==0)) stop(paste(label, ": unrecognized", e1, "name"))
            temp
        }
    } 

    has.qcoef <- !(missing(qcoef) || is.null(qcoef))
    has.rcoef <- !(missing(rcoef) || is.null(rcoef))
    has.pcoef <- !(missing(pcoef) || is.null(pcoef))

    if (has.qcoef) {
        if (!is.data.frame(qcoef)) stop("qcoef must be a data frame")
        index <- match(c("state1", "state2", "term", "coef"),
                       names(qcoef), nomatch=0)
        if (any(index==0)) 
            stop("qcoef must contain variables named state1, state2, term, and coef")
        qcoef$state1 <- qtest(qcoef$state1, statenames, "state", "states")
        qcoef$state2 <- qtest(qcoef$state2, statenames, "state", "states")
        itemp <- as.matrix(qcoef[, index[1:2]])
        if (any(qmatrix[itemp] ==0)) 
            stop("qcoef contains an invalid state1 to state2 transition")
        qmatrix[qmatrix >0] <- 1:ntransitions  # the transition number
        qcoef$lp <- qmatrix[itemp]
        if (is.null(qcoef$init)) qcoef$init <- 0.0
        if (is.character(qcoef$term)) qcoef2$term <- "(Intercept)"
        qcoef <- rbind(qcoef, qcoef2)  # initial values from the qmatrix
        qcoef <- qcoef[!duplicated(qcoef[,1:3]),]
    }
    else qcoef <- qcoef2
    if (has.rcoef) {
        if (!is.data.frame(rcoef)) stop("rcoef must be a data frame")
        index <- match(c("response", "lp", "term", "coef"), 
                       names(rcoef), nomatch=0)
        if (ny > 1 && any(index==0) )
            stop("rcoef must contain variables named response, lp, term, and coef")
        else if (ny ==1 && any(index[-1] ==0))
            stop("rcoef must contain variables named lp, term, and coef")
        if (!is.null(rcoef$response) && !all(rcoef$response %in% 1:ny))
            stop("rcoef$response is out of range")
        if (any(is.na(rcoef$lp))) stop("missing lp value in rcoef")
        if (ny >1) {  # renumber the linear predictors so as to be unique
            maxlp <- max(rcoef$lp)
            newlp <- rcoef$lp + rcoef$response * maxlp
            newlp <- match(newlp, sort(unique(newlp)))
            rcoef$lp <- newlp
            b2map <- vector("list", ny)
            for (i in 1:ny) b2map[[i]] <- unique(newlp[rcoef$response==i])
        }
        else {
            rcoef$lp <- match(rcoef$lp, sort(unique(rcoef$lp)))
            b2map <- list(sort(unique(rcoef$lp)))
        }
    }
    if (has.pcoef) {
        if (!is.data.frame(pcoef)) stop("pcoef must be a data frame")
        index <- match(c("lp", "term", "coef"), 
                       names(pcoef), nomatch=0)
        if (any(index==0)) 
            stop("pcoef must contain variables named lp, term, and coef")
        if (any(is.na(pcoef$lp))) stop("missing lp value in pcoef")
    }
    termvars <- table(xassign)  #max coefficients for a term
    testterm <- function(tcoef, tlab) {
        if (is.numeric(tcoef$term)) {
            if (any(tcoef$term != floor(tcoef$term)) || any(tcoef$term < 0))
                stop(paste(tlab, ": numeric terms must be integers greater >= 0"))
            if (any(tcoef$term > length(termnames)-1)) 
                stop(tlab, ": numeric term that is > number of terms in the model")
        }
        else {
            temp <- match(tcoef$term, termnames, nomatch=0)
            if (any(temp==0)) stop(tlab, ": unrecognized term name")
            tcoef$term <- temp-1
        }
            
        if (!is.numeric(tcoef$coef) || any(tcoef$coef != floor(tcoef$coef)))
            stop(tlab, ": param variable must contain non-negative integers")
        
        # Any intial value columns.  There may be none
        init <- matrix(0, nrow=nrow(tcoef), ncol=max(termvars))
        temp <- as.matrix(tcoef[,grepl("init", names(tcoef))])
        if (ncol(temp) > 0) { #found something
            init[,1:ncol(temp)] <- temp
            init <- ifelse(is.na(init), 0, init) # any not set become 0
        }
            
        tcoef$lp  <- match(tcoef$lp, sort(unique(tcoef$lp)))
        
        beta <- matrix(0., nrow=ncol(X), ncol=max(tcoef$lp))
        cmap <- matrix(0L, nrow(beta), ncol(beta))
        ucoef <- unique(tcoef$coef)
        ucoef <- sort(ucoef[ucoef>0])  # user may not have used 1, 2, 3,..
        # It is not legal to have two rows with the same coef param that
        #  point to different terms.
        # cpar = parameter numbers for each coef param
        cpar <- vector("list", length(ucoef))
        k <- 0
        for (i in seq_along(ucoef)) {
            j <- which(tcoef$coef == ucoef[i])
            if (any(tcoef$term[j] != tcoef$term[j[1]]))
                stop(tlab, ": the same coefficient param points to two terms")
            indx <- 1:termvars[tcoef$term[j[1]] +1] +k
            cpar[[i]] <- indx 
            k <- max(indx)
        }
        for (i in 1:nrow(tcoef)) {
            j <- (xassign == tcoef$term[i]) #columns of X for this term
            k <- 1:sum(j)  #how many coefs for this term
            if (any(init[i, k] !=0)) beta[j, tcoef$lp[i]] <- init[i, k]
            if (tcoef$coef[i] >0) # fixed coefs don't appear in cmap
               cmap[j, tcoef$lp[i]] <- cpar[[match(tcoef$coef[i], ucoef)]]
        }
        list(beta=beta, cmap=cmap)
    }
         
    temp <- testterm(qcoef, "qcoef")
    beta <- temp$beta
    cmap <- temp$cmap
    bcount <- c(q = ncol(cmap), m=0, p=0) # number of columns of beta and cmap
    pcount <- c(length(unique(cmap[cmap>0])), 0L, 0L) #number of parameters
    if (has.rcoef) {
        temp <- testterm(rcoef, "rcoef")
        bcount[2] <- ncol(temp$cmap)
        pcount[2] <- length(unique(temp$cmap[temp$cmap>0]))
        beta <- cbind(beta, temp$beta)
        cmap <- cbind(cmap, temp$cmap + ifelse(temp$cmap>0, max(cmap), 0))
    } 
    if (has.pcoef) {
        temp <- testterm(pcoef, "rcoef")
        bcount[3] <- ncol(temp$cmap)
        pcount[3] <- length(unique(temp$cmap[temp$cmap>0]))
        beta <- cbind(beta, temp$beta)
        cmap <- cbind(cmap, temp$cmap + ifelse(temp$cmap>0, max(cmap), 0))
    }
    b1 <- 1:bcount[1]
    b2 <- seq(bcount[1]+1, length=bcount[2])  #might be nothing
    b3 <- seq(bcount[1] + bcount[2] +1, length=bcount[3])
    nparm <- max(cmap)
    beta <- beta * scale$scale
    if (!missing(icoef)) {
        if (is.matrix(icoef)) {
            bcol <- paste0(row(qmatrix)[qmatrix!=0], ":",
                           col(qmatrix)[qmatrix!=0])
            if (bcount[2]>0) bcol <- c(bcol, paste0("R", 1:bcount[2]))
            if (bcount[3]>0) bcol <- c(bcol, paste0("p", 1:bcount[3]))
            
            dname <- dimnames(icoef)
            cindex <- match(dname[[2]], bcol)
            rindex <- match(dname[[2]], dimnames(X)[[2]])
            if (any(is.na(cindex)) || any(is.na(rindex)))
                stop("row and column names for icoef don't match the model")
            beta[rindex, cindex] <- icoef
        }
        else {
            if (length(icoef) != nparm) 
                stop("icoef is the wrong length")
            beta[cmap>0] <- icoef[cmap]
        }
    }    

    param <- double(nparm)
    temp <- (cmap>0)
    param[cmap[temp]] <- beta[temp]   #load it with the initial parameters
    rindex <- which(qmatrix > 0)
    tempfun <- function(x) length(unique(x[x>0]))
    parmcount <- c(tempfun(cmap[,b1]), tempfun(cmap[,b2]), tempfun(cmap[,b3]))
    if (ny ==1 && !is.list(rfun)) rfun <- list(rfun)  # make it a list of length 1
    if (ny != length(rfun))
        stop("must have one modeling function per respsonse")

    if (bcount[2]) b2map <- lapply(b2map, function(x) x+ bcount[1])
    else b2map <- vector("list", ny)

    doresponse <- (otype==1 | otype==3)  # the rows for a response function
    if (!any(doresponse))
        stop("all observations are exact or censored")
    # Do a dummy call to the response function(s), and make sure they
    #  return an object of the right shape.
    eta <- X%*% beta
    for (i in 1:ny) {
        keep <- (doresponse & !is.na(yobs[,i]))
        fit <- try(rfun[[i]](yobs[keep,i], nstate, 
                             eta[keep,b2map[[i]], drop=FALSE], gradient=FALSE))
        if (is.character(fit)) 
            stop("call failed for response function ", i)
        if (length(i) ==1) fit <- as.matrix(fit)
        if (!is.matrix(fit) || nrow(fit)!= nstate || ncol(fit)!= sum(keep))
            stop("response function must return a matrix with nstate rows and one column per response")
    }
    if (missing(penalty)) pmat <- NULL  #this will be used as a flag
    else {
        pmat <- penalty
        if (!is.matrix(pmat) || !is.numeric(pmat) || nrow(pmat) != ncol(pmat))
            stop("pmat must be a square numeric matrix")
        if (nrow(pmat) != nparm) 
            stop("nrow(pmat) != number of parameters")
        if (!all.equal(pmat, t(pmat))) stop("penalty matrix must be symmetric")
        test <- svd(pmat, nv=0)
        if (any(test$d < 0) || sign(diag(test$u)) != sign(diag(test$v)))
            stop("penalty matrix must be non-negative definite")
        }
    if (missing(constraint)) constraint <- NULL  #this will be used as a flag
    else {
        if (!is.matrix(constraint) || !is.numeric(constraint) || 
            nparm != ncol(constraint))
            stop("constraint must be a numeric matrix with one column per parameter")
        }
    if (bcount[3] ==0) p0fixed <- pfun(nstate)
    else p0fixed <- NULL
    hmm1 <- function(who, beta) {
        rows <- which(id ==uid[who])  # the subjects of interest
        eta <- X[rows,] %*% beta
        # starting probability
        if (is.null(p0fixed))
            alpha <- pfun(nstate, eta[1,b3], gradient=FALSE)
        else alpha <- p0fixed
     
        # Now the response functions for this set
        rneed <- (otype[rows]==1 | otype[rows]==3)
        rlist <- vector("list", ny)
        for (k in 1:ny) {
            j <- b2map[[k]]  #columns of beta for this response
            indx <- rneed & !is.na(yobs[rows, k])
            yy <- yobs[rows[indx], k]
            if (length(yy) >0) {
                if (length(j)==0) rlist[[k]] <-rfun[[k]](yy, nstate, gradient=FALSE)
                else rlist[[k]] <- rfun[[k]](yy, nstate, eta[indx, j, drop=FALSE], 
                                  gradient=FALSE)
            }
        }

        # Compute the collection of matrix exponentials for the subject
        # The upper routine sends back the array of results as a vector
        #  along with the number of times there were tied eigenvalues
        #  we'll send the ties back as an attribute
        # We don't need the last row for each subject.
        # The call to upper uses the Ward approx (nterm=0) rather than
        #  the Higham09.  The former seems to better match my pade routine.
        r2 <- length(rows)
        if (length(rows) > 1) {
            if (any(abs(eta[-r2,]) > .Machine$double.max.exp/2)) {
                # such a bad estimate that it may blow up the matrix exp
                if (debug >1) browser()
                return("underflow")
            }
            myexp <- .Call("upper", nstate, eta[-r2,,drop=FALSE], 
                           ytime[rows[-r2]], rindex, 1e-7, 0)
            ucount <- c(length(rows)-1, myexp$ties)
            Pmat <- array(myexp$P, dim=c(nstate, nstate, length(rows)-1))
        }
        
        # Now walk through the visits one by one
        offset <- 0   # watch out for underflow
        nc <- integer(ny)  # the number of non-censored & non-missing so far
        rmat <- matrix(0., nstate, nstate)

        for (jj in seq_along(rows)) {
            j <- rows[jj] # j is the index in the original data, jj in our subset
            if (otype[j] == 2) {
                # exact event time (death)
                rmat[rindex] <- exp(eta[jj-1, b1]) #covariate just before this point
                dtemp <- rmat[,death]
                alpha <- alpha * dtemp
                 if (debug > 2) cat("A2: j=", j, "alpha=", alpha, "\n")
            }
            else if (otype[j] == 3) {
                # entry to the study
                temp <- sum(alpha*entry)
                if (temp <= 0) 
                  return(paste("subject", uid[who],"enters in an impossible state"))
                alpha <- alpha*entry /temp
                if (debug > 2) cat("A3: j=", j, "alpha=", alpha, "\n")
            }
        
            if (otype[j]==3 || otype[j]==1) {  # an outcome was observed
                temp <- rep(1, nstate)
                for (k in 1:ny) {
                    if (!is.na(yobs[j,k])) {
                        nc[k] <- nc[k] +1
                        alpha <- alpha* rlist[[k]][,nc[k]]
                        temp <- temp * rlist[[k]][,nc[k]]
                    }
                }
                if (!all(is.finite(alpha)) || sum(alpha) <=0) {
                    if (debug > 1) browser()
                    return("underflow") 
                }
                if (debug>2) cat("B: j=", j, "alpha=", alpha, "\n")
            }

            if (jj< r2) {  # if not the last obs
                # transition matrix
                alpha <- alpha %*% Pmat[,,jj]  # transition to next time point
                if (debug > 2) cat("C: j=", j, "alpha=", alpha, "\n")

                if (!all(is.finite(alpha)) || sum(alpha) <=0) {
                    if (debug > 1) browser()
                    return("underflow")
                }
                if (mean(alpha) < exp(-20)) { # beware underflow
                    reset <- min(-20, log(mean(alpha)))
                    if (debug > 2) cat(" offset=", offset,"reset=", reset, "\n")
                    offset <- offset + reset
                    alpha <- alpha * exp(-reset)
                }
            }
        }

        loglik <- offset + log(sum(alpha))
        attr(loglik, "counts") <- ucount
        loglik
    }
    makeindex <- function(cmap, all=cmap) {
        nonzero <- (cmap > 0)
        parms <- sort(unique(all[all>0]))  # the parameter numbers for this group
        p <- ncol(cmap)    # number of linear predictors
        k <- match(cmap[nonzero], parms)  #parameter number
        nparm <- length(parms)
        rr <- row(cmap)[nonzero]  # which X to use
        cc <- col(cmap)[nonzero]  #which eta this is
        list(xindex= rr, tindex= cc + (k-1)*p, dim=c(p, nparm))
        }

    Rtrans <- vector("list", ny)  #one element per response function
    if (bcount[2]) { #if there are response parameters
        tfun <- function(dmat, x, map) {
            tmat <- matrix(0., map$dim[1], map$dim[2])
            tmat[map$tindex] <- x[map$xindex]
            dmat %*% tmat
        }
        for (i in 1:ny) {
            if (length(b2map[[i]]) >0 && any(cmap[, b2map[[i]]] > 0)) {
                formals(tfun)[[3]] <- makeindex(cmap[,b2map[[i]], drop=FALSE],
                                                cmap[,b2])
                Rtrans[[i]] <- tfun
            }
        }
    }
    if (bcount[3]) { #if there are initial probability  parameters
        pitrans <- function(dmat, x, map) {
            tmat <- matrix(0., map$dim[1], map$dim[2])
            tmat[map$tindex] <- x[map$xindex]
            dmat %*% tmat
        }
        formals(pitrans)[[3]] <- makeindex(cmap[,b3, drop=FALSE])
    }
    if (debug > 1) browser()
    cmap.b1 <- makeindex(cmap[,b1, drop=FALSE])
    Ptrans <- function(alpha, dmat, x, map=cmap.b1) {
        tmat <- matrix(0., map$dim[1], map$dim[2])
        tmat[map$tindex] <- x[map$xindex]
        #treat dmat as though it were a matrix with first dim nstate*nstate
        dim(dmat) <- c(nstate*nstate, map$dim[1])
        dmat2 <- dmat %*% tmat  #transform
        t(rowsum(dmat2 * rep(alpha, nstate*map$dim[2]), rep(1:nstate, each=nstate),
                 reorder=FALSE))
    }    
    if (!missing(death)) {  
        dtemp <- col(qmatrix)[rindex]
        deathcol  <- which(dtemp== death)  # list of linear predictors
        deathtrans <- function(R, x, map=cmap.b1) {
            rows <- which(qmatrix[,death] > 0)  #non-zero elements of column d
            dmat <- matrix(0, nstate, map$dim[1])
            for (i in 1:length(rows)) 
                dmat[rows[i], deathcol[i]] <- R[rows[i], death]

            tmat <- matrix(0., map$dim[1], map$dim[2])
            tmat[map$tindex] <- x[map$xindex]
            dmat %*% tmat
        }
    }
    psetup <- function(rmat, rindex) {
        n.eta <- length(rindex)
        out <- array(0., c(nstate, nstate, n.eta))
        temp <- matrix(0., nstate, nstate)
        rr <- row(temp)[rindex]
        for (i in 1:n.eta) {
            temp2 <- temp
            exp.eta <- rmat[rindex[i]]  # elements of rmat are exp(eta)
            temp2[rindex[i]] <-  exp.eta
            temp2[rr[i], rr[i]] <- -exp.eta
            out[,,i] <- temp2
            }
        out
    }
    hmm2 <- function(who, beta) {
        rows <- which(id ==uid[who])  # the subjects of interest
        eta <- X[rows,] %*% beta
        P.d  <- matrix(0., parmcount[1], nstate)
        R.d  <- matrix(0., nstate, parmcount[2])
        pi.d <- matrix(0., nstate, parmcount[3])

        # starting probability
        if (is.null(p0fixed)) alpha <- pfun(nstate, eta[1,b3], gradient=TRUE)
        else alpha <- p0fixed
        if (bcount[3]) {
            pi.d <- pitrans(attr(alpha, 'gradient'), X[rows[1],])
            attr(alpha, 'gradient') <- NULL  # no longer needed
        }
        
        # Execute the response functions, over the uncensored obs
        rlist <- rgrad <- vector("list", ny)
        rneed <- (otype[rows] ==1 | otype[rows]==3)  # non-censored rows
        for (k in 1:ny) {
            index <- rneed & !is.na(yobs[rows,k])
            j <- b2map[[k]]  #linear predictors for this response
            yy <- yobs[rows[index], k]
            if (length(yy) > 0) {
                temp <- rfun[[k]](yy, nstate, eta[index, j, drop=FALSE], 
                                    gradient= TRUE)
                rlist[[k]] <- temp
                rgrad[[k]] <- attr(temp, "gradient")
            }
        }

        # Walk through the observations one by one
        P.d  <- matrix(0., pcount[1], nstate)
        offset <- 0  # watch out for underflow
        ecount <- c(length(rows), 0)
        nc <- integer(ny)    #number non-censored so far
        r2 <- length(rows)
        rmat <- matrix(0., nstate, nstate)

        for (jj in seq_along(rows)) {
            j <- rows[jj]
            if (otype[j] == 2) {
                # exact event time (death)
                dtemp <- rmat[,death]  #rate at this point
                dtemp[death] <- 0      # this line should be redundant
                if (pcount[3]) pi.d <- pi.d * rep(dtemp, pcount[3])
                if (pcount[2]) R.d  <- R.d  * rep(dtemp, pcount[2])
                # Why the j-1 below?  A death density has to depend on covariates
                #  measured prior to the death, not measured at the death
                # dtemp above already has this lag, since rmat is from prior iter
                if (pcount[1]) P.d  <- P.d *  rep(dtemp, each=pcount[1]) +
                                 t(alpha * deathtrans(rmat, X[j-1,]))
                alpha <- alpha * dtemp
                if (debug > 2) {
                    cat("\n death: alpha=", format(alpha), "\n")
                    if (pcount[3]) print(pi.d)
                    if (pcount[2]) print(R.d)
                    print(P.d)
                }
                if (debug > 2) cat("A2: j=", j, "alpha=", alpha, "\n")
            }
            else if (otype[j] == 3) {
                # entry to the study
                temp <- sum(alpha*entry)
                if (temp <= 0) return(paste("subject", uid[who],
                                          "enters in an impossible state"))
                if (pcount[3]) pi.d <- (entry/temp)*( pi.d -
                                   alpha %*% (entry %*% pi.d)/temp )
                if (pcount[2]) R.d <- (entry/temp)* (R.d - 
                                   alpha %*% (entry %*% R.d)/temp)
                # remember that P.d is (nparm, nstate)
                if (pcount[1]) P.d  <- (rep(entry, each=pcount[1])/temp) * (P.d -
                                   outer(c(P.d %*% entry), alpha) /temp)
                alpha <- alpha*entry /temp
                if (debug > 2) {
                    cat("\n entry: alpha=", format(alpha), "\n")
                    if (pcount[3]) print(pi.d)
                    if (pcount[2]) print(R.d)
                    print(P.d)
                }
                if (debug > 2) cat("A3: j=", j, "alpha=", alpha, "\n")
            }

            if (otype[j]==3 || otype[j]==1) {  # an outcome was observed
                for (k in 1:ny) {
                    if (!is.na(yobs[j,k])) {
                        nc[k] <- nc[k] +1
                        temp <- rlist[[k]][,nc[k]]
                        if (pcount[3]) pi.d <- pi.d * temp 
                        if (pcount[1]) P.d  <- P.d * rep(temp, each=pcount[1])
                        if (pcount[2]) R.d  <- R.d * temp
                        if (!is.null(Rtrans[[k]])) { #if there are derivatives
                            dtemp <- Rtrans[[k]](rgrad[[k]][,nc[k],], X[j,])
                            R.d  <- R.d + alpha * dtemp
                         } 
                        alpha <- alpha * temp
                        if (debug > 3) {
                            cat("\n response: alpha=", format(alpha), "\n")
                            if (pcount[3]) print(pi.d)
                            if (pcount[2]) print(R.d)
                            print(P.d)
                        }
                    }
                }
                if (!all(is.finite(alpha)) || sum(alpha) <=0) {
                    if (debug>1) browser()
                    return("underflow")
                }
                if (debug > 2) cat("B: j=", j, "alpha=", alpha, "\n")
            }
            
            if (jj < r2) { # not the last row
                # state matrix transformation P
                rmat[rindex] <- exp(eta[jj,b1])
                if (!all(is.finite(rmat))) {
                    # a horrible beta can overflow
                    if (debug > 1) browser()
                    return("overflow")
                }
                
                diag(rmat) <- diag(rmat) -rowSums(rmat)
                tder <- psetup(rmat, rindex)
                if (any(diff(sort(diag(rmat))) < 1e-6)) {
                    ptemp <- pade(rmat *ytime[j], tder*ytime[j])
                    ecount[2] <- ecount[2] +1
                }
                else ptemp <- derivative(rmat, ytime[j], tder)
                if (pcount[3]) pi.d <- t(ptemp$P) %*% pi.d 
                if (pcount[2]) R.d <-  t(ptemp$P) %*% R.d 
                if (pcount[1]) 
                    P.d <-  P.d %*% ptemp$P + Ptrans(alpha, ptemp$dmat, X[j,])
                alpha <- drop(alpha %*% ptemp$P)   # ditch the dimensions
                if (debug > 2) {
                    cat("\n j=", j, "jj=", jj, "alpha=", format(alpha), "\n")
                    if (pcount[3]) print(pi.d)
                    if (pcount[2]) print(R.d)
                }
                if (debug > 1) cat("C: j=", j, "alpha=", alpha, "\n")

                if (!all(is.finite(alpha)) || sum(alpha) <=0) {
                    if (debug > 1) browser()
                    return("underflow")
                }
                if (mean(alpha) < exp(-20)) {
                    offset <- offset -20
                    alpha <- alpha * exp(20)
                    pi.d <- pi.d * exp(20)
                    R.d  <- R.d  * exp(20)
                    P.d  <- P.d  * exp(20)
                }
            }
        }
        if (debug > 1) cat("alpha=", alpha, "\n")
        list(alpha=alpha, deriv= rbind(P.d, t(R.d), t(pi.d)), ecount=ecount,
             offset = offset)
    }
    uid <- unique(id)
    nid <- length(uid)
    rindex <- which(qmatrix >0)

    #Give this next variable a long name that won't be found in calling
    #  routines.  It is updated farther down the calling chain.
    hmm_count_of_calls <- c(0, 0)  #total calls to expm, number with tied eigens

    if (mc.cores > 1 && makefork)
        hmm_cluster <- makeForkCluster(mc.cores) #start up parallel

    hmmloglik <- function(param, ...) {
        beta[cmap>0] <- param[cmap]
        if (mc.cores > 1) {
            if (makefork)
                mcfit <- parLapply(hmm_cluster, 1:nid, hmm1, beta=beta)
            else mcfit <- mclapply(1:nid, hmm1, beta=beta, 
                                   mc.set.seed=FALSE, mc.cores=mc.cores)
        }
        else mcfit <- lapply(1:nid, hmm1, beta=beta)
        
        if (any(sapply(mcfit, is.character))) {
            # failure
            words <- sapply(mcfit, function(x) ifelse(is.character(x), x, ""))
            if (any(words == "underflow")) {
                if (debug > 1) browser()
                # assume a bad guess from a maximizer, return a bad hit
                return(-2 * abs(initial.loglik))
            }
            words <- words[words!=""]
            stop(words[1])
        }

        total <- sum(unlist(mcfit))
        if (!is.null(pmat)) total <- total - c(param %*% pmat %*% param)/2
        count <- sapply(mcfit, function(x) attr(x, "counts"))
        # this next line reaches back and changes the variable in parent of parent
        #  (the parent called the maximizer, which calls this)
        hmm_count_of_calls <<- hmm_count_of_calls + rowSums(count)
        total
    }

    # hand back everything (debug)
    hmmdb <- function(param, ...) {
        beta[cmap>0] <- param[cmap]
        if (mc.cores > 1) {
            if (makefork)
                mcfit <- parLapply(hmm_cluster, 1:nid, hmm2, beta=beta)
            else mcfit <- mclapply(1:nid, hmm2, beta=beta, 
                                   mc.set.seed=FALSE, mc.cores=mc.cores)
        }
        else mcfit <- lapply(1:nid, hmm2, beta=beta)

        if (any(sapply(mcfit, is.character))) {
            # failure
            words <- sapply(mcfit, function(x) ifelse(is.character(x), x, ""))
            if (any(words == "underflow")) {
                # assume a bad guess from a maximizer, return a bad hit
                return(-2 * abs(initial.loglik))
            }
            words <- words[words!=""]
            stop(words[1])
        }

        ecount <- sapply(mcfit, function(x) x$ecount)
        hmm_count_of_calls <<- hmm_count_of_calls + rowSums(ecount)
        dd <- dim(mcfit[[1]]$deriv)
        alpha <- sapply(mcfit, function(x) x$alpha)
        list(alpha = alpha,
             offset = sapply(mcfit, function(x) x$offset),
             loglik = sum(log(colSums(alpha))),
             deriv = array(unlist(lapply(mcfit, function(x) x$deriv)),
                           dim=c(dd, length(mcfit))),
             ecount=ecount, offset= sapply(mcfit, function(x) x$offset))
    }
    hmmboth <- function(param, ...) {
        beta[cmap>0] <- param[cmap]
        if (mc.cores > 1) {
            if (makefork)
                mcfit <- parLapply(hmm_cluster, 1:nid, hmm2, beta=beta)
            else mcfit <- mclapply(1:nid, hmm2, beta=beta, 
                                   mc.set.seed=FALSE, mc.cores=mc.cores)
        }
        else mcfit <- lapply(1:nid, hmm2, beta=beta)
        
        if (any(sapply(mcfit, is.character))) {
            # failure
            words <- sapply(mcfit, function(x) ifelse(is.character(x), x, ""))
            if (any(words == "underflow")) {
                # assume a bad guess from a maximizer, return a bad hit
                return(list(loglik = -2 * abs(initial.loglik), deriv=NULL, S=NULL))
            }
            words <- words[words!=""]
            stop(words[1])
        }

        ecount <- rowSums(sapply(mcfit, function(x) x$ecount))
        hmm_count_of_calls <<- hmm_count_of_calls + ecount
        
        alpha <- sapply(mcfit, function(x) sum(x$alpha))
        offset <- sapply(mcfit, function(x) x$offset)
        d.alpha <- sapply(mcfit, function(x) rowSums(x$deriv))
        u <- d.alpha * rep(1/alpha, each=nparm)
        # this will be a matrix with nparm rows, one col per subject
        deriv <- rowSums(u)
        S = u %*% t(u)
        u2 <- u - rowMeans(u)
        S2 <- u2 %*% t(u2)
        
        if (is.null(pmat)) loglik <- sum(log(alpha) + offset)
        else {
            loglik <- sum(log(alpha)+offset) - c(param %*% pmat %*% param)/2
            deriv  <- deriv - c(param %*% pmat)
            S <- S + pmat
            S2 <- S2 + pmat
            }
        list(loglik = loglik, deriv= deriv, S=S, S2=S2)
        }

    hmmgrad <- function(param, ...) {
        beta[cmap>0] <- param[cmap]
        if (mc.cores > 1) {
            if (makefork)
                mcfit <- parLapply(hmm_cluster, 1:nid, hmm2, beta=beta)
            else mcfit <- mclapply(1:nid, hmm2, beta=beta, 
                                   mc.set.seed=FALSE, mc.cores=mc.cores)
        }
        else mcfit <- lapply(1:nid, hmm2, beta=beta)
        
        if (any(sapply(mcfit, is.character))) {
            # failure
            words <- sapply(mcfit, function(x) ifelse(is.character(x), x, ""))
            if (any(words == "underflow")) {
                # assume a bad guess from a maximizer, return a bad hit
                return(-2 * abs(initial.loglik))
            }
            words <- words[words!=""]
            stop(words[1])
        }

        ecount <- rowSums(sapply(mcfit, function(x) x$ecount))
        hmm_count_of_calls <<- hmm_count_of_calls + ecount
        
        alpha <- sapply(mcfit, function(x) sum(x$alpha))
        d.alpha <- sapply(mcfit, function(x) rowSums(x$deriv))
        # this will be a matrix with npar rows, one col per subject
        u <- d.alpha * rep(1/alpha, each=nparm)

        if (is.null(pmat)) rowSums(u)
        else rowSums(u) - c(param %*% pmat)
        }
    # Get an initial loglik
    # If the maximizer fails (NA or infinite) it returns the initial;
    #  set it to a dummy value to detect that
    initial.loglik <- numeric(0)
    initial.loglik <- hmmloglik(param)
    if (length(initial.loglik)==0)
        stop("unable to evalutate the likelihood at the intial parameters")

    mpar$par <- param   # the initial parameters
    if (is.null(mpar$fn)) mpar$fn  <- hmmloglik
    else mpar$fn <- get(mpar$fn)
    if (!is.null(mpar$gr) && is.character(mpar$gr)) mpar$gr <- get(mpar$gr)
    if (!is.null(constraint) && is.null(mpar["constraint"]))
        mpar$constraint <- constraint

    #if (is.null(mpar$cmap)) mpar$cmap <- cmap

    time1 <- proc.time()
    fit <- do.call(mfun, mpar)
    if (mc.cores > 1 & makefork) stopCluster(hmm_cluster)
    time2 <- proc.time() 
    bcol <- paste0(row(qmatrix)[qmatrix!=0], ":",
                   col(qmatrix)[qmatrix!=0])
    if (bcount[2]>0) 
        bcol <- c(bcol, paste0("R", paste(rcoef$response, rcoef$lp, sep='.')))
    if (bcount[3]>0) bcol <- c(bcol, paste0("p", 1:bcount[3]))
    dimnames(beta) <- list(dimnames(X)[[2]], bcol)
    dimnames(cmap) <- dimnames(beta)

    nfit <- names(fit)
    indx <- pmatch(c("coef", "par", "log", "value"), nfit, nomatch=0)

    # find the fitted coefs in the output, and the loglik
    fcoef <- if (indx[1] >0) fit[[indx[1]]] else
                 if (indx[2]>0) fit[[indx[2]]] else {
                     zz <- seq_along(nfit)[-indx]
                     fit[[zz[1]]]
                 }

    if (is.matrix(fcoef)) {
        if (ncol(fcoef) == nparm) {
            beta[cmap>0] <- fcoef[nrow(fcoef), cmap]
        } else stop("wrong number of columns in coefficient matrix")
    }
    else beta[cmap>0] <- fcoef[cmap]

    flog <-  if (indx[3] >0) fit[[indx[3]]] else
                 if (indx[4]>0) fit[[indx[4]]] else NULL

    time3 <- proc.time()
    compute.time <- rbind(setup= time1-time0,
                          compute= time2- time1,
                          finish = time3 - time2)

    final <- list(coefficients= fcoef, 
                  loglik = flog,
                  beta=beta,
                  fit=fit, 
                  loglik0 = initial.loglik,
                  evals= hmm_count_of_calls, 
                  time = compute.time,
                  scale =  scale,
                  bcount = bcount,
                  cmap=cmap, rmap=rindex,
                  qmatrix = qmatrix,   # the structure and state names
                  rfun= rfun, pfun=pfun,
                  nstate = nstate,
                  n = nrow(mf),
                  na.action = na.action,
                  call=Call)
    class(final) <- "hmm"
    if (debug) browser()
    final
}

hbind <- function(...) {
    temp <- list(...)
    tlevel <- lapply(temp, function(x) {
                     if (!is.numeric(x)) levels(as.factor(x))
                     else NULL })
    new <- do.call("cbind", lapply(temp, function(x) {
          if (is.numeric(x)) x  else as.numeric(as.factor(x))}))
    attr(new, "levels") <- tlevel
    class(new) <- "hbind"
    new
}
