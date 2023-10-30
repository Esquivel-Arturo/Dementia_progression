print.hmm <- function(x, digits=max(options()$digits - 4, 3), ...) {
     if (!is.null(cl<- x$call)) {
	cat("Call:\n")
	dput(cl)
	cat("\n")
	}
   
     printCoefmat(x$beta, has.Pvalue=FALSE)
     cat("\n")
     cat("Log-likelihood: initial=", x$loglik0, 
         " final=", format(round(x$loglik[length(x$loglik)], 2)), "\n")

     cat("n=", x$n)
     if (length(x$na.action)) 
         cat("\  (", naprint(x$na.action), ")\n", sep="")
     else cat("\n")
     invisible(x)
}
