\name{NEWS}
\title{NEWS file for the hmm package}
\section{Changes in version 1.1-6}{
  \itemize{
    \item Add the hmmdebug routine, to return everything.
    \item Add the constraint argument
    \item Improve the tests
    \item Fix a bug in the deriviatives: if one of the response
    functions had a fixed parameter indexing would be off
    \item Add another check: some coefficients could overflow the
    matrix exponential, giving probabilities > 1
  }
}

\section{Changes in version 1.1-5}{
  \itemize{
    \item Add the hmmemat and hmmesetup routines, to make error matrices
    easier to set up.
}}
\section{Changes in version 1.1-4}{
  \itemize{
    \item Add continuous response functions hmmncut and hmmlcut.
    \item Add backtracking to the hmmscore function, when it gets a bad
    guess.  This speeds recovery.  Change the default shrinkage to 4.
}}
\section{Changes in version 1.1-3}{
  \itemize{
    \item Minor tuning to the update steps in hmmscore and to the
    detection of a bad parameter guess (is.finite instead of is.na).
    \item Add further checks that a penalty matrix is legal.
}}
\section{Changes in version 1.1-1}{
  \itemize{
    \item Better error message if the very first guess overflows exp().
    \item Add the Fisher scoring method.
    \item Add a simple b'Pb penalty
}}

\section{Changes in version 1.1-0}{
  \itemize{
    \item Derivatives now work!
    \item Change the scale argument to be c(TRUE, FALSE) by default,
    which means to center, not scale.  Do not return beta to the
    original covariate measure, but do return the centering and scaling
    factors.
}}
  
\section{Changes in version 1.0-5}{
  \itemize{
    \item Response functions now depend on the current X values
    \item Further work on derivatives
}}

\section{Changes in version 1.0-3}{
  \itemize{
    \item Add processing for the \code{icoef} argument.
    \item Add the initial loglik to the hmm return structure.
    \item Clean up manual pages and the NAMESPACE file so that the
    package passes R CMD check.
}}