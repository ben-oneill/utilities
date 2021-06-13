#' Nonlinear minimisation/maximisation allowing probability vectors as inputs
#'
#' \code{nlm.prob} minimises/maximises a function allowing probability vectors as inputs
#'
#' This is a variation of the \code{stats::nlm} function for nonlinear minimisation.  The present function is designed to
#' minimise an objective function with one or more arguments that are probability vectors.  (The objective function may also
#' have other arguments that are not probability vectors.)  The function uses the same inputs as the \code{stats::nlm} function,
#' except that the user can use the input \code{prob.vectors} to specify which inputs are constrained to be probability vectors.
#' This input is a list where each element in the list specifies a set of indices for the argument of the objective function; the
#' specified set of indices is constrained to be a probability vector (i.e., each corresponding argument is non-negative and the
#' set of these arguments must sum to one).  The input \code{prob.vectors} may list one or more probability vectors, but they must
#' use disjoint elements of the argument (i.e., a variable in the argument cannot appear in more than one probability vector).
#'
#' Optimisation is performed by first converting the objective function into unconstrained form using the softmax transformation
#' and its inverse to convert from unconstrained space to probability space and back.  Optimisation is done on the unconstrained
#' objective function and the results are converted back to probability space to solve the constrained optimisation problem.  For
#' purposes of conversion, this function allows specification of a tuning parameter \code{lambda} for the softmax and inverse-
#' softmax transformations.  (This input can either be a single tuning value used for all conversions, or a vector of values for
#' the respective probability vectors; if the latter, there must be one value for each element of the \code{prob.vector} input.)
#'
#' Most of the input descriptions below are adapted from the corresponding descriptions in \code{stat::nlm}, since our function is
#' a wrapper to that function.  The additional inputs for this function are \code{prob.vectors}, \code{lambda} and \code{eta0max}.
#' The function also adds an option \code{maximise} to conduct maximisation instead of minimisation.
#'
#' @param f The objective function to be minimised; output should be a single numeric value.
#' @param p Starting argument values for the minimisation.
#' @param prob.vectors A list specifying which sets of elements are constrained to be a probability vector (each element in the list
#' should be a vector specifying indices in the argument vector; elements cannot overlap into multiple probability vectors).
#' @param lambda The tuning parameter used in the softmax transformation for the optimisation (a single positive numeric value).
#' @param eta0max The maximum absolute value for the elements of eta0 (the starting value in the unconstrained optimisation problem).
#' @param maximise,maximize Logical value; if \code{TRUE} the function maximises the objective function instead of mimimising.
#' @param hessian Logical; if \code{TRUE} then the output of the function includes the Hessian of \code{f} at the minimising point.
#' @param typsize An estimate of the size of each parameter at the minimum.
#' @param fscale An estimate of the size of \code{f} at the minimum.
#' @param print.level This argument determines the level of printing which is done during the minimisation process. The default value
#' of \code{0} means that no printing occurs, a value of \code{1} means that initial and final details are printed and a value of \code{2}
#' means that full tracing information is printed.
#' @param ndigit The number of significant digits in the function \code{f}.
#' @param gradtol A positive scalar giving the tolerance at which the scaled gradient is considered close enough to zero to terminate
#' the algorithm. The scaled gradient is a measure of the relative change in \code{f} in each direction \code{p[i]} divided by the
#' relative change in \code{p[i]}.
#' @param stepmax A positive scalar which gives the maximum allowable scaled step length.  \code{stepmax} is used to prevent steps which
#' would cause the optimisation function to overflow, to prevent the algorithm from leaving the area of interest in parameter space, or
#' to detect divergence in the algorithm.  \code{stepmax} would be chosen small enough to prevent the first two of these occurrences,
#' but should be larger than any anticipated reasonable step.
#' @param steptol A positive scalar providing the minimum allowable relative step length.
#' @param iterlim A positive integer specifying the maximum number of iterations to be performed before the routine is terminated.
#' @param check.analyticals Logical; if \code{TRUE} then the analytic gradients and Hessians (if supplied) are checked against numerical
#' derivatives at the initial parameter values.  This can help detect incorrectly formulated gradients or Hessians.
#' @param ... Additional arguments to be passed to \code{f} via \code{nlm}
#' @return A list showing the computed minimising point and minimum of \code{f} and other related information.
#'
#' @examples
#' x <- rbinom(100, 1, .2)
#' nlm.prob(function(p) sum(dbinom(x,1,p[2],log=TRUE)), c(.5, .5), maximise = TRUE)

nlm.prob <- function(f, p, prob.vectors = list(1:length(p)), ..., lambda = 1,
                     eta0max = 1e10, maximise = FALSE, maximize = maximise,
                     hessian = FALSE, typsize = rep(1, length(p)),
                     fscale = 1, print.level = 0, ndigit = 12, gradtol = 1e-06,
                     stepmax = max(1000*sqrt(sum((p/typsize)^2)), 1000),
                     steptol = 1e-06, iterlim = 100, check.analyticals = TRUE) {

  #Set printing level
  print.level <- as.integer(print.level)
  if (!(print.level %in% c(0,1,2))) { stop("'print.level' must be in {0,1,2}") }
  MSG <- ifelse(check.analyticals,  c(9, 1, 17)[1+ print.level],
                c(15, 7, 23)[1+ print.level])

  #Check inputs prob.vectors and p (this must be a list of vectors of indices over p)
  #Check that there are no elements in overlapping vectors
  #Check that all probability vectors are non-negative and sum to one
  m <- length(p)
  s <- length(prob.vectors)
  for (i in 1:s) {
    INDICES <- prob.vectors[[i]]
    if (!is.numeric(INDICES))                    { stop('Error: Each element of prob.vectors should be a vector of indices for the input p')  }
    if (!all(as.integer(INDICES) == INDICES))    { stop('Error: Each element of prob.vectors should be a vector of indices for the input p')  }
    if (min(INDICES) < 1)                        { stop('Error: Each element of prob.vectors should be a vector of indices for the input p')  }
    if (max(INDICES) > m)                        { stop('Error: Each element of prob.vectors should be a vector of indices for the input p')  }
    INDICES <- sort(unique(INDICES))
    WARN <- FALSE
    if (min(p[INDICES]) <  0) {
      warning(paste0('Error: prob.vector ', i, ' has a negative element'))
      WARN <- TRUE }
    if (sum(p[INDICES]) != 1) {
      warning(paste0('Error: prob.vector ', i, ' does not sum to one'))
      WARN <- TRUE }
    if (WARN) {
      PVEC <- pmax(1e-10, p[INDICES])
      PVEC <- PVEC/sum(PVEC)
      p[INDICES] <- PVEC
      warning(paste0('We have adjusted the starting values of prob.vector ', i, ' to use a valid probability vector\n')) } }
  if (length(unique(unlist(prob.vectors))) < length(unlist(prob.vectors))) {
                                                   stop('Error: Lists of indices in prob.vectors must not overlap') }

  #Check input lambda
  if (!is.numeric(lambda))                       { stop('Error: Input lambda should be a numeric value') }
  if (min(lambda) <= 0)                          { stop('Error: Input lambda should be a positive value') }
  ss <- length(lambda)
  if ((ss != 1)&(ss != s))                       { stop('Error: Input lambda should either be a single value, or it should have the same length as prob.vectors') }
  if (ss == 1) { lambda <- rep(lambda, s) }

  #Check input maximise
  if (!missing(maximise) && !missing(maximize)) {
    if (maximise != maximize) {
      warning("Specify 'maximise' or 'maximize' but not both") } else {
         stop("Error: specify 'maximise' or 'maximize' but not both") } }
  MAX <- maximize
  if (!isTRUE(MAX) && ! isFALSE(MAX))            { stop('Error: Input maximise/maximize should be a single logical value') }

  #Convert input prob.vector to list of argument-indices
  #The object ARGS.INDEX is a list of argument indices
  #The object ARGS.LENGTH is a vector of their lengths
  ARGS.LENGTH <- rep(0, s+1)
  ARGS.INDEX  <- vector(mode = "list", length = s+1)
  names(ARGS.INDEX)[1:s] <- sprintf('prob.vector.%s', 1:s)
  names(ARGS.INDEX)[s+1] <- 'other.args'
  OTHER <- rep(TRUE, m)
  for (i in 1:s) {
    IND <- prob.vectors[[i]]
    OTHER[IND] <- FALSE
    ARGS.INDEX[[i]] <- IND
    ARGS.LENGTH[i]  <- length(IND) }
  ARGS.INDEX[[s+1]] <- which(OTHER)
  ARGS.LENGTH[s+1]  <- sum(OTHER)

  #Create mapping from eta to p
  #Eta is the unconstrained input vector
  #p is the constrained input vector (with probability vectors)
  eta_to_p <- function(eta) {

    #Check input
    if (!is.numeric(eta))      { stop('Error: Input eta must be numeric') }
    if (length(eta) != m-s)    { stop('Error: Input eta must have length m-s') }

    #Generate piecewise mapping and derivatives
    ARGS <- vector(mode = "list", length = s+1)
    DD1  <- vector(mode = "list", length = s+1)
    DD2  <- vector(mode = "list", length = s+1)
    t <- 0
    for (i in 1:s) {
      r <- ARGS.LENGTH[i]
      if (r == 1) {
        ARGS[[i]] <- 1 }
      if (r > 1) {
        SOFT <- softmax(eta[(t+1):(t+r-1)], lambda = lambda[i], gradient = TRUE, hessian = TRUE)
        ARGS[[i]] <- c(SOFT)/sum(c(SOFT))
        DD1[[i]]  <- attributes(SOFT)$gradient
        DD2[[i]]  <- attributes(SOFT)$hessian }
      t <- t+r-1 }
    r <- ARGS.LENGTH[s+1]
    ARGS[[s+1]] <- eta[(t+1):(t+r)]
    DD1[[s+1]]  <- diag(r)
    DD2[[s+1]]  <- array(0, dim = c(r, r, r))

    #Generate output vector (labelled PPP)
    PPP <- rep(NA, m)
    D1  <- array(0, dim = c(m, m-s))
    D2  <- array(0, dim = c(m, m-s, m-s))
    t <- 0
    for (i in 1:s) {
      r <- ARGS.LENGTH[i]
      IND <- ARGS.INDEX[[i]]
      PPP[IND] <- ARGS[[i]]
      if (r > 1) {
        D1[IND, (t+1):(t+r-1)] <- DD1[[i]]
        D2[IND, (t+1):(t+r-1), (t+1):(t+r-1)] <- DD2[[i]] }
      t <- t+r-1 }
    r <- ARGS.LENGTH[s+1]
    if (r > 0) {
      IND <- ARGS.INDEX[[s+1]]
      PPP[IND] <- ARGS[[s+1]]
      D1[IND, (t+1):(t+r)] <- DD1[[s+1]]
      D2[IND, (t+1):(t+r), (t+1):(t+r)] <- DD2[[s+1]] }
    attr(PPP, 'gradient') <- D1
    attr(PPP, 'hessian')  <- D2

    #Give output
    PPP }

  #Create mapping from p to eta
  #Eta is the unconstrained input vector
  #p is the constrained input vector (with probability vectors)
  p_to_eta <- function(p) {

    #Generate piecewise mapping and derivatives
    ARGS <- vector(mode = "list", length = s+1)
    DD1  <- vector(mode = "list", length = s+1)
    DD2  <- vector(mode = "list", length = s+1)
    for (i in 1:s) {
      IND <- ARGS.INDEX[[i]]
      SOFTINV <- softmaxinv(p[IND], lambda = lambda[i], gradient = TRUE, hessian = TRUE)
      ARGS[[i]] <- c(SOFTINV)
      DD1[[i]]  <- attributes(SOFTINV)$gradient
      DD2[[i]]  <- attributes(SOFTINV)$hessian }
    r <- ARGS.LENGTH[s+1]
    if (r > 0) {
      IND <- ARGS.INDEX[[s+1]]
      ARGS[[s+1]] <- p[IND]
      DD1[[s+1]]  <- diag(r)
      DD2[[s+1]]  <- array(0, dim = c(r, r, r)) }

    #Generate output vector (labelled EEE)
    EEE <- unlist(ARGS)
    D1  <- array(0, dim = c(m-s, m))
    D2  <- array(0, dim = c(m-s, m, m))
    t <- 0
    for (i in 1:s) {
      r <- ARGS.LENGTH[i]
      if (r > 1) {
        IND <- ARGS.INDEX[[i]]
        D1[(t+1):(t+r-1), IND]      <- DD1[[i]]
        D2[(t+1):(t+r-1), IND, IND] <- DD2[[i]] }
      t <- t+r-1 }
    r <- ARGS.LENGTH[s+1]
    if (r > 0) {
      IND <- ARGS.INDEX[[s+1]]
      D1[(t+1):(t+r), IND]      <- DD1[[s+1]]
      D2[(t+1):(t+r), IND, IND] <- DD2[[s+1]] }
    attr(EEE, 'gradient') <- D1
    attr(EEE, 'hessian')  <- D2

    #Give output
    EEE }

  #Create unconstrained objective function
  SGN <- ifelse(MAX, -1, 1)
  OBJ <- function(eta, ...) {

    #Compute output (labelled GG)
    PP <- eta_to_p(eta)
    GG <- SGN*f(c(PP), ...)

    #Compute gradient (if available)
    GRAD.f <- SGN*attributes(GG)$gradient
    if (is.matrix(GRAD.f)) {
      GRAD.p <- attributes(PP)$gradient
      D1 <- GRAD.f %*% GRAD.p
      attr(GG, 'gradient') <- D1 }

    #Compute Hessian (if available)
    HESS.f <- SGN*attributes(GG)$hessian
    if ((is.matrix(GRAD.f))&&(is.matrix(HESS.f))) {
      HESS.p <- attributes(PP)$hessian
      T1 <- (t(GRAD.p) %*% (HESS.f %*% GRAD.p))
      T2 <- matrix(0, m-1, m-1)
      for (i in 1:(m-1)) {
        for (j in 1:(m-1)) {
          T2[i,j] <- sum(GRAD.f*HESS.p[ ,i,j]) } }
      attr(GG, 'hessian') <- T1 + T2 }

    #Give output
    GG }

  #Compute nonlinear minimisation (via unconstrained objective function)
  eta0 <- c(p_to_eta(p))
  eta0 <- pmin(pmax(eta0, -eta0max), eta0max)
  NLM <- stats::nlm(f = OBJ, p = eta0, hessian = hessian, typsize = p_to_eta(typsize),
                    fscale = fscale, print.level = print.level, ndigit = ndigit,
                    gradtol = gradtol, stepmax = stepmax, steptol = steptol, iterlim = iterlim,
                    check.analyticals = TRUE, ...)

  #Convert back to probability space
  ESTIMATE <- c(eta_to_p(NLM$estimate))
  OPT      <- SGN*NLM$minimum
  GRAD.OPT <- attributes(f(ESTIMATE))$gradient
  HESS.OPT <- attributes(f(ESTIMATE))$hessian

  #Generate output list
  if (MAX) {
    if (hessian) {
      OUT <- list(maximum = OPT, estimate = ESTIMATE,
                  gradient = GRAD.OPT, hessian = HESS.OPT,
                  code = NLM$code, iterations = NLM$iterations)
      } else {
      OUT <- list(maximum = OPT, estimate = ESTIMATE, gradient = GRAD.OPT,
                  code = NLM$code, iterations = NLM$iterations) }
  } else {
    if (hessian) {
      OUT <- list(minimum = OPT, estimate = ESTIMATE,
                  gradient = GRAD.OPT, hessian = HESS.OPT,
                  code = NLM$code, iterations = NLM$iterations)
      } else {
      OUT <- list(minimum = OPT, estimate = ESTIMATE, gradient = GRAD.OPT,
                  code = NLM$code, iterations = NLM$iterations) } }

  #Give output
  OUT }
