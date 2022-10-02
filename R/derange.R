#' Generate derangement of an input vector
#'
#' \code{derange} returns one or more pseudo-random derangements of the input vector
#' \code{derange.int} returns one or more pseudo-random derangements of the values 1,...,n
#'
#' A "derangement" is a permutation with no element mapped to itself (i.e., no fixed points) and a "generalised
#' derangement" is a permutation with a specified number of fixed points.  The present funciton generates pseudo-
#' random derangements or generalised derangements of the input vector.  The input vector for the function should
#' be either numeric, integer or character type.  The function will generate the desired number of derangements
#' taken over the input vector.  For a single derangement the output is a vector and for multiple derangements
#' the output is a matrix with each row representing one derangement.
#'
#' @usage \code{derange(set, size = 1, fixed.points = 0)}
#' @usage \code{derange.int(n, size = 1, fixed.points = 0)}
#' @param n Number of elements to derange (must be at least two)
#' @param set A vector of elements to derange (must have at least two elements to derange)
#' @param size A non-negative integer specifying the number of derangements to generate
#' @param fixed.points Number of fixed points for a generalised derangement
#' @return A vector/matrix of derangements

derange <- function(set, size = 1, fixed.points = 0) {

  #Check input set
  if (!is.vector(set))                          stop('Error: Input set should be a vector')
  n <- length(set)
  if (length(set) < 2)                          stop('Error: Input set must have at least two elements')
  TYPE <- class(set)
  if (!(('numeric' %in% TYPE)|('integer' %in% TYPE)|('character' %in% TYPE))) {
                                                stop('Error: Input set should be a numeric, integer or character vector') }

  #Check input size
  if (!is.vector(size))                         stop('Error: Input size should be a numeric value')
  if (!is.numeric(size))                        stop('Error: Input size should be numeric')
  if (length(size) != 1)                        stop('Error: Input size should be a single numeric value')
  if (size != as.integer(size))                 stop('Error: Input size should be a non-negative integer')
  if (size == 0) { return(NULL) }
  if (size < 0)                                 stop('Error: Input size should be a non-negative integer')

  #Check input fixed.points
  if (!is.vector(fixed.points))                 stop('Error: Input fixed.points should be a numeric value/vector')
  if (!(length(fixed.points) %in% c(1, size)))  stop('Error: Input fixed.points does not match size')
  if (length(fixed.points) == 1) {
    FIXED <- rep(fixed.points, size)
    if (!(FIXED[1] %in% 0:n)|(FIXED[1] == n-1))  stop(paste0('Error: Input fixed.points is not a possible value'))
  } else {
    FIXED <- fixed.points
    for (k in 1:size) {
      if (!(FIXED[k] %in% 0:n)|(FIXED[k] == n-1))  stop(paste0('Error: Element ', k,' of input fixed.points is not a possible value')) } }

  #Determine type of input set and generate output
  TYPE <- class(set)
  if (('numeric' %in% TYPE)|('integer' %in% TYPE)) {
    OUT <- matrix(0,  nrow = size, ncol = n)
  } else {
    OUT <- matrix('', nrow = size, ncol = n) }
  rownames(OUT) <- sprintf('D[%s]', 1:size)
  colnames(OUT) <- set[1:n]

  #Generate generalised derangements
  for (k in 1:size) {

    #Add fixed points
    FF <- FIXED[k]
    if (FF > 0) {
      FIXED.VALS <- sort(sample.int(n, size = FF, replace = FALSE))
      OUT[k, FIXED.VALS] <- set[FIXED.VALS] }

    #Generate derangement over remaining elements
    if (FF < n) {
      nn <- n-FF
      DERANGE  <- sample.int(nn, size = nn, replace = FALSE)
      EXCESS.FIXED <- sum(DERANGE == 1:nn)
      while (EXCESS.FIXED > 0) {
        i <- which(DERANGE == 1:nn)[1]
        j <- (1:nn)[-i][sample.int(nn-1, size = 1)]
        SWAP    <- DERANGE[j]
        DERANGE[j] <- DERANGE[i]
        DERANGE[i] <- SWAP
        EXCESS.FIXED <- sum(DERANGE == 1:nn) }
      if (FF == 0) { OUT[k, ] <- set[DERANGE]
            } else { OUT[k, -FIXED.VALS] <- set[-FIXED.VALS][DERANGE] } } }

  #Return output
  OUT[1:k, ] }


derange.int <- function(n, size = 1, fixed.points = 0) {

  #Check input n
  if (!is.vector(n))                            stop('Error: Input n should be a numeric value')
  if (!is.numeric(n))                           stop('Error: Input n should be a numeric value')
  if (length(n) != 1)                           stop('Error: Input n should be a single numeric value')
  if (n < 2)                                    stop('Error: Input n must be at least two')

  #Check input size
  if (!is.vector(size))                         stop('Error: Input size should be a numeric value')
  if (!is.numeric(size))                        stop('Error: Input size should be numeric')
  if (length(size) != 1)                        stop('Error: Input size should be a single numeric value')
  if (size != as.integer(size))                 stop('Error: Input size should be a non-negative integer')
  if (size == 0) { return(NULL) }
  if (size < 0)                                 stop('Error: Input size should be a non-negative integer')

  #Check input fixed.points
  if (!is.vector(fixed.points))                 stop('Error: Input fixed.points should be a numeric value/vector')
  if (!(length(fixed.points) %in% c(1, size)))  stop('Error: Input fixed.points does not match size')
  if (length(fixed.points) == 1) {
    FIXED <- rep(fixed.points, size)
    if (!(FIXED[1] %in% 0:n)|(FIXED[1] == n-1))  stop(paste0('Error: Input fixed.points is not a possible value'))
  } else {
    FIXED <- fixed.points
    for (k in 1:size) {
      if (!(FIXED[k] %in% 0:n)|(FIXED[k] == n-1))  stop(paste0('Error: Element ', k,' of input fixed.points is not a possible value')) } }

  #Determine type of input set and generate output
  OUT <- matrix(0,  nrow = size, ncol = n)
  rownames(OUT) <- sprintf('D[%s]', 1:size)
  colnames(OUT) <- 1:n

  #Generate generalised derangements
  for (k in 1:size) {

    #Add fixed points
    FF <- FIXED[k]
    if (FF > 0) {
      FIXED.VALS <- sort(sample.int(n, size = FF, replace = FALSE))
      OUT[k, FIXED.VALS] <- (1:n)[FIXED.VALS] }

    #Generate derangement over remaining elements
    if (FF < n) {
      nn <- n-FF
      DERANGE  <- sample.int(nn, size = nn, replace = FALSE)
      EXCESS.FIXED <- sum(DERANGE == 1:nn)
      while (EXCESS.FIXED > 0) {
        i <- which(DERANGE == 1:nn)[1]
        j <- (1:nn)[-i][sample.int(nn-1, size = 1)]
        SWAP    <- DERANGE[j]
        DERANGE[j] <- DERANGE[i]
        DERANGE[i] <- SWAP
        EXCESS.FIXED <- sum(DERANGE == 1:nn) }
      if (FF == 0) { OUT[k, ] <- (1:n)[DERANGE]
      } else { OUT[k, -FIXED.VALS] <- (1:n)[-FIXED.VALS][DERANGE] } } }

  #Return output
  OUT[1:k, ] }
