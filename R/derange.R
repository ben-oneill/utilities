#' Generate derangement of an input vector
#'
#' \code{derange} returns one or more pseudo-random derangements of the input vector
#'
#' A "derangement" is a permutation with no element mapped to itself (i.e., no fixed points).  This function
#' generates pseudo-random derangement of the input vector.  The input vector for the function should be
#' either numeric, integer or character type.  The function will generate the desired number of derangements
#' taken over the input vector.  For a single derangement the output is a vector and for multiple derangements
#' the output is a matrix with each row representing one derangement.
#'
#' @usage \code{derange(set, size = 1)}
#' @param set A vector of elements to derange (ust have at least two elements to derange)
#' @param size A non-negative integer specifying the number of derangements to generate
#' @return A vector/matrix of derangements of the input set

derange <- function(set, size = 1) {

  #Check input set
  if (!is.vector(set))                    stop('Error: Input set should be a vector')
  n <- length(set)
  if (length(set) < 2)                    stop('Error: Input set must have at least two elements')
  TYPE <- class(set)
  if (!(('numeric' %in% TYPE)|('integer' %in% TYPE)|('character' %in% TYPE))) {
                                          stop('Error: Input set should be a numeric, integer or character vector') }

  #Check input size
  if (!is.vector(size))                   stop('Error: Input size should be a numeric value')
  if (!is.numeric(size))                  stop('Error: Input size should be numeric')
  if (length(size) != 1)                  stop('Error: Input size should be a single numeric value')
  if (size != as.integer(size))           stop('Error: Input size should be a non-negative integer')
  if (size == 0) { return(NULL) }
  if (size < 0)                           stop('Error: Input size should be a non-negative integer')

  #Determine type of input set and generate output
  TYPE <- class(set)
  if (('numeric' %in% TYPE)|('integer' %in% TYPE)) {
    OUT <- matrix(0,  nrow = size, ncol = n)
  } else {
    OUT <- matrix('', nrow = size, ncol = n) }
  rownames(OUT) <- sprintf('D[%s]', 1:size)
  colnames(OUT) <- set[1:n]

  #Generate derangements
  for (k in 1:size) {
    PERM  <- sample.int(n, size = n, replace = FALSE)
    FIXED <- sum(PERM == 1:n)
    while (FIXED > 0) {
      i <- which(PERM == 1:n)[1]
      j <- (1:n)[-i][sample.int(n-1, size = 1)]
      SWAP    <- PERM[j]
      PERM[j] <- PERM[i]
      PERM[i] <- SWAP
      FIXED <- sum(PERM == 1:n) }
    OUT[k, ] <- set[PERM] }

  #Return output
  OUT[1:k, ] }

