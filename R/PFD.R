#' Prime Factor Decomposition (PFD)
#'
#' \code{PFD} converts a positive integer to its prime-factor decomposition or *vice versa*
#'
#' This function converts a vector of integers to a corresponding character vector giving the prime-factor decomposition in a condensed form.
#' The input can be a vector of integers or a 'bigz' vector containing large integers.  In either case the function returns the corresponding
#' vector of the prime-factor decomposition (PFD) values, written in a condensed character form.  The function also converts back from the PFD
#' form to an integer/bigz vector.
#'
#' @usage \code{PFD}
#' @param x An input vector (can be a vector of integers/bigz or PFDs)
#' @return If the input is an integer/bigz vector, the output is the PFD vector; if the input is a PFD vector, the output is an integer/bigz vector

PFD <- function(x) {

  #Convert PFD to integer/bigz
  if ('prime.factor.decomposition' %in% class(x)) {

    #Check input
    if (!is.character(x))                        { stop('Input x must either be a PFD or a vector of positive integers') }

    #Convert to vector
    xx <- as.vector(x)

    #Create output vector
    n   <- length(xx)
    OUT <- as.bigz(integer(n))
    for (i in 1:n) {
      if (xx[i] == '|0|') { OUT[i] <- as.bigz(0) }
      if (xx[i] == '|1|') { OUT[i] <- as.bigz(1) }
      if (!(xx[i] %in% c('|0|', '|1|'))) {
        TABLE  <- read.table(text = gsub("\\|", "\n", xx[i]), sep = "^", col.names = c("Prime", "Multiplicity"))
        PRIMES <- gmp::as.bigz(TABLE$Prime)
        MULTS  <- gmp::as.bigz(TABLE$Multiplicity)
        OUT[i] <- prod(PRIMES^MULTS) } }

    #Convert to integer vector (if able)
    if (max(abs(OUT)) <= .Machine$integer.max) { OUT <- as.integer(OUT) }

    #Convert to matrix/array
    if (length(dim(x)) > 1) {
      DIM <- dim(x)
      OUT <- array(OUT, dim = DIM) } }

  #Convert integer/bigz to PFD
  if (!('prime.factor.decomposition' %in% class(x))) {

    #Convert to bigz vector
    xx <- as.vector(as.bigz(x))

    #Create output vector
    n   <- length(xx)
    OUT <- character(n)
    class(OUT) <- 'prime.factor.decomposition'
    for (i in 1:n) {
      if (xx[i] == 0) {
        OUT[i]  <- '|0|'
      } else {
      if (xx[i] == 1) {
        OUT[i]  <- '|1|'
      } else {
        FACTORS <- gmp::factorize(xx[i])
        PRIMES  <- unique(FACTORS)
        K       <- length(PRIMES)
        MULTS   <- integer(K)
        STRING  <- character(K)
        for (k in 1:K) {
          MULTS[k]  <- sum(FACTORS == PRIMES[k])
          STRING[k] <- paste0(PRIMES[k], '^', MULTS[k]) }
        OUT[i]  <- paste0('|', paste(STRING, collapse = '|'), '|') } } }

    #Convert to matrix/array
    if (length(dim(x)) > 1) {
      attr(OUT, 'dim') <- dim(x) } }

  #Give output
  OUT }


print.prime.factor.decomposition <- function (OBJECT, quote = FALSE, ...) {

  #Get object types
  n      <- length(OBJECT)
  DIM    <- attr(OBJECT, 'dim')
  DIMLAB <- paste(DIM, collapse = ' x ')

  if (n > 0) {
    if (is.null(DIM)) {
      cat('Prime Factor Decomposition (\'PFD\') vector of length', length(OBJECT), '\n')
      print(as.character(OBJECT), quote = quote, ...) } else {
    if (length(DIM) == 1) {
      cat('Prime Factor Decomposition (\'PFD\') matrix of dimensions', DIMLAB, '\n')
      print(array(as.character(OBJECT), dim = DIM), quote = quote, ...) }
    if (length(DIM) > 1) {
      cat('Prime Factor Decomposition (\'PFD\') array of dimensions', DIMLAB, '\n')
      print(array(as.character(OBJECT), dim = DIM), quote = quote, ...) } }
  } else {
    cat('Prime Factor Decomposition (\'PFD\') with no elements\n') }

  invisible(x) }
