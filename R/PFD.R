#' Prime Factor Decomposition (PFD)
#'
#' \code{PFD} converts a positive integer to its prime-factor decomposition or *vice versa*
#'
#' This function converts a vector of integers to a corresponding character vector giving the prime-factor decomposition in a condensed form.
#' The input can be a vector of integers or a 'bigz' vector containing large integers.  In either case the function returns the corresponding
#' vector of the prime-factor decomposition (PFD) values, written in a condensed character form.  The function also converts back from the PFD
#' form to an integer/bigz vector.
#'
#' This function depends on the gmp package.
#'
#' @param x An input vector/matrix/array (can be a vector of integers/bigz or PFDs)
#' @return If the input is integer/bigz then the output is the PFD; if the input is PFD then the output is integer/bigz
#'
#' @examples
#' PFD(1:10)
#' stopifnot(all.equal(1:100, PFD(PFD(1:100))))
PFD <- function(x) {

  stopifnot("PFD() depends on the `gmp` package."=requireNamespace("gmp"))

  #Convert PFD to integer/bigz
  if ('prime.factor.decomposition' %in% class(x)) {

    #Check input
    if (!is.character(x))                        { stop('Input x must either be a PFD or a vector of positive integers') }

    #Convert to vector
    xx <- as.vector(x)

    #Create output vector
    n   <- length(xx)
    OUT <- gmp::as.bigz(integer(n))
    for (i in 1:n) {
      if (xx[i] == '|0|') { OUT[i] <- gmp::as.bigz(0) }
      if (xx[i] == '|1|') { OUT[i] <- gmp::as.bigz(1) }
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
    xx <- as.vector(gmp::as.bigz(x))

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
        PRIMES  <- sort(unique(FACTORS))
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

#' @rdname PFD
#' @param quote logical, indicating whether or not strings should be printed with surrounding quotes.
#' @param ... further arguments passed to or from other methods.
print.prime.factor.decomposition <- function (x, quote = FALSE, ...) {

  #Get object types
  n      <- length(x)
  DIM    <- attr(x, 'dim')
  DIMLAB <- paste(DIM, collapse = ' x ')

  if (n > 0) {
    if (is.null(DIM)) {
      cat('Prime Factor Decomposition (\'PFD\') vector of length', length(x), '\n')
      print(as.character(x), quote = quote, ...) } else {
    if (length(DIM) == 1) {
      cat('Prime Factor Decomposition (\'PFD\') matrix of dimensions', DIMLAB, '\n')
      print(array(as.character(x), dim = DIM), quote = quote, ...) }
    if (length(DIM) > 1) {
      cat('Prime Factor Decomposition (\'PFD\') array of dimensions', DIMLAB, '\n')
      print(array(as.character(x), dim = DIM), quote = quote, ...) } }
  } else {
    cat('Prime Factor Decomposition (\'PFD\') with no elements\n') }

  invisible(x) }
