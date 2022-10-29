#' Benford Analysis
#'
#' \code{benford} returns analysis of the leading digits of a set of numbers
#'
#' It is known that under broad conditions, the leading digits in a set of positive values should follow
#' Benford's distribution.  This tends to occur when those values have a smooth distribution spanning
#' several orders of magnitude.  In order to identify unusual numeric patterns one can test the frequencies
#' of one or more leading digits in a set of values against the expected outcome under Benford's distribution.
#' Large deviations from Benford's distribution for values spanning over multiple orders of magnitude can
#' be indicative of unusual numeric patterns warranting further review (e.g., for use in fraud detection).
#'
#' The present function takes a set of values and determines the leading digits of those values under a
#' stipulated base and number of digits for examination.  It then computes the frequency and proportion of
#' every possible combination of leading digits and compares this to Benford's distribution.  The function
#' then computes the chi-squared statistic for comparison to the Benford distribution and computes the
#' corresponding p-value.  The output is a list containing information about the analysis, including a data
#' table with the leading digits of each value, a frequency table showing the frequency and proportion of
#' the leading digits, and the results of the chi-squared test against Benford's distribution.
#'
#' @usage \code{benford(x, base = 10, no.digits = 1)}
#' @param x A vector of values to analyse
#' @param base A positive integer representing the base for analysis
#' @param no.digits A positive integer representing the number of leading digits for analysis
#' @return A list of class \code{benford} containing information on the analysis

benford <- function(x, base = 10, no.digits = 1) {

  #Check input x
  if (!is.vector(x))             stop('Error: Input x should be a numeric vector')
  if (!is.numeric(x))            stop('Error: Input x should be a numeric vector')
  if (min(x) <= 0) {
    warning('Benford analysis is designed for positive values\n---analysis is conducted on absolute values, ignoring zeros, to proceed')
    XX <- abs(x) } else { XX <- x }

  #Check input base
  if (!is.vector(base))             stop('Error: Input base should be an integer')
  if (!is.numeric(base))            stop('Error: Input base should be an integer')
  if (length(base) != 1)            stop('Error: Input base should be an integer')
  BASE <- as.integer(base)
  if (base != BASE)                 stop('Error: Input base should be an integer')
  if (base < 2)                     stop('Error: Input base should be at least two')

  #Check input no.digits
  if (!is.vector(no.digits))        stop('Error: Input no.digits should be an integer')
  if (!is.numeric(no.digits))       stop('Error: Input no.digits should be an integer')
  if (length(no.digits) != 1)       stop('Error: Input no.digits should be an integer')
  DD <- as.integer(no.digits)
  if (no.digits != DD)              stop('Error: Input no.digits should be an integer')
  if (no.digits < 1)                stop('Error: Input no.digits should be at least one')

  #Generate Data Table
  ORDER <- floor(log(XX, base = BASE))
  DATA  <- data.frame(Values = XX, Order = ORDER)

  #Add leading digits
  DATA$Digit  <- ifelse(XX > 0, floor(XX/(BASE^(ORDER))), 0)
  DIGIT.NAMES <- sprintf('Digit[%s]', 1:DD)
  names(DATA)[3] <- DIGIT.NAMES[1]
  if (DD > 1) {
  for (k in 2:DD) {
    XX <- (XX - DATA[, 1+k]*(BASE^(ORDER)))
    ORDER <- ORDER-1
    DATA$new <- ifelse(XX > 0, floor(XX/(BASE^(ORDER))), 0)
    names(DATA)[2+k] <- DIGIT.NAMES[k] } }
  DATA$Digit.String <- ''
  for (i in 1:nrow(DATA)) { DATA$Digit.String[i] <- paste0(DATA[i, 3:(2+DD)], collapse = '-') }

  #Generate digit combinations and labels
  VALS0 <- 0:(BASE-1)
  VALS1 <- 1:(BASE-1)
  WWW <- matrix(0, nrow = (BASE-1)*(BASE^(DD-1)), ncol = DD)
  WWW[, 1] <- rep(VALS1, each = (BASE^(DD-1)))
  if (DD > 1) {
  for (k in 2:DD) {
    WWW[, k] <- rep(rep(VALS0, each = (BASE^(DD-k))), times = (BASE-1)*(BASE^(k-2))) } }
  VALSTRING <- rep('', nrow(WWW))
  VALNUMBER <- rep(0,  nrow(WWW))
  for (i in 1:nrow(WWW)) {
    VALSTRING[i] <- paste0(WWW[i, ], collapse = '-')
    VALNUMBER[i] <- sum(WWW[i, ]*rev(BASE^(0:(DD-1)))) }

  #Generate Benford data
  DATA2 <- data.frame(Digits = VALSTRING, Frequency = 0, Benford.Freq = 0, Proportion = 0, Benford.Probs = 0)
  N <- sum(DATA$Values != 0)
  for (i in 1:nrow(WWW)) {
    DATA2$Frequency[i]     <- sum(DATA$Digit.String == VALSTRING[i])
    DATA2$Benford.Probs[i] <- log((1+1/VALNUMBER[i]), base = BASE) }
  DATA2$Proportion   <- DATA2$Frequency/N
  DATA2$Benford.Freq <- N*DATA2$Benford.Probs

  #Generate output
  TEST <- suppressWarnings(chisq.test(x = DATA2$Frequency, p = DATA2$Benford.Probs))
  OUT  <- list(base = BASE, no.digits = DD, data.table = DATA, frequency.table = DATA2,
               chisq.stat = c(TEST$statistic), p.value = c(TEST$p.value))
  class(OUT) <- 'benford'

  #Return output
  OUT }


print.benford <- function(object) {

  #Check object class
  if (!('benford' %in% class(object)))       stop('Error: This print methd is for objects of class \'benford\'')

  #Extract information
  N      <- nrow(object$data.table)
  DD     <- object$no.digits
  BASE   <- object$base
  STAT   <- object$chisq.stat
  PVAL   <- object$p.value
  FTABLE <- object$frequency.table

  #Print title
  cat('\n  Benford Analysis\n\n')
  if (DD == 1) {
    cat(paste0('Analysis of ', N, ' values using base = ', BASE, '\n'))
  } else {
    cat(paste0('Analysis of ', N, ' values using base = ', BASE, ' and digits = ', DD, '\n')) }
  cat('Null hypothesis: Digits follow Benford\'s distribution\n')
  cat(paste0('Chi-squared statistic = ', format(round(STAT, 4), nsmall = 4),
             ', p-value = ', format(round(PVAL, 4), nsmall = 4), '\n\n'))

  #Print table
  cat('-------------------------------------------------------\n')
  print(FTABLE, row.names = FALSE)
  cat('\n') }

