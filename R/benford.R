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
#' @usage \code{benford(x, base = 10, no.digits = 1, simulate.p.value = TRUE, simulations = 10000)}
#' @param x A vector of values to analyse
#' @param base A positive integer representing the base for analysis
#' @param no.digits A positive integer representing the number of leading digits for analysis
#' @param simulate.p.value Logical value; if \code{TRUE} the p-value for the chi-squared test is simulated if any frequencies are below five
#' @param simulations A positive integer representing the number of simulations for the p-value if simulation is used
#' @return A list of class \code{benford} containing results of the analysis

benford <- function(x, base = 10, no.digits = 1, simulate.p.value = TRUE, simulations = 10000) {

  #Check input x
  if (!is.vector(x))                     stop('Error: Input x should be a numeric vector')
  if (!is.numeric(x))                    stop('Error: Input x should be a numeric vector')
  if (min(x) <= 0)                       stop('Error: Input x should contain only positive values')
  VALUES <- x

  #Check input base
  if (!is.vector(base))                  stop('Error: Input base should be an integer')
  if (!is.numeric(base))                 stop('Error: Input base should be an integer')
  if (length(base) != 1)                 stop('Error: Input base should be an integer')
  BASE <- as.integer(base)
  if (base != BASE)                      stop('Error: Input base should be an integer')
  if (base < 2)                          stop('Error: Input base should be at least two')

  #Check input no.digits
  if (!is.vector(no.digits))             stop('Error: Input no.digits should be an integer')
  if (!is.numeric(no.digits))            stop('Error: Input no.digits should be an integer')
  if (length(no.digits) != 1)            stop('Error: Input no.digits should be an integer')
  DD <- as.integer(no.digits)
  if (no.digits != DD)                   stop('Error: Input no.digits should be an integer')
  if (no.digits < 1)                     stop('Error: Input no.digits should be at least one')

  #Check input simulate.p.value
  if (!is.vector(simulate.p.value))      stop('Error: Input simulate.p.value should be a logical value')
  if (!is.logical(simulate.p.value))     stop('Error: Input simulate.p.value should be a logical value')
  if (length(simulate.p.value) != 1)     stop('Error: Input simulate.p.value should be a single logical value')

  #Check input simulations
  if (!is.vector(simulations))           stop('Error: Input simulations should be an integer')
  if (!is.numeric(simulations))          stop('Error: Input simulations should be an integer')
  if (length(simulations) != 1)          stop('Error: Input simulations should be an integer')
  SIMS <- as.integer(simulations)
  if (simulations != SIMS)               stop('Error: Input simulations should be an integer')
  if (simulations < 2000)                stop('Error: Input simulations should be at least 2000')

  #Generate Data Table
  ORDER <- floor(log(VALUES, base = BASE))
  DATA  <- data.frame(Values = VALUES, Order = ORDER)

  #Add leading digits
  DATA$Digit  <- ifelse(VALUES > 0, floor(VALUES/(BASE^(ORDER))), 0)
  DIGIT.NAMES <- sprintf('Digit[%s]', 1:DD)
  names(DATA)[3] <- DIGIT.NAMES[1]
  if (DD > 1) {
  for (k in 2:DD) {
    VALUES   <- (VALUES - DATA[, 1+k]*(BASE^(ORDER)))
    ORDER    <- ORDER-1
    DATA$new <- ifelse(VALUES > 0, floor(VALUES/(BASE^(ORDER))), 0)
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
  FTABLE <- data.frame(Digits = VALSTRING, Frequency = 0, Benford.Freq = 0, Proportion = 0, Benford.Probs = 0)
  N <- length(VALUES)
  for (i in 1:nrow(WWW)) {
    FTABLE$Frequency[i]     <- sum(DATA$Digit.String == VALSTRING[i])
    FTABLE$Benford.Probs[i] <- log((1+1/VALNUMBER[i]), base = BASE) }
  FTABLE$Proportion   <- FTABLE$Frequency/N
  FTABLE$Benford.Freq <- N*FTABLE$Benford.Probs
  FTABLE$Digits <- factor(FTABLE$Digits, levels = FTABLE$Digits)

  #Run chi-squared test
  SIMULATE <- FALSE
  if (simulate.p.value) {
  if (sum(DATA2$Frequency < 5) > 0) {
    SIMULATE <- TRUE } }
  TEST <- suppressWarnings(chisq.test(x = FTABLE$Frequency, p = FTABLE$Benford.Probs,
                                      simulate.p.value = SIMULATE, B = SIMS))

  #Generate output
  OUT  <- list(base = BASE, no.digits = DD, data.table = DATA, frequency.table = FTABLE,
               chisq.stat = c(TEST$statistic), p.value = c(TEST$p.value),
               simulated = SIMULATE, simulations = simulations)
  class(OUT) <- 'benford'

  #Return output
  OUT }


print.benford <- function(object, max.rows = 20) {

  #Check object class
  if (!('benford' %in% class(object)))       stop('Error: This print method is for objects of class \'benford\'')

  #Check max.rows
  if (!is.vector(max.rows))                  stop('Error: Input max.rows should be a vector')
  if (!is.numeric(max.rows))                 stop('Error: Input max.rows should be a numeric vector')
  if (length(max.rows) != 1)                 stop('Error: Input max.rows should be a single integer value')
  MAX <- as.integer(max.rows)
  if (max.rows != MAX)                       stop('Error: Input max.rows should be an integer')
  if (max.rows < 1)                          stop('Error: Input max.rows should be a positive integer')

  #Extract information
  N        <- nrow(object$data.table)
  DD       <- object$no.digits
  BASE     <- object$base
  STAT     <- object$chisq.stat
  PVAL     <- object$p.value
  FTABLE   <- object$frequency.table
  SIMULATE <- object$simulated
  SIMS     <- object$simulations

  #Print title
  cat('\n  Benford Analysis\n\n')
  if (DD == 1) {
    cat(paste0('Analysis of ', N, ' values using base = ', BASE, '\n'))
  } else {
    cat(paste0('Analysis of ', N, ' values using base = ', BASE, ' and digits = ', DD, '\n')) }
  cat('Null hypothesis: Digits follow Benford\'s distribution\n')
  cat(paste0('Chi-Sq statistic = ', format(round(STAT, 4), nsmall = 4),
             ' (', nrow(FTABLE)-1, ' DF), p-value = ', format(round(PVAL, 4), nsmall = 4), '\n'))
  if (SIMULATE) {
    cat('Test used simulated p-value with', format(SIMS, nsmall = 0, big.mark = ','), 'simulations\n')
  } else {
    cat('Test used chi-squared distribution to compute p-value\n')
    SS <- sum(FTABLE$Frequency < 5)
    if (SS == 1) {
      cat(paste0('\n--- There was ', SS, ' frequency below five\n',
                 '--- Consider switching to simulated p-value instead\n')) }
    if (SS > 1) {
      cat(paste0('\n--- There were ', SS, ' frequencies below five\n',
                 '--- Consider switching to simulated p-value instead\n')) } }
  cat('\n')

  #Print table
  cat('-------------------------------------------------------\n')
  print(FTABLE, row.names = FALSE, max = 5*max.rows)
  cat('\n') }


plot.benford <- function(object, conf.level = 0.95,
                         colour.density = 'blue', colour.bar = 'blue', colour.point = 'green') {

  #Check object class
  if (!('benford' %in% class(object)))       stop('Error: This plot method is for objects of class \'benford\'')

  #Check if required packages are installed
  GGPLOT2   <- requireNamespace('ggplot2',   quietly = TRUE)
  SCALES    <- requireNamespace('scales',    quietly = TRUE)
  GRID      <- requireNamespace('grid',      quietly = TRUE)
  GRIDEXTRA <- requireNamespace('gridExtra', quietly = TRUE)
  if (!GGPLOT2)    stop('Error: This plot method requires the ggplot2 package')
  if (!SCALES)     stop('Error: This plot method requires the scales package')
  if (!GRIDEXTRA)  stop('Error: This plot method requires the gridExtra package')
  if (!GRID)       stop('Error: This plot method requires the grid package')

  #Extract information
  N        <- nrow(object$data.table)
  DD       <- object$no.digits
  BASE     <- object$base
  STAT     <- object$chisq.stat
  PVAL     <- object$p.value
  DATA     <- object$data.table
  FTABLE   <- object$frequency.table
  SIMULATE <- object$simulated
  SIMS     <- object$simulations

  #Add bound for error bars
  FTABLE$Lower <- 0
  FTABLE$Upper <- 1
  for (i in 1:nrow(FTABLE)) {
    CONF.INT <- stat.extend::CONF.prop(alpha = 1-conf.level, n = N, sample.prop = FTABLE$Proportion[i])
    FTABLE$Lower[i] <- N*min(CONF.INT)
    FTABLE$Upper[i] <- N*max(CONF.INT) }

  #Set plot title and theme
  TITLE <- paste0('Benford Plot (Base ', BASE, ')')
  SUBTITLE <- paste0('Chi-Sq statistic = ', format(round(STAT, 4), nsmall = 4), ' (', nrow(FTABLE)-1, ' DF), p-value = ', format(round(PVAL, 4), nsmall = 4), '\n',
                     '(Points show expected frequency; Error bars show ',
                     format(100*conf.level, nsmall = 0), '% confidence intervals)')
  TITLE    <- grid::textGrob(TITLE,    gp = grid::gpar(fontsize = 14, fontface = 'bold'))
  SUBTITLE <- grid::textGrob(SUBTITLE, gp = grid::gpar(fontsize = 10))
  THEME <- ggplot2::theme(plot.title    = ggplot2::element_text(hjust = 0.5, size = 14, face = 'bold'),
                          plot.subtitle = ggplot2::element_text(hjust = 0.5, face = 'bold'),
                          axis.title.x  = ggplot2::element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                          axis.title.y  = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)))

  #Create figures
  #UNRESOLVED PROBLEM: The labeling for the x-axis in Figure 2 is not working correctly
  #Presently showing the label with the word "BASE" instead of the value BASE
  LOGB <- function(x) { log(x, base = BASE) }
  FIGURE1 <- ggplot2::ggplot(ggplot2::aes(x = log(10, base = BASE)*Values), data = DATA) +
               ggplot2::geom_density(fill = colour.density) +
               ggplot2::scale_x_log10(breaks = scales::trans_breaks("LOGB", function(x) BASE^x),
                                      labels = scales::trans_format("LOGB", scales::math_format(BASE^.x))) +
               THEME + ggplot2::xlab('Value') + ggplot2::ylab('Density')
  FIGURE2 <- ggplot2::ggplot(ggplot2::aes(x = Digits, y = Frequency), data = FTABLE) +
               ggplot2::geom_bar(stat = 'identity', fill = colour.bar, position = ggplot2::position_dodge()) +
               ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper), width = 0.2,
                                      position = ggplot2::position_dodge(0.9)) +
               ggplot2::geom_point(ggplot2::aes(y = Benford.Freq),
                                   stat = 'identity', colour = colour.point, size = 3) +
               THEME + ggplot2::xlab(ifelse(DD == 1, 'Leading Digit', 'Leading Digits')) + ggplot2::ylab('Frequency')
  FIGURES <- gridExtra::grid.arrange(TITLE, SUBTITLE, FIGURE1, FIGURE2,
                                     ncol = 1, heights = c(1.4, 2, 14, 14))


  #Return plot
  FIGURES }

