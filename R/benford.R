#' Benford Analysis
#'
#' \code{benford} returns analysis of the leading digits of a set of numbers
#'
#' It is known that under broad conditions, the leading digits in a set of positive values should follow
#' Benford's distribution.  This tends to occur when those values have a smooth distribution spanning
#' several orders of magnitude, though it also occurs in some special cases described in the statistical
#' literature on the topic.  In order to identify unusual numeric patterns one can test the frequencies
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
#' the leading digits, another frequency table showing the frequency and proportion of combinations of
#' leading digits and orders-of-magnitude, and Pearson's chi-squared test against Benford's distribution.
#'
#' In order to accommodate "second-order tests", "third-order tests" and so on, the function allows the
#' user to first transform the data by taking differences of unique sorted values a given number of times.
#' The output of the function keeps the original data values and each level of this "diffsort" operation,
#' using only the last set of values for the Benford analysis.
#'
#' **Note:** In some cases (e.g., when taking a high-order test or using data with extreme values) the values
#' used for the test may hit the limits of the lower and upper bounds of representable floating-point numbers
#' in the machine.  If this occurs it can lead to serious inaccuracy in the test due to miscalculation of the
#' leading digits.  In this case the print and plot output for the analysis will give a warning to the user
#' recommending that they change the minimum/maximum representable numbers in the machine.  It is strongly
#' recommended that you do not rely on any analysis that gives either of these warnings.
#'
#' @usage \code{benford(x, base = 10, no.digits = 1, diffsort = 0, conf.level = 0.95, remove.duplicates = FALSE, simulate.p.value = TRUE, simulations = 10000)}
#' @param x A vector of values to analyse
#' @param base A positive integer representing the base for analysis
#' @param no.digits A positive integer representing the number of leading digits for analysis
#' @param diffsort A non-negative integer giving the number of times (if any) to take the sorted-difference of the data
#' @param conf.level The confidence level for confidence intervals for the probability of leading digits
#' @param remove.duplicates Logical; if \code{TRUE} then duplicate values are removed before testing
#' @param simulate.p.value Logical value; if \code{TRUE} the p-value for the chi-squared test is simulated if any frequencies are below five
#' @param simulations A positive integer representing the number of simulations for the p-value if simulation is used
#' @return A list of class \code{benford} containing results of the analysis

benford <- function(x, base = 10, no.digits = 1, diffsort = 0,
                    conf.level = 0.95, remove.duplicates = FALSE,
                    simulate.p.value = TRUE, simulations = 10000) {

  #Get data name
  DATA.NAME <- deparse(substitute(x))

  #Check input x
  if (!is.vector(x))                     stop('Error: Input x should be a numeric vector')
  if (!is.numeric(x))                    stop('Error: Input x should be a numeric vector')
  if (min(x) <= 0)                       stop('Error: Input x should contain only positive values')
  VALUES <- x

  #Check input base
  if (!is.vector(base))                  stop('Error: Input base should be an integer')
  if (!is.numeric(base))                 stop('Error: Input base should be an integer')
  if (length(base) != 1)                 stop('Error: Input base should be an integer')
  BASE <- floor(base)
  if (base != BASE)                      stop('Error: Input base should be an integer')
  if (base < 2)                          stop('Error: Input base should be at least two')

  #Check input no.digits
  if (!is.vector(no.digits))             stop('Error: Input no.digits should be an integer')
  if (!is.numeric(no.digits))            stop('Error: Input no.digits should be an integer')
  if (length(no.digits) != 1)            stop('Error: Input no.digits should be an integer')
  DD <- floor(no.digits)
  if (no.digits != DD)                   stop('Error: Input no.digits should be an integer')
  if (no.digits < 1)                     stop('Error: Input no.digits should be at least one')

  #Check input diffsort
  if (!is.vector(diffsort))              stop('Error: Input diffsort should be an integer')
  if (!is.numeric(diffsort))             stop('Error: Input diffsort should be an integer')
  if (length(diffsort) != 1)             stop('Error: Input diffsort should be an integer')
  DIFFSORT <- floor(diffsort)
  if (diffsort != DIFFSORT)              stop('Error: Input diffsort should be an integer')
  if (diffsort < 0)                      stop('Error: Input diffsort should be non-negative')

  #Check input conf.level
  if (!is.vector(conf.level))            stop('Error: Input conf.level should be a number')
  if (!is.numeric(conf.level))           stop('Error: Input conf.level should be a number')
  if (length(conf.level) != 1)           stop('Error: Input conf.level should be a single number')
  if (min(conf.level) < 0)               stop('Error: Input conf.level should be between zero and one')
  if (max(conf.level) > 1)               stop('Error: Input conf.level should be between zero and one')

  #Check input remove.duplicates
  if (!is.vector(remove.duplicates))     stop('Error: Input remove.duplicates should be a logical value')
  if (!is.logical(remove.duplicates))    stop('Error: Input remove.duplicates should be a logical value')
  if (length(remove.duplicates) != 1)    stop('Error: Input remove.duplicates should be a single logical value')

  #Check input simulate.p.value
  if (!is.vector(simulate.p.value))      stop('Error: Input simulate.p.value should be a logical value')
  if (!is.logical(simulate.p.value))     stop('Error: Input simulate.p.value should be a logical value')
  if (length(simulate.p.value) != 1)     stop('Error: Input simulate.p.value should be a single logical value')

  #Check input simulations
  if (!is.vector(simulations))           stop('Error: Input simulations should be an integer')
  if (!is.numeric(simulations))          stop('Error: Input simulations should be an integer')
  if (length(simulations) != 1)          stop('Error: Input simulations should be an integer')
  SIMS <- floor(simulations)
  if (simulations != SIMS)               stop('Error: Input simulations should be an integer')
  if (simulations < 2000)                stop('Error: Input simulations should be at least 2000')

  #Check if required packages are installed
  if (!requireNamespace('stat.extend', quietly = TRUE)) {
    stop('Error: This function requires the stat.extend package') }

  #############################################################################################################

  #Apply diffsort (if required)
  DATA.INITIAL <- vector(diffsort + 1, mode= 'list')
  DATA.INITIAL[[1]] <- VALUES
  names(DATA.INITIAL)[[1]] <- 'Data'
  if (diffsort > 0) {
    names(DATA.INITIAL)[1+1:diffsort] <- sprintf('Diffsort[%s]', 1:diffsort)
    for (i in 1:diffsort) {
      DATA.INITIAL[[i+1]] <- diff(sort(unique(DATA.INITIAL[[i]])), decreasing = TRUE) } }
  VALUES <- DATA.INITIAL[[diffsort+1]]

  #Determine floating-point warning
  WARNING.MIN <- FALSE
  WARNING.MAX <- FALSE
  if (min(VALUES) == .Machine$double.eps)  { WARNING.MIN <- TRUE }
  if (max(VALUES) == .Machine$double.xmax) { WARNING.MAX <- TRUE }

  #Generate Data Table
  ORDER <- floor(log(VALUES, base = BASE))
  DATA  <- data.frame(Values = VALUES, Remove = FALSE, Base = BASE, Order = ORDER)
  if (remove.duplicates) { DATA$Remove <- duplicated(VALUES) }

  #Add leading digits
  DATA$Digit  <- ifelse(VALUES > 0, floor(VALUES/(BASE^(ORDER))), 0)
  DIGIT.NAMES <- sprintf('D[%s]', 1:DD)
  names(DATA)[5] <- DIGIT.NAMES[1]
  if (DD > 1) {
  for (k in 2:DD) {
    VALUES   <- (VALUES - DATA[, 3+k]*(BASE^(ORDER)))
    ORDER    <- ORDER-1
    DATA$new <- ifelse(VALUES > 0, floor(VALUES/(BASE^(ORDER))), 0)
    names(DATA)[4+k] <- DIGIT.NAMES[k] } }
  DATA$Digit.String <- ''
  for (i in 1:nrow(DATA)) { DATA$Digit.String[i] <- paste0(DATA[i, 5:(4+DD)], collapse = '-') }

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
  FTABLE <- data.frame(Digits = VALSTRING, Frequency = 0, Benford.Freq = 0,
                       Lower = 0, Upper = 1,
                       Proportion = 0, Benford.Probs = 0)
  N <- sum(!DATA$Remove)
  for (i in 1:nrow(WWW)) {
    FTABLE$Frequency[i]     <- sum(DATA$Digit.String[!DATA$Remove] == VALSTRING[i])
    FTABLE$Benford.Probs[i] <- log((1+1/VALNUMBER[i]), base = BASE) }
  FTABLE$Proportion   <- FTABLE$Frequency/N
  FTABLE$Benford.Freq <- N*FTABLE$Benford.Probs
  FTABLE$Digits       <- factor(FTABLE$Digits, levels = FTABLE$Digits)
  FTABLE$Pearson      <- (FTABLE$Frequency - FTABLE$Benford.Freq)^2/FTABLE$Benford.Freq
  if (DD == 1) { names(FTABLE)[1] <- 'Digit' }

  #Add confidence intervals (using Wilson score method)
  for (i in 1:nrow(FTABLE)) {
    CONF.INT <- stat.extend::CONF.prop(alpha = 1-conf.level, n = N, sample.prop = FTABLE$Proportion[i])
    FTABLE$Lower[i] <- N*min(CONF.INT)
    FTABLE$Upper[i] <- N*max(CONF.INT) }
  names(FTABLE)[4] <- paste0('Lower[', conf.level, ']')
  names(FTABLE)[5] <- paste0('Upper[', conf.level, ']')

  #Generate Benford data (full)
  TABLE.ORD  <- table(ORDER)
  ORDER.ALL  <- as.integer(names(TABLE.ORD))
  ORDER.FREQ <- unname(c(TABLE.ORD))
  FTABLE2 <- data.frame(Digits = rep(VALSTRING, times = length(ORDER.ALL)),
                        Order  = rep(ORDER.ALL, each  = length(VALSTRING)),
                        Number = rep(VALNUMBER, times = length(ORDER.ALL)),
                        Total  = rep(ORDER.FREQ, each = length(VALSTRING)),
                        Frequency = 0, Benford.Freq = 0,
                        Lower = 0, Upper = 1,
                        Proportion = 0, Benford.Probs = 0)
  for (i in 1:nrow(FTABLE2)) {
    DIG <- FTABLE2$Digits[i]
    ORD <- FTABLE2$Order[i]
    NUM <- FTABLE2$Number[i]
    TOT <- FTABLE2$Total[i]
    FTABLE2$Frequency[i]     <- sum((DATA$Digit.String[!DATA$Remove] == DIG)*(DATA$Order[!DATA$Remove] == ORD))
    FTABLE2$Benford.Probs[i] <- log((1+1/NUM), base = BASE)
    FTABLE2$Proportion[i]    <- FTABLE2$Frequency[i]/TOT
    FTABLE2$Benford.Freq[i]  <- TOT*FTABLE2$Benford.Probs[i] }
  FTABLE2$Digits  <- factor(FTABLE2$Digits, levels = VALSTRING)
  FTABLE2$Pearson <- (FTABLE2$Frequency - FTABLE2$Benford.Freq)^2/FTABLE2$Benford.Freq
  if (no.digits == 1) { names(FTABLE2)[1] <- 'Digit' }

  #Add confidence intervals (using Wilson score method)
  for (i in 1:nrow(FTABLE2)) {
    TOT <- FTABLE2$Total[i]
    CONF.INT <- stat.extend::CONF.prop(alpha = 1-conf.level, n = TOT,
                                       sample.prop = FTABLE2$Proportion[i])
    FTABLE2$Lower[i] <- TOT*min(CONF.INT)
    FTABLE2$Upper[i] <- TOT*max(CONF.INT) }
  names(FTABLE2)[7] <- paste0('Lower[', conf.level, ']')
  names(FTABLE2)[8] <- paste0('Upper[', conf.level, ']')
  FTABLE2 <- FTABLE2[, -(3:4)]

  #Generate output
  PTABLE <- FTABLE[order(FTABLE$Pearson, decreasing = TRUE), c(1, 8)]
  rownames(PTABLE) <- 1:nrow(PTABLE)

  #Run chi-squared test
  SIMULATE <- FALSE
  if (simulate.p.value) {
  if (sum(FTABLE$Frequency < 5) > 0) {
    SIMULATE <- TRUE } }
  TEST <- suppressWarnings(chisq.test(x = FTABLE$Frequency, p = FTABLE$Benford.Probs,
                                      simulate.p.value = SIMULATE, B = SIMS))

  #Generate output
  OUT  <- list(parameters = list(base = BASE, no.digits = DD, diffsort = diffsort,
                                 remove.duplicates = remove.duplicates, conf.level = conf.level),
               data = list(data.name = DATA.NAME, data.input = DATA.INITIAL),
               data.table = DATA,
               frequency.table = FTABLE,
               frequency.table.full = FTABLE2,
               pearson.table = PTABLE,
               test = list(chisq.stat = c(TEST$statistic), p.value = c(TEST$p.value),
                           simulated = SIMULATE, simulations = simulations),
               warnings = list(warning.min = WARNING.MIN, warning.max = WARNING.MAX))
  class(OUT) <- 'benford'

  #Return output
  OUT }


benford.data <- function(object) {

  #Check object class
  if (!('benford' %in% class(object)))       stop('Error: This plot method is for objects of class \'benford\'')

  #Extract information
  DATA   <- object$data.table
  REMOVE <- DATA$Remove
  MM     <- ncol(DATA)
  BASE   <- object$parameters$base

  #Generate output
  OUT <- DATA[!REMOVE, c(1, 3:(MM-1))]
  attributes(OUT)$base <- BASE

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
  NAME     <- object$data$data.name
  N        <- sum(!object$data.table$Remove)
  RD       <- object$parameters$remove.duplicates
  DD       <- object$parameters$no.digits
  BASE     <- object$parameters$base
  DS       <- object$parameters$diffsort
  STAT     <- object$test$chisq.stat
  PVAL     <- object$test$p.value
  SIMULATE <- object$test$simulated
  SIMS     <- object$test$simulations
  FTABLE   <- object$frequency.table
  PTABLE   <- object$pearson.table

  #Give warnings
  WARNING.MIN <- object$warnings$warning.min
  WARNING.MAX <- object$warnings$warning.max
  if (WARNING.MIN)  {
    warning(paste0('The smallest value used in this test is at the smallest representable floating-point number\n',
                   '  This can lead to serious inaccuracy in the test due to miscalculation of the leading digits\n',
                   '  Recommendation: Reduce .machine$double.eps and recompute to allow accurate calculation of leading digits')) }
  if (WARNING.MAX) {
    warning(paste0('The largest value used in this test is at the largest representable floating-point number\n',
                   '  This can lead to serious inaccuracy in the test due to miscalculation of the leading digits\n',
                   '  Recommendation: Increase .machine$double.xmax and recompute to allow accurate calculation of leading digits')) }

  #Print title
  cat('\n  Benford Analysis\n\n')
  if (DS == 0) {
    cat(paste0('Dataset: ', NAME, '\n'))
  } else {
    cat(paste0('Dataset: ', NAME, ' (diffsort = ', DS, ')\n')) }
  if (DD == 1) {
    cat(paste0('Analysis of ', N, ' values using base ', BASE, '\n'))
  } else {
    cat(paste0('Analysis of ', N, ' values using ', DD, ' digits with base ', BASE, '\n')) }
  cat('Null hypothesis: Digits follow Benford\'s distribution\n')
  cat(paste0('Pearson statistic = ', format(round(STAT, 4), nsmall = 4),
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

  #Print frequency table
  cat('------------------------------------------------------\n')
  cat(' Frequency Table\n\n')
  print(FTABLE[, 1:5], row.names = FALSE, max = 5*max.rows)
  cat('\n')

  #Print Pearson table
  PTABLE.PRINT <- PTABLE
  PTABLE.PRINT$a <- ''
  PTABLE.PRINT <- PTABLE.PRINT[, c(1, 3, 2)]
  names(PTABLE.PRINT)[2] <- ' '
  cat('------------------------------------------------------\n')
  cat(' Pearson Table\n\n')
  print(PTABLE.PRINT, row.names = FALSE, max = 5*max.rows)
  cat('\n') }


plot.benford <- function(object, conf.level = NULL, by.order = FALSE, order.legend = FALSE,
                         colour.density = 'blue', colour.bar = 'blue',
                         colour.point = 'green', colour.pearson = 'darkgreen') {

  #Check object class
  if (!('benford' %in% class(object)))       stop('Error: This plot method is for objects of class \'benford\'')

  #Check conf.level
  if (!is.null(conf.level)) {
    if (!is.vector(conf.level))              stop('Error: Input conf.level should be NULL or a numeric value')
    if (!is.numeric(conf.level))             stop('Error: Input conf.level should be NULL or a numeric value')
    if (length(conf.level) != 1)             stop('Error: Input conf.level should be NULL or a single numeric value')
    if (min(conf.level) <= 0)                stop('Error: Input conf.level should be between zero and one')
    if (max(conf.level) >= 1)                stop('Error: Input conf.level should be between zero and one') }

  #Check input by.order
  if (!is.vector(by.order))                  stop('Error: Input by.order should be a vector')
  if (!is.logical(by.order))                 stop('Error: Input by.order should be a logical value')
  if (length(by.order) != 1)                 stop('Error: Input by.order should be a single logical value')

  #Check input order.legend
  if (!is.vector(order.legend))              stop('Error: Input order.legend should be a vector')
  if (!is.logical(order.legend))             stop('Error: Input order.legend should be a logical value')
  if (length(order.legend) != 1)             stop('Error: Input order.legend should be a single logical value')

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
  N        <- sum(!object$data.table$Remove)
  RD       <- object$parameters$remove.duplicates
  DD       <- object$parameters$no.digits
  BASE     <- object$parameters$base
  CL       <- object$parameters$conf.level
  STAT     <- object$test$chisq.stat
  PVAL     <- object$test$p.value
  SIMULATE <- object$test$simulated
  SIMS     <- object$test$simulations
  DATA     <- object$data.table
  FTABLE   <- object$frequency.table
  FTABLE2  <- object$frequency.table.full
  PTABLE   <- object$pearson.table

  #Give warnings
  WARNING.MIN <- object$warnings$warning.min
  WARNING.MAX <- object$warnings$warning.max
  if (WARNING.MIN)  {
    warning(paste0('The smallest value used in this test is at the smallest representable floating-point number\n',
                   '  This can lead to serious inaccuracy in the test due to miscalculation of the leading digits\n',
                   '  Recommendation: Reduce .machine$double.eps and recompute to allow accurate calculation of leading digits')) }
  if (WARNING.MAX) {
    warning(paste0('The largest value used in this test is at the largest representable floating-point number\n',
                   '  This can lead to serious inaccuracy in the test due to miscalculation of the leading digits\n',
                   '  Recommendation: Increase .machine$double.xmax and recompute to allow accurate calculation of leading digits')) }

  #Recalculate confidence intervals (if conf.level is specified)
  if (!is.null(conf.level)) {
    if (conf.level != CL) {

      #Remove previous confidence intervals and add new confidence intervals (using Wilson score method)
      names(FTABLE)[4] <- paste0('Lower[', conf.level, ']')
      names(FTABLE)[5] <- paste0('Upper[', conf.level, ']')
      for (i in 1:nrow(FTABLE)) {
        CONF.INT <- stat.extend::CONF.prop(alpha = 1-conf.level, n = N, sample.prop = FTABLE$Proportion[i])
        FTABLE[i, 4] <- N*min(CONF.INT)
        FTABLE[i, 5] <- N*max(CONF.INT) }
      CL <- conf.level } }

  #Set plot title and theme
  TITLE <- paste0('Benford Plot')
  if (DD == 1) {
    SUBTITLE <- paste0('Analysis of ', N, ' values with base ', BASE, '\n',
                       'Pearson Statistic = ', format(round(STAT, 4), nsmall = 4), ' (', nrow(FTABLE)-1, ' DF), p-value = ', format(round(PVAL, 4), nsmall = 4), '\n',
                       '(Points show expected frequency; Error bars show ',
                       format(100*CL, nsmall = 0), '% confidence intervals)')
  } else {
    SUBTITLE <- paste0('Analysis of ', N, ' values using ', DD, ifelse(DD == 1, ' digit ', ' digits '),
                       'with base ', BASE, '\n',
                       'Pearson Statistic = ', format(round(STAT, 4), nsmall = 4), ' (', nrow(FTABLE)-1, ' DF), p-value = ', format(round(PVAL, 4), nsmall = 4), '\n',
                       '(Points show expected frequency; Error bars show ',
                       format(100*CL, nsmall = 0), '% confidence intervals)') }
  TITLE    <- grid::textGrob(TITLE,    gp = grid::gpar(fontsize = 14, fontface = 'bold'))
  SUBTITLE <- grid::textGrob(SUBTITLE, gp = grid::gpar(fontsize = 10))
  THEME <- ggplot2::theme(plot.title    = ggplot2::element_text(hjust = 0.5, size = 14, face = 'bold'),
                          plot.subtitle = ggplot2::element_text(hjust = 0.5, face = 'bold'),
                          axis.title.x  = ggplot2::element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                          axis.title.y  = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)))

  #Set label function for density axis
  #This function is a variation of scales::math_format which substitutes .base = BASE
  label_format <- function (expr = .base^.x, format = force) {
    .x <- NULL
    quoted <- substitute(expr)
    subs <- function(x) { do.call("substitute", list(quoted, list(.x = x, .base = BASE))) }
    function(x) {
      x <- format(x)
      ret <- lapply(x, subs)
      ret <- as.expression(ret)
      ret[is.na(x)] <- NA
      names(ret) <- names(x)
      ret } }

  #Preliminary changes for figures
  LOGB <- function(x) { log(x, base = BASE) }
  names(FTABLE)[c(1, 4:5)]  <- c('Digits', 'Lower', 'Upper')
  names(FTABLE2)[c(1, 5:6)] <- c('Digits', 'Lower', 'Upper')
  xxx <- seq(from = 1, to = (BASE-1)*BASE^(DD-1), by = BASE^(DD-1))
  LABEL.DIGITS <- rep('', nrow(FTABLE))
  LABEL.DIGITS[xxx] <- as.character(FTABLE[xxx, 1])

  #Generate Figure 1 - Density Plot (logarithmic scale)
  FIGURE1 <- ggplot2::ggplot(ggplot2::aes(x = Values), data = DATA[!DATA$Remove, ]) +
             ggplot2::geom_density(fill = colour.density) +
             ggplot2::scale_x_log10(breaks = scales::trans_breaks("LOGB", function(x, base = BASE) base^x),
                                    labels = scales::trans_format("LOGB", label_format(.base^.x))) +
             THEME + ggplot2::xlab('Value') + ggplot2::ylab('Density')

  #Generate Figure 2 - Frequency Plot
  if (!by.order) {
    FIGURE2 <- ggplot2::ggplot(ggplot2::aes(x = Digits, y = Frequency), data = FTABLE) +
      ggplot2::geom_bar(stat = 'identity', fill = colour.bar, position = ggplot2::position_dodge()) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper), width = 0.2,
                             position = ggplot2::position_dodge(0.9)) +
      ggplot2::geom_point(ggplot2::aes(y = Benford.Freq),
                          stat = 'identity', colour = colour.point, size = 3) +
      ggplot2::scale_x_discrete(labels = LABEL.DIGITS) +
      THEME + ggplot2::xlab(ifelse(DD == 1, 'Leading Digit', 'Leading Digits')) + ggplot2::ylab('Frequency')
  } else {
    FIGURE2 <- ggplot2::ggplot(ggplot2::aes(x = Digits, y = Frequency, fill = Order), data = FTABLE2) +
      ggplot2::geom_bar(stat = 'identity', position = 'stack') +
      ggplot2::scale_x_discrete(labels = LABEL.DIGITS) +
      THEME + ggplot2::xlab(ifelse(DD == 1, 'Leading Digit', 'Leading Digits')) + ggplot2::ylab('Frequency')
    FIGURE2 <- FIGURE2 + ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper, fill = NULL), width = 0.2,
                                                position = ggplot2::position_dodge(0.9), data = FTABLE) +
                         ggplot2::geom_point(mapping = ggplot2::aes(x = Digits, y = Benford.Freq, fill = NULL),
                                             data = FTABLE, stat = 'identity', colour = colour.point, size = 3) +
                         ggplot2::scale_fill_gradient2(low = 'black', mid = colour.bar, high = 'lightgrey') +
                         ggplot2::theme(legend.position = 'bottom') +
                         ggplot2::labs(fill = 'Order-of-Magnitude\n') +
      { if (!order.legend) ggplot2::theme(legend.position = 'none') } }

  #Generate Figure 3 - Pearson Statistic
  names(PTABLE)[1] <- 'Digits'
  FIGURE3 <- ggplot2::ggplot(ggplot2::aes(x = Digits, y = Pearson), data = PTABLE) +
             ggplot2::geom_bar(stat = 'identity', fill = colour.pearson, position = ggplot2::position_dodge()) +
             ggplot2::scale_y_continuous(limits = c(0, 1.2*max(PTABLE$Pearson))) +
             THEME +
             ggplot2::xlab(ifelse(DD == 1, 'Leading Digit', 'Leading Digits')) +
             ggplot2::ylab('Pearson Statistic')

  #Combine figures into single plot
  #Line up vertical axes in combined plot
  G1 <- ggplot2::ggplotGrob(FIGURE1)
  G2 <- ggplot2::ggplotGrob(FIGURE2)
  G3 <- ggplot2::ggplotGrob(FIGURE3)
  MAXWIDTH <- grid::unit.pmax(G1$widths[2:5], G2$widths[2:5], G3$widths[2:5])
  G1$widths[2:5] <- as.list(MAXWIDTH)
  G2$widths[2:5] <- as.list(MAXWIDTH)
  G3$widths[2:5] <- as.list(MAXWIDTH)
  FIGURES <- gridExtra::grid.arrange(TITLE, SUBTITLE, G1, G2, G3, ncol = 1,
                                     heights = c(1.4, 4, 12, 12+4*by.order*order.legend, 12))

  #Return plot
  FIGURES }

