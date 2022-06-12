#' Generate Zipf plot
#'
#' \code{zipfplot} generates the Zipf plot for the input data
#'
#' The Zipf plot for a dataset shows the ranks of outcomes versus their frequency on a log-log scale.
#' It is used to determine how closely a dataset follows "Zipf's law".  The present function takes in
#' a vector of values and produces the Zipf plot.  The data input can be either a vector, a matrix or
#' a data frame.  If the input data is a vector then the output will be a Zipf plot for that data vector.
#' If the input data is a matrix or data frame then each column will be treated as a separate variable
#' and the output will be a single Zipf plot showing each of the variables.  The user can control
#' whether the variables are shown on a single plot or separate plots.
#'
#' @param x Data vector, matrix or data-frame
#' @param relative.freq Logical; if \code{TRUE} the plot shows the relative frequency on vertical axis
#' @param smooth.line Logical; if \code{TRUE} the plot shows a smoothed line through the data using LOESS method
#' @param smooth.conf Logical; if \code{TRUE} the plot shows confidence bands on the smoothed line (only shown if smoothed line is shown)
#' @param conf.level The confidence level for the confidence bands on the smoothed line
#' @param separate.plots Logical; if \code{TRUE} the plot shows
#' @param data.name Logical; if \code{TRUE} the subtitle will state the name of the input data
#' @param point.size Size of the points in the plot
#' @param point.alpha Alpha-transparency of the points in the plot
#' @return Zipf plot for the input data
#' @examples
#' try(zipfplot(sample(LETTERS, 300, replace = TRUE)))

zipfplot <- function(x, relative.freq = TRUE, smooth.line = TRUE, smooth.conf = TRUE,
                     conf.level = 0.99, separate.plots = FALSE,
                     data.name = FALSE, point.size = 3, point.alpha = 0.4) {

  #Check input x
  DATA.NAME <- deparse(substitute(x))
  if ((!is.vector(x))&(!is.matrix(x))&(!is.data.frame(x))) {
    stop('Input x must be a vector, matrix or data frame') }
  DATA <- as.data.frame(x)

  #Check input logical inputs
  if (!is.logical(relative.freq))               stop('Input relative.freq should be a logical value')
  if (length(relative.freq) != 1)               stop('Input relative.freq should be a single logical value')
  if (!is.logical(smooth.line))                 stop('Input smooth.line should be a logical value')
  if (length(smooth.line) != 1)                 stop('Input smooth.line should be a single logical value')
  if (!is.logical(smooth.conf))                 stop('Input smooth.conf should be a logical value')
  if (length(smooth.conf) != 1)                 stop('Input smooth.conf should be a single logical value')
  if (!is.logical(separate.plots))              stop('Input separate.plots should be a logical value')
  if (length(separate.plots) != 1)              stop('Input separate.plots should be a single logical value')
  if (!is.logical(data.name))                   stop('Input data.name should be a logical value')
  if (length(data.name) != 1)                   stop('Input data.name should be a single logical value')

  #Check other graphical inputs
  if (!is.numeric(conf.level))                  stop('Input conf.level must be numeric')
  if (length(conf.level) != 1)                  stop('Input conf.level must be a single numeric value')
  if (conf.level <= 0)                          stop('Input conf.level must be above zero')
  if (conf.level >= 1)                          stop('Input conf.level must be below one')
  if (!is.numeric(point.size))                  stop('Input point.size must be numeric')
  if (length(point.size) != 1)                  stop('Input point.size must be a single numeric value')
  if (point.size <= 0)                          stop('Input point.size must be positive')
  if (!is.numeric(point.alpha))                 stop('Input point.alpha must be numeric')
  if (length(point.alpha) != 1)                 stop('Input point.alpha must be a single numeric value')
  if (point.alpha <= 0)                         stop('Input point.alpha must be positive')
  if (point.alpha > 1)                          stop('Input point.alpha cannot be greater than one')

  #Check installed packages and load them
  GGPLOT2   <- requireNamespace('ggplot2', quietly = TRUE)
  SCALES    <- requireNamespace('scales',  quietly = TRUE)
  if (GGPLOT2) {  } else { stop('Error: Zipf plot requires the ggplot2 package') }
  if (SCALES)  {  }  else { stop('Error: Zipf plot requires the scales package')  }

  #Set theme and colours
  THEME <- ggplot2::theme(plot.title    = ggplot2::element_text(hjust = 0.5, size = 14, face = 'bold'),
                          plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8, face = 'bold', colour = 'darkred'),
                          axis.title.x  = ggplot2::element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                          axis.title.y  = ggplot2::element_text(margin = ggplot2::margin(t = 0,  r = 5, b = 0, l = 0)))

  ###############################################################################################
  #########################################  ZIPF PLOT  #########################################
  ###############################################################################################

  #Compute rank and frequency statistics
  n <- nrow(DATA)
  m <- ncol(DATA)
  RANKS <- vector(mode = 'list', length = m)
  for (k in 1:m) {
    TABLE <- table(DATA[, k])
    DF    <- data.frame(VAR = colnames(DATA)[k],
                        RR  = rank(-c(TABLE), ties.method = 'min'),
                        FF  = c(TABLE))
    if (relative.freq) { DF$FF <- DF$FF/n }
    RANKS[[k]] <- DF[order(DF$RR), ] }
  PLOTDATA <- do.call('rbind', RANKS)
  rownames(PLOTDATA) <- 1:nrow(PLOTDATA)

  #Set subtitle
  nn <- format(n, big.mark = ',', scientific = FALSE)
  if (data.name) {
    if (m == 1) {
      SUBTITLE <- paste0('Data vector ', DATA.NAME, ' contains ', nn, ' values') } else {
      SUBTITLE <- paste0('Data-frame ', DATA.NAME, ' contains ', m, ' variables each with ', nn, ' values') } }
  if (!data.name) {
    if (m == 1) {
      SUBTITLE <- paste0('Data vector contains ', nn, ' values') } else {
        SUBTITLE <- paste0('Data-frame contains ', m, ' variables each with ', nn, ' values') } }
  if ((smooth.line)&(smooth.conf)) {
    SUBTITLE <- paste0(SUBTITLE, '\n(Bands around the smoothing line show ', round(100*conf.level, 2), '% CI)') }

  #Generate plot
  ZIPFPLOT  <- ggplot2::ggplot(ggplot2::aes(x = !!quote(RR), y = !!quote(FF), colour = !!quote(VAR), fill = !!quote(VAR)), data = PLOTDATA) +
               ggplot2::geom_point(size = point.size) +
               { if (smooth.line) ggplot2::geom_smooth(formula = y ~ x, method = 'loess',
                                                       se = smooth.conf, level = conf.level) } +
               ggplot2::scale_x_log10() +
               ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                             labels = scales::trans_format("log10", function(.x) scales::label_math(10^.x)(.x))) +
               ggplot2::expand_limits(y = ifelse(relative.freq, 1, n)) +
               { if (separate.plots) ggplot2::facet_wrap(~ VAR) } +
               THEME + ggplot2::theme(legend.title = ggplot2::element_blank(),
                             legend.position = ifelse(m == 1, 'none', 'bottom'),
                             legend.spacing.x = grid::unit(0.5, 'cm')) +
               ggplot2::ggtitle('Zipf Plot') +
               ggplot2::labs(subtitle = SUBTITLE) +
               ggplot2::xlab('Rank') +
               ggplot2::ylab(ifelse(relative.freq, 'Relative Frequency', 'Frequency'))

  #Print the plot
  ZIPFPLOT }

