#' Generate tail plots for a data vector
#'
#' \code{tailplot} generates the tail plot and Hill plot for the input data
#'
#' The tail plot for a dataset shows the rate of decay of the tails in the data.  It is used to
#' diagnose whether certain moments of the underlying distribution exist (e.g., the variance),
#' which is used in turn to determine whether certain statistical laws apply to the distribution
#' (e.g., the central limit theorem).  The Hill plot shows the adjusted Hill estimator for the tail
#' index of the distribution.  The DSM (De-Sousa-Michailidis) plot shows the adjusted DSM estimator
#' for the tail index of the distribution.  These latter estimators are very similar; more details
#' are found in the reference below.  Our adjusted versions of these estimators are computed using
#' data for the left and right deviations from the extreme values of the dataset; this is done to
#' ensure that the estimators are shift-invariant (the standard Hill and DSM estimators are not).
#'
#' The present function produces tail plots and Hill/DSM plots for the input data vector to show
#' the rate of decay in the tails.  By default, both the tails are shown, but the user can show
#' the plots only for one tail if preferred.  The user can turn any of the plots on or off using
#' the inputs to the function.
#'
#' By default, the tail plot includes lines showing cubic decay in the tails --- if the tails of
#' the distribution decay faster than cubic decay then the variance of the distribution exists.
#' The user can change the order of the lines in the plot to show other decay orders; this can be
#' used to diagnose the existence of moments of other orders.
#'
#' De Sousa, B. and Michailidis, G. (2004) A Diagnostic Plot for Estimating the Tail Index of a
#' Distribution.  Journal of Computational and Graphical Statistics 13(4), pp. 1-22.
#'
#' @usage tailplot(x, tail.prop = 0.05, left = TRUE, right = TRUE,
#'        show.lines = TRUE, lines = 16, line.order = 3,
#'        tail.plot = TRUE, hill.plot = FALSE, dsm.plot = FALSE,
#'        point.size = 3, point.alpha = 0.4,
#'        point.color = NULL, point.colour = point.color,
#'        line.color  = NULL, line.colour  = line.color,
#'        ytop.mult.left = 3, ytop.mult.right = 3)
#' @param x Data vector (numeric)
#' @param tail.prop The proportion of values to use in each tail
#' @param left Logical; if \code{TRUE} the tail plot includes a plot for the left tail
#' @param right Logical; if \code{TRUE} the tail plot includes a plot for the upper tail
#' @param show.lines Logical; if \code{TRUE} the plots include lines for fixed logarithmic decay
#' @param lines Number of lines to include in the plot (included if \code{show.lines} is \code{TRUE})
#' @param line.order Order of the lines in the plot (included if \code{show.lines} is \code{TRUE})
#' @param tail.plot Logical; if \code{TRUE} the output includes the tail plot
#' @param hill.plot Logical; if \code{TRUE} the output includes the Hill plot
#' @param dsm.plot Logical; if \code{TRUE} the output includes the DSM plot
#' @param point.size Size of the points in the plots
#' @param point.alpha Alpha-transparency of the points in the plots
#' @param point.color,point.colour Colour of the points in the plots (default is blue)
#' @param line.color,line.colour Colour of the lines in the plots (default is darkred)
#' @param ytop.mult.left Multiplier used to determine the height of axis for left Hill/DSM plots
#' @param ytop.mult.right Multiplier used to determine the height of axis for right Hill/DSM plots
#' @return Tail plots for the input data (and Hill plots or DSM plots if requested)
#' @examples
#' try(tailplot(rnorm(500)))

tailplot <- function(x, tail.prop = 0.05, left = TRUE, right = TRUE,
                     show.lines = TRUE, lines = 16, line.order = 3,
                     tail.plot = TRUE, hill.plot = FALSE, dsm.plot = FALSE,
                     point.size = 3, point.alpha = 0.4,
                     point.color = NULL, point.colour = point.color,
                     line.color  = NULL, line.colour  = line.color,
                     ytop.mult.left = 3, ytop.mult.right = 3) {

  #Check inputs x and tail.prop
  if (!is.vector(x))                            stop('Error: Input x must be a numeric vector')
  if (!is.numeric(x))                           stop('Error: Input x must be a numeric vector')
  if (length(x) == 0)                           stop('Error: Input x is empty --- no tail plot to print')
  if (!is.numeric(tail.prop))                   stop('Error: Input tail.prop must be numeric')
  if (length(tail.prop) != 1)                   stop('Error: Input tail.prop must be a single numeric value')
  if (tail.prop <= 0)                           stop('Error: Input tail.prop must be greater than zero')
  if (tail.prop > 1)                            stop('Error: Input tail.prop cannot be greater than one')

  #Check inputs left and right
  if (!is.logical(left))                        stop('Error: Input left must be a logical value')
  if (!is.logical(right))                       stop('Error: Input right must be a logical value')
  if (length(left) != 1)                        stop('Error: Input left must be a single logical value')
  if (length(right) != 1)                       stop('Error: Input right must be a single logical value')
  if ((!left)&(!right))                         stop('Error: Input left or right must be TRUE')

  #Check inputs for types of plots to print
  if (!is.logical(tail.plot))                   stop('Error: Input tail.plot must be a logical value')
  if (length(tail.plot) != 1)                   stop('Error: Input tail.plot must be a single logical value')
  if (!is.logical(hill.plot))                   stop('Error: Input hill.plot must be a logical value')
  if (length(hill.plot) != 1)                   stop('Error: Input hill.plot must be a single logical value')
  if (!is.logical(dsm.plot))                    stop('Error: Input dsm.plot must be a logical value')
  if (length(dsm.plot) != 1)                    stop('Error: Input dsm.plot must be a single logical value')
  if ((!tail.plot)&&(!hill.plot)&&(!dsm.plot))  stop('No plot returned')

  #Check inputs for lines
  if (!is.logical(show.lines))                  stop('Error: Input show.lines must be a logical value')
  if (length(show.lines) != 1)                  stop('Error: Input show.lines must be a single logical value')
  if (!is.numeric(lines))                       stop('Error: Input lines must be numeric')
  if (length(lines) != 1)                       stop('Error: Input lines must be a single numeric value')
  if (as.integer(lines) != lines)               stop('Error: Input lines must be an integer')
  if (lines <= 1)                               stop('Error: Input lines must be larger than one')
  if (!is.numeric(line.order))                  stop('Error: Input line.order must be numeric')
  if (length(line.order) != 1)                  stop('Error: Input line.order must be a single numeric value')
  if (line.order < 1)                           stop('Error: Input line.order must be at least one')

  #Check graphical inputs
  if (!is.numeric(point.size))                  stop('Error: Input point.size must be numeric')
  if (length(point.size) != 1)                  stop('Error: Input point.size must be a single numeric value')
  if (point.size <= 0)                          stop('Error: Input point.size must be positive')
  if (!is.numeric(point.alpha))                 stop('Error: Input point.alpha must be numeric')
  if (length(point.alpha) != 1)                 stop('Error: Input point.alpha must be a single numeric value')
  if (point.alpha <= 0)                         stop('Error: Input point.alpha must be positive')
  if (point.alpha > 1)                          stop('Error: Input point.alpha cannot be greater than one')
  if ((!missing(point.color))&(!missing(point.colour))) {
    stop('Error: Specify point.color or point.colour but not both') }
  if (!is.null(point.colour)) {
    if (!is.character(point.colour))            stop('Error: point.colour must be in \'colours()\'')
    if (length(point.colour) != 1)              stop('Error: point.colour must be in \'colours()\'')
    if (!(point.colour %in% colours()))         stop('Error: point.colour must be in \'colours()\'') }
  if ((!missing(line.color))&(!missing(line.colour))) {
    stop('Error: Specify line.color or line.colour but not both') }
  if (!is.null(line.colour)) {
    if (!is.character(line.colour))             stop('Error: line.colour must be in \'colours()\'')
    if (length(line.colour) != 1)               stop('Error: line.colour must be in \'colours()\'')
    if (!(line.colour %in% colours()))          stop('Error: line.colour must be in \'colours()\'') }
  if(!is.numeric(ytop.mult.left))               stop('Error: ytop.mult.left should be numeric')
  if(length(ytop.mult.left) != 1)               stop('Error: ytop.mult.left should be a single numeric value')
  if(ytop.mult.left <= 0)                       stop('Error: ytop.mult.left should be positive')
  if(!is.numeric(ytop.mult.right))              stop('Error: ytop.mult.right should be numeric')
  if(length(ytop.mult.right) != 1)              stop('Error: ytop.mult.right should be a single numeric value')
  if(ytop.mult.right <= 0)                      stop('Error: ytop.mult.right should be positive')

  #Set sizes for tail plots
  n <- length(x)
  m <- ceiling(n*tail.prop)

  #Check installed packages and load them
  GGPLOT2   <- requireNamespace('ggplot2', quietly = TRUE)
  SCALES    <- requireNamespace('scales',  quietly = TRUE)
  if (GGPLOT2) {  } else { stop('Error: Tail plot requires the ggplot2 package') }
  if (SCALES)  {  }  else { stop('Error: Tail plot requires the scales package')  }
  if ((left)&(right)) {
    GRIDEXTRA <- requireNamespace('gridExtra',  quietly = TRUE)
    if (GRIDEXTRA)  {  } else {
                      stop('Error: Tail plot requires the gridExtra package') } }

  #Set theme and colours
  THEME <- ggplot2::theme(plot.title    = ggplot2::element_text(hjust = 0.5, size = 14, face = 'bold'),
                          plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8, face = 'bold', colour = 'darkred'),
                          axis.title.x  = ggplot2::element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                          axis.title.y  = ggplot2::element_text(margin = ggplot2::margin(t = 0,  r = 5, b = 0, l = 0)))
  if (!is.null(point.colour)) { POINT.COLOUR <- point.colour } else { POINT.COLOUR <- 'blue' }
  if (!is.null(line.colour))  { LINE.COLOUR  <- line.colour  } else { LINE.COLOUR  <- 'darkred' }

  #########################################  TAIL PLOT  #########################################

  #Generate the plot data and ranges for each tail
  XORD <- sort(x, decreasing = FALSE)
  DEV.LEFT  <- max(x) - XORD
  DEV.RIGHT <- XORD - min(x)
  PLOTDATA <- data.frame(deviation.left  = DEV.LEFT,
                         deviation.right = DEV.RIGHT,
                         tail.prob.left  = (2*(1:n)-1)/(2*n),
                         tail.prob.right = (2*(n-1:n)+1)/(2*n))
  LEFT.RANGE   <- 1:m
  RIGHT.RANGE  <- (n-m+1):n

  #Set subtitle for lines
  if (show.lines) {
    ORDER.LABELS <- c('linear', 'quadratic', 'cubic', 'quartic', 'quintic',
                      'sextic', 'septic', 'octic', 'nontic', 'decic')
    ORDER.TYPE   <- c('probability', 'mean', 'variance', 'skewness', 'kurtosis',
                      'hyperskewness', 'hyperkurtosis', 'second hyperskewness',
                      'second hyperkurtosis', 'third hyperskewness')
    if (line.order %in% 1:10) {
      SUBTITLE <- paste0('(Diagonal lines show ', ORDER.LABELS[line.order],
                         ' tail decay)', '\n',
                         '(Faster decay than this gives finite ',
                         ORDER.TYPE[line.order], ')')
    } else {
      SUBTITLE <- paste0('(Diagonal lines show tail decay of order = ',
                         round(line.order, 4), ')', '\n',
                         '(Faster decay than this gives finite moment of order ',
                         round(line.order-1, 4), ')') } }


  ###### Create the left plot #######
  if (left) {

    #Set cubic lines
    if (show.lines) {

      #Compute intercepts for cubic lines
      #We add two cubic lines on each side (beyond cubic.lines) for padding
      MIN.X <- log10(min(PLOTDATA$deviation.left[LEFT.RANGE]))
      MAX.X <- log10(max(PLOTDATA$deviation.left[LEFT.RANGE]))
      MIN.Y <- log10(min(PLOTDATA$tail.prob.left[LEFT.RANGE]))
      MAX.Y <- log10(max(PLOTDATA$tail.prob.left[LEFT.RANGE]))
      DIST  <- (MAX.Y - MIN.Y) + line.order*(MAX.X - MIN.X)
      INTS  <- MIN.Y + (DIST/(lines-1))*((-2):(lines+1)) + line.order*MIN.X

      #Generate plot with cubic lines
      TAILPLOT.LEFT  <- ggplot2::ggplot(ggplot2::aes(x = !!quote(deviation.left), y = !!quote(tail.prob.left)),
                                        data = PLOTDATA[LEFT.RANGE, ]) +
        ggplot2::geom_abline(intercept = INTS, slope = -line.order,
                             color = LINE.COLOUR, linetype = 'dashed') +
        ggplot2::geom_point(size = point.size, alpha = point.alpha, colour = POINT.COLOUR) +
        ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                               labels = scales::trans_format("log10", function(.x) scales::label_math(10^.x)(.x))) +
        ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                               labels = scales::trans_format("log10", function(.x) scales::label_math(10^.x)(.x))) +
        THEME +
        ggplot2::ggtitle('Tail Plot (Left)') +
        ggplot2::labs(subtitle = SUBTITLE) +
        ggplot2::xlab('Left-Deviation from Maximum Value') +
        ggplot2::ylab('Empirical Tail Probability')

    } else {

      #Generate plot without cubic lines
      TAILPLOT.LEFT  <- ggplot2::ggplot(ggplot2::aes(x = !!quote(deviation.left), y = !!quote(tail.prob.left)),
                                        data = PLOTDATA[LEFT.RANGE, ]) +
        ggplot2::geom_point(size = point.size, alpha = point.alpha, colour = POINT.COLOUR) +
        ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                               labels = scales::trans_format("log10", function(.x) scales::label_math(10^.x)(.x))) +
        ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                               labels = scales::trans_format("log10", function(.x) scales::label_math(10^.x)(.x))) +
        THEME +
        ggplot2::ggtitle('Tail Plot (Left)') +
        ggplot2::labs(subtitle = SUBTITLE) +
        ggplot2::xlab('Left-Deviation from Maximum Value') +
        ggplot2::ylab('Empirical Tail Probability') } }

  ###### Create the right plot #####
  if (right) {

    #Set cubic lines
    if (show.lines) {

      #Compute intercepts for cubic lines
      #We add two cubic lines on each side (beyond cubic.lines) for padding
      MIN.X <- log10(min(PLOTDATA$deviation.right[RIGHT.RANGE]))
      MAX.X <- log10(max(PLOTDATA$deviation.right[RIGHT.RANGE]))
      MIN.Y <- log10(min(PLOTDATA$tail.prob.right[RIGHT.RANGE]))
      MAX.Y <- log10(max(PLOTDATA$tail.prob.right[RIGHT.RANGE]))
      DIST  <- (MAX.Y - MIN.Y) + line.order*(MAX.X - MIN.X)
      INTS  <- MIN.Y + (DIST/(lines-1))*((-2):(lines+1)) + line.order*MIN.X

      #Generate plot with cubic lines
      TAILPLOT.RIGHT <- ggplot2::ggplot(ggplot2::aes(x = !!quote(deviation.right), y = !!quote(tail.prob.right)),
                                       data = PLOTDATA[RIGHT.RANGE, ]) +
                        ggplot2::geom_abline(intercept = INTS, slope = -line.order,
                             color = LINE.COLOUR, linetype = 'dashed') +
                        ggplot2::geom_point(size = point.size, alpha = point.alpha, colour = POINT.COLOUR) +
                        ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                                    labels = scales::trans_format("log10", function(.x) scales::label_math(10^.x)(.x))) +
                        ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                                    labels = scales::trans_format("log10", function(.x)scales::label_math(10^.x)(.x))) +
                        THEME +
                        ggplot2::ggtitle('Tail Plot (Right)') +
                        ggplot2::labs(subtitle = SUBTITLE) +
                        ggplot2::xlab('Right-Deviation from Minimum Value') +
                        if (!left) {
                          ggplot2::ylab('Empirical Tail Probability') } else {
                          ggplot2::ylab(NULL) }

    } else {

      #Generate plot without cubic lines
      TAILPLOT.RIGHT <- ggplot2::ggplot(ggplot2::aes(x = !!quote(deviation.right), y = !!quote(tail.prob.right)),
                                        data = PLOTDATA[RIGHT.RANGE, ]) +
                        ggplot2::geom_point(size = point.size, alpha = point.alpha, colour = POINT.COLOUR) +
                        ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                                    labels = scales::trans_format("log10", function(.x) scales::math_format(10^.x)(.x))) +
                        ggplot2::scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                                    labels = scales::trans_format("log10", function(.x) scales::math_format(10^.x)(.x))) +
                        THEME +
                        ggplot2::ggtitle('Tail Plot (Right)') +
                        ggplot2::labs(subtitle = SUBTITLE) +
                        ggplot2::xlab('Right-Deviation from Minimum Value') +
                        if (!left) {
                          ggplot2::ylab('Empirical Tail Probability') } else {
                          ggplot2::ylab(NULL) } } }

  ##############################  HILL/DE-SOUSA-MICHAILIDIS PLOTS  ##############################

  #Create plotting data
  PLOTDATA2 <- data.frame(Index = 1:m, H.LEFT = 0, H.RIGHT = 0, S.LEFT = 0, S.RIGHT = 0)
  for (k in 1:m) {
    PLOTDATA2$H.LEFT[k]  <- ((n-2)/n)*(k/sum(base::log(DEV.LEFT[1:k]/DEV.LEFT[k+1])))
    PLOTDATA2$H.RIGHT[k] <- ((n-2)/n)*(k/sum(base::log(DEV.RIGHT[n+1-(1:k)]/DEV.RIGHT[n-k])))
    PLOTDATA2$S.LEFT[k]  <- ((n-2)/n)*(k/sum((1:k)*base::log(DEV.LEFT[1:k]/DEV.LEFT[2:(k+1)])))
    PLOTDATA2$S.RIGHT[k] <- ((n-2)/n)*(k/sum((1:k)*base::log(DEV.RIGHT[n:(n-k+1)]/DEV.RIGHT[(n-1):(n-k)]))) }

  #Set maximum vertical axis values
  MM.H.LEFT   <- mean(PLOTDATA2$H.LEFT)
  MM.S.LEFT   <- mean(PLOTDATA2$S.LEFT)
  MAX.H.LEFT  <- min(max(PLOTDATA2$H.LEFT), ytop.mult.left*MM.H.LEFT)
  MAX.S.LEFT  <- min(max(PLOTDATA2$S.LEFT), ytop.mult.left*MM.S.LEFT)
  MM.H.RIGHT  <- mean(PLOTDATA2$H.RIGHT)
  MM.S.RIGHT  <- mean(PLOTDATA2$S.RIGHT)
  MAX.H.RIGHT <- min(max(PLOTDATA2$H.RIGHT), ytop.mult.right*MM.H.RIGHT)
  MAX.S.RIGHT <- min(max(PLOTDATA2$S.RIGHT), ytop.mult.right*MM.S.RIGHT)

  #Set subtitles
  SUBTITLE2 <- '(Adjusted Hill estimator for the tail index)'
  SUBTITLE3 <- '(Adjusted DSM estimator for the tail index)'

  #Create the left plots
  if (left) {

    #Generate plot
    HILLPLOT.LEFT  <- ggplot2::ggplot(ggplot2::aes(x = !!quote(Index), y = !!quote(H.LEFT)), data = PLOTDATA2) +
      ggplot2::geom_point(size = point.size, alpha = point.alpha, colour = POINT.COLOUR) +
      ggplot2::scale_y_continuous(limits = c(0, MAX.H.LEFT)) +
      THEME +
      ggplot2::ggtitle('Hill Plot (Left)') +
      ggplot2::labs(subtitle = SUBTITLE2) +
      ggplot2::xlab('Points in Tail') +
      ggplot2::ylab('Adjusted Hill Estimator')
    DSMPLOT.LEFT   <- ggplot2::ggplot(ggplot2::aes(x = !!quote(Index), y = !!quote(S.LEFT)), data = PLOTDATA2) +
      ggplot2::geom_point(size = point.size, alpha = point.alpha, colour = POINT.COLOUR) +
      ggplot2::scale_y_continuous(limits = c(0, MAX.S.LEFT)) +
      THEME +
      ggplot2::ggtitle('DSM Plot (Left)') +
      ggplot2::labs(subtitle = SUBTITLE3) +
      ggplot2::xlab('Points in Tail') +
      ggplot2::ylab('Adjusted DSM Estimator') }

  #Create the right plots
  if (right) {

    #Generate plots
    HILLPLOT.RIGHT <- ggplot2::ggplot(ggplot2::aes(x = !!quote(Index), y = !!quote(H.RIGHT)), data = PLOTDATA2) +
      ggplot2::geom_point(size = point.size, alpha = point.alpha, colour = POINT.COLOUR) +
      ggplot2::scale_y_continuous(limits = c(0, MAX.H.RIGHT)) +
      THEME +
      ggplot2::ggtitle('Hill Plot (Right)') +
      ggplot2::labs(subtitle = SUBTITLE2) +
      ggplot2::xlab('Points in Tail') +
      if (!left) {
        ggplot2::ylab('Adjusted Hill Estimator') } else {
        ggplot2::ylab(NULL) }
    DSMPLOT.RIGHT  <- ggplot2::ggplot(ggplot2::aes(x = !!quote(Index), y = !!quote(S.RIGHT)), data = PLOTDATA2) +
      ggplot2::geom_point(size = point.size, alpha = point.alpha, colour = POINT.COLOUR) +
      ggplot2::scale_y_continuous(limits = c(0, MAX.S.RIGHT)) +
      THEME +
      ggplot2::ggtitle('DSM Plot (Right)') +
      ggplot2::labs(subtitle = SUBTITLE3) +
      ggplot2::xlab('Points in Tail') +
      if (!left) {
        ggplot2::ylab('Adjusted DSM Estimator') } else {
        ggplot2::ylab(NULL) } }

  #########################################  FULL PLOT  #########################################

  #Set spacing object
  SPACE <- ggplot2::theme(plot.margin = ggplot2::margin(t = 5.5, r = 5.5, b = 16.5, l = 5.5))

  #Create overall plot
  if ((left)&(right)) {

    #Set list of plots
    PLOTS <- list(TAILPLOT.LEFT + SPACE, TAILPLOT.RIGHT + SPACE,
                  HILLPLOT.LEFT + SPACE, HILLPLOT.RIGHT + SPACE,
                  DSMPLOT.LEFT  + SPACE, DSMPLOT.RIGHT  + SPACE)

    #Set combined plot
    if ((tail.plot)&&(hill.plot)&&(dsm.plot)) {
      PLOT <- gridExtra::grid.arrange(PLOTS[[1]], PLOTS[[2]], PLOTS[[3]],
                                      PLOTS[[4]], PLOTS[[5]], PLOTS[[6]],
                                      nrow = 3, ncol = 2) }
    if ((tail.plot)&&(hill.plot)&&(!dsm.plot)) {
      PLOT <- gridExtra::grid.arrange(PLOTS[[1]], PLOTS[[2]], PLOTS[[3]], PLOTS[[4]],
                                      nrow = 2, ncol = 2) }
    if ((tail.plot)&&(!hill.plot)&&(dsm.plot)) {
      PLOT <- gridExtra::grid.arrange(PLOTS[[1]], PLOTS[[2]], PLOTS[[5]], PLOTS[[6]],
                                      nrow = 2, ncol = 2) }
    if ((tail.plot)&&(!hill.plot)&&(!dsm.plot)) {
      PLOT <- gridExtra::grid.arrange(PLOTS[[1]], PLOTS[[2]],
                                      nrow = 1, ncol = 2) }

    if ((!tail.plot)&&(hill.plot)&&(dsm.plot)) {
      PLOT <- gridExtra::grid.arrange(PLOTS[[3]], PLOTS[[4]], PLOTS[[5]], PLOTS[[6]],
                                      nrow = 2, ncol = 2) }
    if ((!tail.plot)&&(hill.plot)&&(!dsm.plot)) {
      PLOT <- gridExtra::grid.arrange(PLOTS[[3]], PLOTS[[4]],
                                      nrow = 1, ncol = 2) }
    if ((!tail.plot)&&(!hill.plot)&&(dsm.plot)) {
      PLOT <- gridExtra::grid.arrange(PLOTS[[5]], PLOTS[[6]],
                                      nrow = 1, ncol = 2) } }

  if ((left)&(!right)) {

    #Set list of plots
    PLOTS <- list(TAILPLOT.LEFT + SPACE,
                  HILLPLOT.LEFT + SPACE,
                  DSMPLOT.LEFT  + SPACE)

    #Set combined plot
    if ((tail.plot)&&(hill.plot)&&(dsm.plot)) {
      PLOT <- gridExtra::grid.arrange(PLOTS[[1]], PLOTS[[2]], PLOTS[[3]],
                                      nrow = 3, ncol = 1) }
    if ((tail.plot)&&(hill.plot)&&(!dsm.plot)) {
      PLOT <- gridExtra::grid.arrange(PLOTS[[1]], PLOTS[[2]],
                                      nrow = 2, ncol = 1) }
    if ((tail.plot)&&(!hill.plot)&&(dsm.plot)) {
      PLOT <- gridExtra::grid.arrange(PLOTS[[1]], PLOTS[[3]],
                                      nrow = 2, ncol = 1) }
    if ((tail.plot)&&(!hill.plot)&&(!dsm.plot)) {
      PLOT <- gridExtra::grid.arrange(PLOTS[[1]],
                                      nrow = 1, ncol = 1) }
    if ((!tail.plot)&&(hill.plot)&&(dsm.plot)) {
      PLOT <- gridExtra::grid.arrange(PLOTS[[2]], PLOTS[[3]],
                                      nrow = 2, ncol = 1) }
    if ((!tail.plot)&&(hill.plot)&&(!dsm.plot)) {
      PLOT <- gridExtra::grid.arrange(PLOTS[[2]],
                                      nrow = 1, ncol = 1) }
    if ((!tail.plot)&&(!hill.plot)&&(dsm.plot)) {
      PLOT <- gridExtra::grid.arrange(PLOTS[[3]],
                                      nrow = 1, ncol = 1) } }

  if ((!left)&(right)) {

    #Set list of plots
    PLOTS <- list(TAILPLOT.RIGHT + SPACE,
                  HILLPLOT.RIGHT + SPACE,
                  DSMPLOT.RIGHT  + SPACE)

    #Set combined plot
    if ((tail.plot)&&(hill.plot)&&(dsm.plot)) {
      PLOT <- gridExtra::grid.arrange(PLOTS[[1]], PLOTS[[2]], PLOTS[[3]],
                                      nrow = 3, ncol = 1) }
    if ((tail.plot)&&(hill.plot)&&(!dsm.plot)) {
      PLOT <- gridExtra::grid.arrange(PLOTS[[1]], PLOTS[[3]],
                                      nrow = 2, ncol = 1) }
    if ((tail.plot)&&(!hill.plot)&&(dsm.plot)) {
      PLOT <- gridExtra::grid.arrange(PLOTS[[1]], PLOTS[[3]],
                                      nrow = 2, ncol = 1) }
    if ((tail.plot)&&(!hill.plot)&&(!dsm.plot)) {
      PLOT <- gridExtra::grid.arrange(PLOTS[[1]],
                                      nrow = 1, ncol = 1) }
    if ((!tail.plot)&&(hill.plot)&&(dsm.plot)) {
      PLOT <- gridExtra::grid.arrange(PLOTS[[2]], PLOTS[[3]],
                                      nrow = 2, ncol = 1) }
    if ((!tail.plot)&&(hill.plot)&&(!dsm.plot)) {
      PLOT <- gridExtra::grid.arrange(PLOTS[[3]],
                                      nrow = 1, ncol = 1) }
    if ((!tail.plot)&&(!hill.plot)&&(dsm.plot)) {
      PLOT <- gridExtra::grid.arrange(PLOTS[[3]],
                                      nrow = 1, ncol = 1) } }

  #Print the plot
  PLOT }

