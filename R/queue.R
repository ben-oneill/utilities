#' Generate queuing information from arrival and use times
#'
#' \code{queue} returns queuing information for users and service facilities.
#'
#' This function computes takes inputs giving the arrival times and (intended) use times for a set of users at
#' an amenity, plus the number of service facilities at the amenity.  The function computes full information on
#' the use of the facilities by the users, including their waiting time, actual use time, leaving time, and the
#' facility that was used by each user.
#'
#' In addition to the required inputs, the function also accepts inputs for a maximum-waiting time for each user;
#' if the user waits up to this time then the user will leave without service.  The user can also impose closure
#' times on new arrivals, new services, or termination of services.
#'
#' **Note:** Service facilities are assumed to be allocated to users on a "first-come, first-served" basis; in the
#' event that more than one service facility is available for a user then the user is allocated to facilities
#' first-to-last based on the facility number (i.e., the allocation favours the earlier facilities and it is
#' not exchangeable with respect to the facility number).
#'
#' @param n Number of service facilities at the amenity (positive integer)
#' @param arrive Vector of arrival-times for the users (non-negative numeric values)
#' @param use.full Vector of (intended) use-times for the users (non-negative numeric values)
#' @param wait.max Vector of maximum-waiting-times for the users (non-negative numeric values)
#' @param revive Revival-time for service facilities
#' @param close.arrive Closure-time for new arrivals (no new arrivals allowed)
#' @param close.service Closure-time for new services (no new services allowed)
#' @param close.full Closure-time for all services (all existing services are terminated)
#' @return If all inputs are correctly specified then the function will return a list of class \code{queue}
#' containing queuing information for the users and service facilities
#'
#' @examples
#' q <- queue(2, 4:6, 7:9)
#' summary(q)
#' plot(q)
#' plot(summary(q))

queue <- function(n, arrive, use.full, wait.max = NULL, revive = 0,
                  close.arrive = Inf, close.service = Inf, close.full = Inf) {

  #Check main inputs
  if (!is.numeric(arrive))           stop('Error: Input arrive should be numeric')
  if (!is.numeric(use.full))         stop('Error: Input use.full should be numeric')
  if (!is.numeric(n))                stop('Error: Input n should be numeric')
  if (length(arrive) != length(use.full))   stop('Error: Inputs arrive and use.full must have same length')
  if (length(n) != 1)                stop('Error: Input n should be a single value')
  if (min(arrive) < 0)               stop('Error: Input arrive should contain non-negative values')
  if (min(use.full) < 0)             stop('Error: Input use.full should contain non-negative values')
  nn <- as.integer(n)
  if (n != nn)                       stop('Error: Input n should be a positive integer')
  if (min(n) <= 0)                   stop('Error: Input n should be a positive integer')
  if (max(n) == Inf)                 stop('Error: Input n should be a positive integer (not infinity)')

  #Check input wait.max and create maximum-waiting-time matrix
  K <- length(arrive)
  WW <- .convert.wait.max(wait.max, K)

  #Check inputs for service facility
  if (!is.numeric(revive))           stop('Error: Input revive should be numeric')
  if (!is.numeric(close.arrive))     stop('Error: Input close.arrive should be numeric')
  if (!is.numeric(close.service))    stop('Error: Input close.service should be numeric')
  if (!is.numeric(close.full))       stop('Error: Input close.full should be numeric')
  if (length(revive) != 1)           stop('Error: Input revive should be a single value')
  if (length(close.arrive) != 1)     stop('Error: Input close.arrive should be a single value')
  if (length(close.service) != 1)    stop('Error: Input close.service should be a single value')
  if (length(close.full) != 1)       stop('Error: Input close.full should be a single value')
  if (min(revive) < 0)               stop('Error: Input revive should be a non-negative value')
  if (min(close.arrive) < 0)         stop('Error: Input close.arrive should be a non-negative value')
  if (min(close.service) < 0)        stop('Error: Input close.service should be a non-negative value')
  if (min(close.full) < 0)           stop('Error: Input close.full should be a non-negative value')
  close.service <- min(close.service, close.full)
  close.arrive  <- min(close.arrive, close.service, close.full)

  #Set the parameters
  t   <- arrive
  uu  <- use.full
  ww  <- rep(0, K)
  ORD <- order(t)
  t   <- t[ORD]
  uu  <- uu[ORD]
  r   <- revive
  T1  <- close.arrive
  T2  <- close.service
  T3  <- close.full

  #Set output vectors
  w  <- numeric(K)
  u  <- numeric(K)
  F  <- numeric(K)
  for (k in 1:K) { F[k] <- NA }
  e  <- numeric(K)
  L  <- t

  #Set delay matrix
  D  <- matrix(0, nrow = K+1, ncol = n)
  rownames(D) <- sprintf('Arrival[%s]', 0:K)
  colnames(D) <- sprintf('Delay[%s]', 1:n)

  #Set queue-priority matrix
  QP <- matrix(0, nrow = K, ncol = K)
  QP[1,1] <- 1
  rownames(QP) <- sprintf('User[%s]', 1:K)
  colnames(QP) <- sprintf('Priority[%s]', 1:K)

  #Compute user outputs
  TIME <- 0
  for (k in 1:K) {

    #Update delays/time
    D[k+1,] <- pmax(D[k,]-t[k]+TIME, 0)
    TIME    <- t[k]

    #Update user statistics
    if (t[k] < T1) {
      w.hat  <- min(D[k+1, ])
      if (k > 1) {
        end.queue.times <- t[1:(k-1)] + w[1:(k-1)]
        for (kk in 1:k) {
          QP[k, kk] <- 1 + sum(end.queue.times > t[k] + WW[k,kk]) } }
      kk.hat <- max(which(QP[k, ] >= 1:K))
      ww[k] <- WW[k, kk.hat]
      w.ss   <- min(ww[k], T2-TIME)
      w[k]   <- min(w.hat, w.ss)
      SERVED <- (w.hat < w.ss)
      if (SERVED) {
        F[k] <- which.min(D[k+1, ])
        u[k] <- min(w[k] + uu[k], T3-TIME) - min(w[k], T3-TIME)
        L[k] <- t[k] + w[k] + u[k]
        e[k] <- u[k]/uu[k]
        D[k+1, F[k]] <- D[k+1, F[k]] + u[k] + r } } }

  #Generate user information
  USER <- data.frame(arrive = t, wait = w, use = u, leave = L,
                     unserved = uu-u, F = F,
                     use.full = uu, use.prop = e, wait.max = ww)
  rownames(USER) <- sprintf('user[%s]', 1:K)

  #Generate user-facility matrix
  MATRIX <- matrix(0, nrow = K, ncol = n)
  for (k in 1:K) {
  for (i in 1:n) {
    if (!is.na(F[k])) { if (F[k] == i) { MATRIX[k,i] <- 1 } } } }
  rownames(MATRIX) <- sprintf('user[%s]', 1:K)
  colnames(MATRIX) <- sprintf('F[%s]',    1:n)

  #Generate facility information
  U <- c(u %*% MATRIX)
  R <- c(rep(r, K) %*% MATRIX)
  BEGIN <- rep(0, n)
  END   <- rep(T, n)
  for (i in 1:n) {
    if (sum(MATRIX[,i] == 1) > 0) { END[i] <- max(L[MATRIX[,i] == 1]) } }
  FACILITY <- data.frame(open = BEGIN, end.service = END, use = U, revive = R)
  rownames(FACILITY) <- sprintf('F[%s]', 1:n)

  #Generate output
  OUT  <- list(users = USER, facilities = FACILITY, users.facilities = MATRIX, delay = D, n = n,
               revival = r, close.arrive = T1, close.service = T2, close.full = T3)
  class(OUT) <- c('queue', 'list')

  #Return output
  OUT }


.convert.wait.max <- function(wait.max = NULL, K) {

  #Check input wait.max
  if (!is.null(wait.max)) {

    #If wait.max is a vector then generate numeric matrix WW
    if (is.vector(wait.max)) {
      if (is.numeric(wait.max)) {
        if (length(wait.max) != K)     stop('Error: If input wait.max is a vector it should the same length as the arrival vector')
        WW <- matrix(rep(wait.max, K), nrow = K, ncol = K, byrow = TRUE) } }

    #If wait.max is a function then generate numeric matrix WW
    if ('function' %in% class(wait.max)) {
      WW <- matrix(0, nrow = K, ncol = K)
      for (k in 1:K) {
        WW[, k] <- wait.max(k) } }

    #If wait.max is a matrix then it is numeric matrix WW
    if (is.vector(wait.max)) {
      WW <- wait.max
      if (is.numeric(WW))            stop('Error: If input wait.max is a matrix it should be numeric')
      if (nrow(WW) != K)             stop('Error: If input wait.max is a matrix it should have one row per user')
      if (ncol(WW) != K)             stop('Error: If input wait.max is a matrix it should have one row per user') }

  } else {

    #If wait.max is null then create matrix WW with infinite maximum waiting-times
    WW <- matrix(Inf, nrow = K, ncol = K) }

    #Check that rows of WW are non-increasing
    for (k in 1:K) {
      if (k > 1) {
      for (kk in 2:k) {
        if (WW[kk] > WW[kk-1])       stop(paste0('Error: Row ', k, ' of wait.max matrix is decreasing')) } } }

  WW }

#' @rdname queue
#' @param x,object a \code{queue} object
#' @param ... further arguments passed to or from other methods.
print.queue <- function(x, ...) {

  #Check input
  if (!('queue' %in% class(x)))     stop('Error: Object must have class \'queue\'' )

  #Extract object information
  USER     <- x$users
  FACILITY <- x$facilities
  MATRIX   <- x$users.facilities
  n        <- x$n
  r        <- x$revival
  T1       <- x$close.arrive
  T2       <- x$close.service
  T3       <- x$close.full

  #Print heading
  PLURAL <- ifelse(n > 1, 'facilities', 'facility')
  cat('\n    Queue Information \n \n')
  if (r == 0) {
    cat(paste0('Model of an amenity with ', n, ' service ', PLURAL, ' \n'))
  } else {
    cat(paste0('Model of an amenity with ', n, ' service ', PLURAL, ' with revival-time ',
                 round(r, 4), '\n')) }
  if (T1 < Inf) { cat('Service facilities close to new arrivals at closure-time =', T1, '\n') }
  if (T1 < Inf) { cat('Service facilities close to new services at closure-time =', T2, '\n') }
  if (T1 < Inf) { cat('Service facilities end existing services at closure-time =', T3, '\n') }
  cat('\n')
  cat('Users are allocated to facilities on a \'first-come first-served\' basis \n \n')

  #Print user data
  cat('---------------------------------------------------------------------- \n \n')
  cat('User information \n \n')
  print(USER[, 1:6])
  cat('\n')

  #Print facility data
  cat('---------------------------------------------------------------------- \n \n')
  cat('Facility information \n \n')
  print(FACILITY)
  cat('\n') }

#' @rdname queue
#' @param print,gap,line.width,line.colors,line.colours plotting paramaters
plot.queue <- function(x, print = TRUE, gap = NULL, line.width = 2,
                       line.colors = NULL, line.colours = line.colors, ...) {

  #Check inputs
  if (!inherits(x, 'queue'))     stop('Error: Object must have class \'queue\'' )
  if (!is.logical(print))                stop('Error: Input print should be a logical value')
  if (length(print) != 1)                stop('Error: Input print should be a single logical value')
  if (!is.null(gap)) {
    if(!is.numeric(gap))                 stop('Error: Input gap must be a non-negative integer')
    if(length(gap) != 1)                 stop('Error: Input gap must be a single non-negative integer')
    gap.int <- as.integer(gap)
    if(gap != gap.int)                   stop('Error: Input gap must be a non-negative integer') }
  if (!is.numeric(line.width))           stop('Error: Input line.width must be a positive value')
  if (length(line.width) != 1)           stop('Error: Input line.width must be a single positive value')
  if (min(line.width) <= 0)              stop('Error: Input line.width must be a positive value')
  if (is.null(line.colours)) {
    COLOURS <- c('Waiting' = 'Red', 'Service' = 'Blue', 'Unserved' = 'Black', 'Revival' = 'Turquoise3')
  } else {
    if (!is.character(line.colours))     stop('Error: Input line.colours must be a character vector with four colours')
    if (length(line.colours) != 4)       stop('Error: Input line.colours must be a character vector with four colours')
    if (!(all(tolower(line.colours) %in% colours())))
      stop('Error: Input line.colours must be a character vector with four colours')
    COLOURS <- c('Waiting' = line.colours[1], 'Service' = line.colours[2],
                 'Unserved' = line.colours[3], 'Revival' = line.colours[4]) }

  #Check required packages
  GGPLOT2 <- requireNamespace('ggplot2', quietly = TRUE)
  if (!GGPLOT2) stop('Error: Plotting a \'queue\' object requires the ggplot2 package')
  GRID <- requireNamespace('grid', quietly = TRUE)
  if (!GRID) stop('Error: Plotting a \'queue\' object requires the grid package')
  GRIDEXTRA <- requireNamespace('gridExtra', quietly = TRUE)
  if (!GRIDEXTRA) stop('Error: Plotting a \'queue\' object requires the gridExtra package')

  #Extract object information
  USER     <- x$users
  FACILITY <- x$facilities
  MATRIX   <- x$users.facilities
  K        <- nrow(USER)
  n        <- x$n
  r        <- x$revival
  T1       <- x$close.arrive
  T2       <- x$close.service
  T3       <- x$close.full

  #Create plot data
  if (is.null(gap)) { GAP <- max(3, ceiling(K/10)) } else { GAP <- gap }
  ARRIVE    <- USER$arrive
  SERVICE   <- USER$arrive + USER$wait
  LEAVE     <- USER$arrive + USER$wait + USER$use
  UNSERVE   <- USER$arrive + USER$wait + USER$use.full
  REVIVE    <- USER$arrive + USER$wait + USER$use + r*(USER$use > 0)
  FACILITY  <- USER$F
  PLOTDATA  <- data.frame(user    = 1:K,   facility = FACILITY,
                          arrive  = ARRIVE, service = SERVICE,
                          leave   = LEAVE,  unserve = UNSERVE, revive = REVIVE)
  PLOTDATA.REDUCED <- PLOTDATA[!is.na(USER$F), ]

  #Create subtitle
  PLURAL   <- ifelse(n > 1, 'facilities', 'facility')
  if (r == 0) {
    SUBTITLE <- paste0('(Queuing process with ', n, ' service ', PLURAL, ')')
  } else {
    SUBTITLE <- paste0('(Queuing process with ', n, ' service ', PLURAL,
                    ' with revival-time ', round(r, 4), ')') }

  #Create the plot
  SIZE <- line.width
  PLOT <- ggplot2::ggplot(data = PLOTDATA) +
    ggplot2::geom_segment(ggplot2::aes(x = !!quote(arrive),  xend = !!quote(service),
                                       y = !!quote(user),    yend = !!quote(user), colour = 'Waiting'), size = SIZE) +
    ggplot2::geom_segment(ggplot2::aes(x = !!quote(service), xend = !!quote(leave),
                                       y = !!quote(user),    yend = !!quote(user), colour = 'Service'), size = SIZE) +
    ggplot2::geom_segment(ggplot2::aes(x = !!quote(leave),   xend = !!quote(unserve),
                                       y = !!quote(user),    yend = !!quote(user), colour = 'Unserved'), size = SIZE) +
    ggplot2::geom_segment(ggplot2::aes(x = !!quote(service), xend = !!quote(leave),
                                       y = !!quote(K + GAP + facility), yend = !!quote(K + GAP + facility),
                                       colour = 'Service'), size = SIZE, data = PLOTDATA.REDUCED) +
    ggplot2::geom_segment(ggplot2::aes(x = !!quote(leave), xend = !!quote(revive),
                                       y = !!quote(K + GAP + facility), yend = !!quote(K + GAP + facility),
                                       colour = 'Revival'), size = SIZE, data = PLOTDATA.REDUCED) +
    ggplot2::scale_colour_manual(name = '', values = COLOURS) +
    ggplot2::scale_y_reverse(breaks = c(1:K, (K+GAP+1):(K+GAP+n)),
                             labels = c(as.character(1:K), sprintf('F[%s]', 1:n)),
                             limits = c(K+GAP+n, 1)) +
    ggplot2::theme(plot.title    = ggplot2::element_text(hjust = 0.5, size = 14, face = 'bold',
                                     margin = ggplot2::margin(t = 10, r = 0, b = 10, l = 0)),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = 'bold',
                                     margin = ggplot2::margin(t = 0, r = 0, b = 10, l = 0)),
                   legend.position = 'bottom', legend.box = 'vertical',
                   legend.title = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_text(
                                    margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                   panel.grid.minor.y = ggplot2::element_blank()) +
    ggplot2::labs(title = 'Queuing Plot', subtitle = SUBTITLE, x = 'Time', y = NULL)

  #Add lines for closure-times (unless these are infinite)
  if (T1 < Inf) {
    LINE1   <- data.frame(xintercept = T1, label = factor('Close-time (new arrivals)'))
    PLOT <- PLOT + ggplot2::geom_vline(ggplot2::aes(xintercept = !!quote(xintercept),
                                                    linetype = !!quote(label)), LINE1) }
  if (T2 < Inf) {
    LINE2   <- data.frame(xintercept = T2, label = factor('Close-time (new service)'))
    PLOT <- PLOT + ggplot2::geom_vline(ggplot2::aes(xintercept = !!quote(xintercept),
                                                    linetype = !!quote(label)), LINE2) }
  if (T3 < Inf) {
    LINE3   <- data.frame(xintercept = T3, label = factor('Close-time (all service)'))
    PLOT <- PLOT + ggplot2::geom_vline(ggplot2::aes(xintercept = !!quote(xintercept),
                                                    linetype = !!quote(label)), LINE3) }

  #Print/return the plot
  if (print) { print(PLOT) }
  PLOT }

#' @rdname queue
#' @param probs summary quantiles to be included in output.
#' @param probs.decimal.places rounds the output to specified number of decimal places.
summary.queue <- function(object, probs = NULL, probs.decimal.places = 2, ...) {

  #Check input
  if (!inherits(object, 'queue'))     stop('Error: Object must have class \'queue\'' )
  if (is.null(probs)) {
    probs <- c(0, 0.25, 0.5, 0.75, 1)
  } else {
    if (!is.numeric(probs))              stop('Error: Input probs must be a vector of probabilities')
    if (min(probs) < 0)                   stop('Error: Input probs must be a vector of probabilities')
    if (max(probs) > 1)                   stop('Error: Input probs must be a vector of probabilities') }
  if (!is.numeric(probs.decimal.places)) stop('Error: Input probs.decimal.places should be a non-negative integer')
  if (length(probs.decimal.places) != 1) stop('Error: Input probs.decimal.places should be a non-negative integer')
  if (as.integer(probs.decimal.places) != probs.decimal.places)
                                         stop('Error: Input probs.decimal.places should be a non-negative integer')
  if (min(probs.decimal.places) < 0)     stop('Error: Input probs.decimal.places should be a non-negative integer')

  #Extract object information
  USER     <- object$users
  FACILITY <- object$facilities
  n        <- object$n
  r        <- object$revival
  T1       <- object$close.arrive
  T2       <- object$close.service
  T3       <- object$close.full

  #Compute user summary statistics
  mean.wait      <- mean(USER$wait)
  mean.use       <- mean(USER$use)
  mean.unserved  <- mean(USER$unserved)
  mean.use.full  <- mean(USER$use.full)
  mean.use.prop  <- mean(USER$use.prop)
  sd.wait        <- sd(USER$wait)
  sd.use         <- sd(USER$use)
  sd.unserved    <- sd(USER$unserved)
  sd.use.full    <- sd(USER$use.full)
  sd.use.prop    <- sd(USER$use.prop)
  quant.wait     <- quantile(USER$wait,     probs = probs)
  quant.use      <- quantile(USER$use,      probs = probs)
  quant.unserved <- quantile(USER$unserved, probs = probs)
  quant.use.full <- quantile(USER$use.full, probs = probs)
  quant.use.prop <- quantile(USER$use.prop, probs = probs)
  WAIT.STATS     <- c(mean.wait,     sd.wait,     quant.wait)
  USE.STATS      <- c(mean.use,      sd.use,      quant.use)
  UNSERVED.STATS <- c(mean.unserved, sd.unserved, quant.unserved)
  USE.FULL.STATS <- c(mean.use.full, sd.use.full, quant.use.full)
  USE.PROP.STATS <- c(mean.use.prop, sd.use.prop, quant.use.prop)

  #Generate service facility statistics
  mean.useF    <- mean(FACILITY$use)
  mean.revive  <- mean(FACILITY$revive)
  sd.useF      <- sd(FACILITY$use)
  sd.revive    <- sd(FACILITY$revive)
  quant.useF   <- quantile(FACILITY$use,    probs = probs)
  quant.revive <- quantile(FACILITY$revive, probs = probs)
  USEF.STATS   <- c(mean.useF, sd.useF, quant.useF)
  REVIVE.STATS <- c(mean.revive, sd.revive, quant.revive)

  #Generate summary statistics table
  STATS <- data.frame(wait = WAIT.STATS, use = USE.STATS, unserved = UNSERVED.STATS,
                      use.full = USE.FULL.STATS, use.prop = USE.PROP.STATS,
                      useF = USEF.STATS, revive = REVIVE.STATS)
  rownames(STATS) <- c('Mean', 'Std.Dev', sprintf('Quantile[%s]',
                     format(round(probs, probs.decimal.places),
                     nsmall = probs.decimal.places)))

  #Generate output
  OUT <- list(summary.stats = STATS, users = USER, facilities = FACILITY, n = n, revival = r,
              close.arrive = T1, close.service = T2, close.full = T3)
  class(OUT) <- c('summary.queue', 'list')

  #Return output
  OUT }

#' @rdname queue
print.summary.queue <- function(x, ...) {

  #Check input
  if (!('summary.queue' %in% class(x))) stop('Error: Object must have class \'summary.queue\'' )

  #Extract object information
  STATS <- x$summary.stats
  n     <- x$n
  r     <- x$revival
  T1    <- x$close.arrive
  T2    <- x$close.service
  T3    <- x$close.full

  #Print heading
  PLURAL <- ifelse(n > 1, 'facilities', 'facility')
  cat('\n    Summary Statistics for a Queuing Process \n \n')
  if (r == 0) {
    cat(paste0('Model of an amenity with ', n, ' service ', PLURAL, ' \n'))
  } else {
    cat(paste0('Model of an amenity with ', n, ' service ', PLURAL, ' with revival-time ',
               round(r, 4), '\n')) }
  if (T1 < Inf) { cat('Service facilities close to new arrivals at closure-time =', T1, '\n') }
  if (T1 < Inf) { cat('Service facilities close to new services at closure-time =', T2, '\n') }
  if (T1 < Inf) { cat('Service facilities end existing services at closure-time =', T3, '\n') }
  cat('\n')
  cat('Users are allocated to facilities on a \'first-come first-served\' basis \n \n')

  #Print summary statistics
  cat('---------------------------------------------------------------------- \n \n')
  cat('\n')
  print(STATS[, 1:5])
  cat('\n') }

#' @rdname queue
#' @param count absolute or relative frequencies
#' @param bar.colors,bar.colours plotting parameters
plot.summary.queue <- function(x, print = TRUE, count = FALSE,
                               bar.colors = NULL, bar.colours = bar.colors, ...) {

  #Check input
  if (!('summary.queue' %in% class(x))) stop('Error: Object must have class \'summary.queue\'')
  if (!is.logical(print))                stop('Error: Input print should be a logical value')
  if (length(print) != 1)                stop('Error: Input print should be a single logical value')
  if (!is.logical(count))                stop('Error: Input count should be a logical value')
  if (length(count) != 1)                stop('Error: Input count should be a single logical value')
  if (is.null(bar.colours)) {
    COLOURS <- c('Waiting-time' = 'Red', 'Use-time' = 'Blue', 'Unserved-time' = 'Black')
  } else {
    if (!is.character(bar.colours))     stop('Error: Input bar.colours must be a character vector with three colours')
    if (length(bar.colours) != 4)       stop('Error: Input bar.colours must be a character vector with three colours')
    if (!(all(tolower(bar.colours) %in% colours())))
      stop('Error: Input bar.colours must be a character vector with three colours')
    COLOURS <- c('Waiting-time'  = bar.colours[1],
                 'Use-time' = bar.colours[2],
                 'Unserved-time' = bar.colours[3]) }

  #Extract object information
  USER     <- x$users
  K        <- nrow(USER)
  n        <- x$n
  r        <- x$revival

  #Create subtitle
  PLURAL   <- ifelse(n > 1, 'facilities', 'facility')
  if (r == 0) {
    SUBTITLE <- paste0('(Queuing process with ', n, ' service ', PLURAL, ')')
  } else {
    SUBTITLE <- paste0('(Queuing process with ', n, ' service ', PLURAL,
                       ' with revival-time ', round(r, 4), ')') }

  #Generate the plot data
  LABELS <- c('Waiting-time', 'Use-time', 'Unserved-time')
  TYPE   <- rep(c('Waiting', 'Service', 'Unserved'), each = K)
  PLOTDATA <- data.frame(Value = c(USER$wait, USER$use, USER$unserved),
                         Type  = factor(TYPE, labels = LABELS,
                                        levels = c('Waiting', 'Service', 'Unserved')))

  #Generate plot
  PLOT <- ggplot2::ggplot(ggplot2::aes(x = !!quote(Value), fill = !!quote(Type)), data = PLOTDATA) +
          { if (count)  ggplot2::geom_histogram(ggplot2::aes(y = !!quote(..count..)), bins = 50) } +
          { if (!count) ggplot2::geom_histogram(ggplot2::aes(y = !!quote(..count../sum(..count..))), bins = 50) } +
          ggplot2::facet_wrap(~ Type, ncol = 1) +
          ggplot2::scale_fill_manual(name = '', values = COLOURS) +
          ggplot2::theme(plot.title    = ggplot2::element_text(hjust = 0.5, size = 14, face = 'bold',
                                           margin = ggplot2::margin(t = 10, r = 0, b = 10, l = 0)),
                         plot.subtitle = ggplot2::element_text(hjust = 0.5, face = 'bold',
                                           margin = ggplot2::margin(t = 0, r = 0, b = 10, l = 0)),
                         legend.position = 'none',
                         axis.title.x = ggplot2::element_text(
                           margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                         axis.title.y = ggplot2::element_text(
                           margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0))) +
          ggplot2::labs(title = 'Queuing Summary Plot', subtitle = SUBTITLE,
                        x = 'Time', y = ifelse(count, 'Count', 'Proportion'))

  #Print/return the plot
  if (print) { print(PLOT) }
  PLOT }

