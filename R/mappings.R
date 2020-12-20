#' Examine mappings between variables in a data-frame
#'
#' \code{mappings} determines the mappings between variables in a data-frame
#'
#' In preliminary data analysis prior to statistical modelling, it is often useful to investigate whether there are mappings between
#' variables in a data-frame in order to see if any of the variables are redundant (i.e., fully determined by other variables).  This
#' function takes an input data-frame \code{data} and examines whether there are any mappings between the variables.  The output is a
#' list showing the uniqueness of the binary relations between the variables (a logical matrix showing left-uniqueness in the binary
#' relations), the mappings as character strings, the redundant and non-redundant variables, and the directed acyclic graph (DAG) of
#' the mappings (the last element requires the user to have the \code{ggdag} package installed; it is omitted if the package is not
#' installed).  If \code{plot = TRUE} the function also returns a plot of the DAG (if \code{ggdag} and \code{ggplot2} packages are installed).
#'
#' @usage \code{mappings}
#' @param data A data-frame (or an object coercible to a data-frame)
#' @param na.rm Logical value; if \code{TRUE} the function removes \code{NA} values from consideration
#' @param plot Logical value; if \code{TRUE} the function plots the DAG for the mappings (requires \code{ggplot2} and \code{ggdag} to work)
#' @return A list object of class 'mappings' giving information on the mappings between the variables

mappings <- function(data, na.rm = TRUE, plot = TRUE) {

  #Extract data name
  DATANAME <- deparse(substitute(data))

  #Check input data
  data <- as.data.frame(data)
  if (!is.data.frame(DATA))                                       { stop('Error: Input is not a data frame') }

  #Check input na.rm
  if (!is.logical(na.rm))                                         { stop('Error: Input na.rm should be a logical value') }
  if (length(na.rm) != 1)                                         { stop('Error: Input na.rm should be a single logical value') }

  #Check input plot
  if (!is.logical(plot))                                          { stop('Error: Input plot should be a logical value') }
  if (length(plot) != 1)                                          { stop('Error: Input plot should be a single logical value') }

  #Examine uniqueness relations and create output list
  NAMES <- colnames(data)
  n     <- length(NAMES)
  UU    <- array(NA, dim = c(n,n))
  rownames(UU) <- NAMES
  colnames(UU) <- NAMES
  for (i in 1:n) {
  for (j in 1:n) {
    II <- data[,i]
    JJ <- data[,j]
    if (na.rm) { REMOVE <- which(is.na(II)|is.na(JJ))
    if (length(REMOVE) > 0) { II <- II[-REMOVE]; JJ <- JJ[-REMOVE] } }
    LU <- TRUE
    while ((length(II) > 0) & LU) {
      IND <- which(II == II[1]);
      if (length(unique(JJ[IND])) > 1) { LU <- FALSE }
      II <- II[-IND]; JJ <- JJ[-IND] }
    UU[i,j] <- LU } }
  OUT <- list(uniqueness = UU)
  for (i in 1:n) { UU[i,i] <- FALSE }

  #Create vector of mappings
  MAP <- character(0)
  k   <- 0
  for (i in 1:n) {
  for (j in 1:n) {
    if (UU[i,j]) {
      MAP[k+1] <- paste0(NAMES[i], " \U2192 ", NAMES[j])
      k <- k+1 } } }
  OUT$mappings <- MAP

  #Compute redundancies
  RED <- rep(FALSE, n)
  for (i in 1:n) { RED[i] <- (sum(UU[,i]) > 0) }
  OUT$redundant    <- colnames(UU)[RED]
  OUT$nonredundant <- colnames(UU)[!RED]

  #Add class and attributes
  class(OUT) <- c('data.mappings', 'list')
  attr(OUT, 'na.rm')    <- na.rm
  attr(OUT, 'plot')     <- plot
  attr(OUT, 'dataname') <- DATANAME

  #Print plot if requested
  if (plot) { plot.data.mappings(OUT) }

  #Return output
  OUT }


print.data.mappings <- function(object) {

  #Extract information
  DATANAME <- attributes(object)$dataname
  RED  <- object$redundant
  NRED <- object$nonredundant
  MM   <- object$mappings
  m    <- length(MM)
  n0   <- length(RED)
  n1   <- length(NRED)
  n    <- n0 + n1

  #Print results of mapping function
  cat('\n    Mapping analysis for data-frame', DATANAME, 'containing', n, 'variables \n \n')
  if (m == 0)   { cat('There were no mappings between the variables \n \n') } else {
    if (m == 1) { cat('There was one mappings between the variables \n \n') } else {
                  cat('There were', m, 'mappings between the variables \n \n') } }

  #Print the mappings
  if (m > 0){ for (i in 1:m) { cat('    ', MM[i], '\n') } }

  #Print the redundant variables
  if (n0 > 0) {
    cat('\n', 'Redundant variables: \n \n')
    for (i in 1:n0) { cat('    ', RED[i], '\n') } }

  #Print the nonredundant variables
  if (n1 > 0) {
    cat('\n', 'Non-redundant variables: \n \n')
    for (i in 1:n1) { cat('    ', NRED[i], '\n') }
    cat('\n') } }


plot.data.mappings <- function(object) {

  if (length(object$redundant) == 0) {

    message('There were no mappings to print')

  } else {

  #Check installed packages and load them
  GGDAG   <- requireNamespace('ggdag',   quietly = TRUE)
  GGPLOT2 <- requireNamespace('ggplot2', quietly = TRUE)
  if (GGDAG)   { library(ggdag)   } else { stop('Error: Plotting a data.mappings object requires the ggdag package')   }
  if (GGPLOT2) { library(ggplot2) } else { stop('Error: Plotting a data.mappings object requires the ggplot2 package') }

  #Create directed acyclic graph (DAG)
  UU    <- object$uniqueness
  NAMES <- colnames(UU)
  n     <- length(NAMES)
  for (i in 1:n) { UU[i,i] <- FALSE }
  FORMULAS <- list()
  k <- 0
  for (i in 1:n) {
  for (j in 1:n) {
    if (UU[i,j]) {
      k <- k+1
      FORMULAS[[k]] <- as.formula(paste0(NAMES[j], " ~ ", NAMES[i])) } } }
  DAG <- do.call(ggdag::dagify, FORMULAS)
  DAG <- ggdag::tidy_dagitty(DAG)
  nnn <- DAG$data$name
  IND <- rep(FALSE, length(nnn))
  for (i in 1:length(nnn)) { IND[i] <- (nnn[i] %in% object$redundant) }
  DAG$data$Redundant      <- 'Non Redundant Variable'
  DAG$data$Redundant[IND] <- 'Redundant Variable'

  #Create the plot
  PLOT <- ggplot2::ggplot(ggplot2::aes(x = x, y = y, xend = xend, yend = yend, colour = Redundant), data = DAG) +
          ggdag::geom_dag_point() +
          ggplot2::scale_colour_manual(values = c('Black', 'DarkGrey')) +
          ggdag::geom_dag_edges() +
          ggdag::geom_dag_text(colour = 'white') +
          ggdag::theme_dag() +
          ggplot2::theme(legend.position = 'none') +
          ggplot2::theme(plot.title    = ggplot2::element_text(hjust = 0.5, size = 14, face = 'bold'),
                         plot.subtitle = ggplot2::element_text(hjust = 0.5, face = 'bold')) +
          ggplot2::ggtitle('Variable Mappings Plot') +
          ggplot2::labs(subtitle = '(redundant variables are shown in grey)')

  #Print the plot
  plot(PLOT) } }


