#' Examine mappings between factor variables in a data-frame
#'
#' \code{mappings} determines the mappings between factor variables in a data-frame
#'
#' In preliminary data analysis prior to statistical modelling, it is often useful to investigate whether there are mappings between factor
#' variables in a data-frame in order to see if any of these factor variables are redundant (i.e., fully determined by other factor variables).
#' This function takes an input data-frame \code{data} and examines whether there are any mappings between the factor variables.  (Note that
#' the function will interpret all character variables as factors but will not interpret numeric or logical variables as factors.) The output
#' is a list showing the uniqueness of the binary relations between the factor variables (a logical matrix showing left-uniqueness in the binary
#' relations), the mappings between factor variables, the redundant and non-redundant factor variables, and the directed acyclic graph (DAG) of
#' these mappings (the last element requires the user to have the \code{ggdag} package installed; it is omitted if the package is not installed).
#' If \code{plot = TRUE} the function also returns a plot of the DAG (if \code{ggdag} and \code{ggplot2} packages are installed).
#'
#' Note that the function also allows the user to examine mappings between all variables in the data-frame (i.e., not just the factor variables)
#' by setting \code{all.vars = TRUE}.  The output from this analysis should be interpreted with caution; one-to-one mappings between non-factor variables
#' are common (e.g., when two variables are continuous it is almost certain that they will be in a one-to-one mapping), and so the existence of a
#' mapping may not be indicative of variable redundancy.
#'
#' Note on operation: If \code{na.rm = FALSE} then the function analyses the mappings between the factors/variables without removing NA values.  In
#' this case an \code{NA} value is treated as a missing value that could be any outcome.  Consequently, for purposes of determining whether there
#' is a mapping between the variables, an \code{NA} value is treated as if it were every possible value.  The mapping is falsified if there are at
#' least two identical values in the domain (which may include one or more \code{NA} values) that map to different values in the codomain (which
#' may include one or more \code{NA} values).
#'
#' @aliases print.data.mappings
#'
#' @param data A data-frame (or an object coercible to a data-frame)
#' @param na.rm Logical value; if \code{TRUE} the function removes \code{NA} values from consideration
#' @param all.vars Logical value; if \code{TRUE} the function only examines factor variables in the data-frame; if \code{FALSE} the function
#' examines all variables in the data-frame (caution is required in interpretation of output)
#' @param plot Logical value; if \code{TRUE} the function plots the DAG for the mappings (requires \code{ggplot2} and \code{ggdag} to work)
#' @return A list object of class 'mappings' giving information on the mappings between the variables
#'
#' @examples
#'
#' DATA <- data.frame(
#'   VAR1 = c(0,1,2,2,0,1,2,0,0,1),
#'   VAR2 = c('A','B','B','B','A','B','B','A','A','B'),
#'   VAR3 = 1:10,
#'   VAR4 = c('A','B','C','D','A','B','D','A','A','B'),
#'   VAR5 = c(1:5,1:5)
#' )
#'
#' # Apply mappings
#' mappings(DATA, all.vars = TRUE, plot = FALSE)

mappings <- function(data, na.rm = TRUE, all.vars = FALSE, plot = TRUE) {

  #Extract data name
  DATANAME <- deparse(substitute(data))

  #Check input data
  DATA <- data <- as.data.frame(data)
  if (!is.data.frame(DATA))                                       { stop('Error: Input is not a data frame') }

  #Check input na.rm
  if (!is.logical(na.rm))                                         { stop('Error: Input na.rm should be a logical value') }
  if (length(na.rm) != 1)                                         { stop('Error: Input na.rm should be a single logical value') }

  #Check input all.vars
  if (!is.logical(all.vars))                                      { stop('Error: Input all.vars should be a logical value') }
  if (length(all.vars) != 1)                                      { stop('Error: Input all.vars should be a single logical value') }

  #Check input plot
  if (!is.logical(plot))                                          { stop('Error: Input plot should be a logical value') }
  if (length(plot) != 1)                                          { stop('Error: Input plot should be a single logical value') }

  #Record the factors in the data-frame and reduce to factors (if needed)
  n <- ncol(data)
  FACTORS <- rep(FALSE, n)
  for (i in 1:n) { FACTORS[i] <- ((is.factor(data[, i]))|(is.character(data[, i]))) }
  if (!all.vars) {
    data    <- data[, FACTORS]
    FACTORS <- FACTORS[FACTORS]
    n       <- ncol(data) }

  #Examine uniqueness relations and create output list
  NAMES <- colnames(data)
  UU    <- array(NA, dim = c(n,n))
  rownames(UU) <- NAMES
  colnames(UU) <- NAMES
  for (i in 1:n) {
  for (j in 1:n) {
    II <- data[,i]
    JJ <- data[,j]
    LU <- TRUE
    if (na.rm) {
      REMOVE <- which(is.na(II)|is.na(JJ))
      if (length(REMOVE) > 0) { II <- II[-REMOVE]; JJ <- JJ[-REMOVE] } }
    while ((length(II) > 0) & LU) {
      IND    <- which(II == II[1])
      IND.NA <- which(is.na(II))
      if (length(unique(c(JJ[IND], JJ[IND.NA]))) > 1) { LU <- FALSE }
      II <- II[-IND]; JJ <- JJ[-IND] }
    UU[i,j] <- LU } }
  attr(UU, 'factor') <- FACTORS
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
  for (i in 1:n) { RED[i] <- (sum(UU[FACTORS ,i]) > 0) }
  OUT$factor.redundant    <- colnames(UU)[FACTORS&RED]
  OUT$factor.nonredundant <- colnames(UU)[FACTORS&(!RED)]
  OUT$nonfactor           <- colnames(UU)[!FACTORS]

  #Add class and attributes
  class(OUT) <- c('data.mappings', 'list')
  attr(OUT, 'na.rm')    <- na.rm
  attr(OUT, 'all.vars') <- all.vars
  attr(OUT, 'dataname') <- DATANAME

  #Print plot if requested
  if (plot) { plot.data.mappings(OUT) }

  #Return output
  OUT }

print.data.mappings <- function(x, ...) {

  object <- x
  #Extract information
  DATANAME <- attributes(object)$dataname
  FACTORS  <- attributes(object$uniqueness)$factor
  ALLFACS  <- all(FACTORS)
  NA.RM    <- attributes(object)$na.rm
  RED  <- object$factor.redundant
  NRED <- object$factor.nonredundant
  VARS <- object$nonfactor
  MM   <- object$mappings
  m    <- length(MM)
  n0   <- length(RED)
  n1   <- length(NRED)
  n2   <- length(VARS)
  n    <- n0 + n1

  #Print title
  cat('\n    Mapping analysis for data-frame', DATANAME, 'containing ')
  TYPE1 <- min(n,  2)
  TYPE2 <- min(n2, 2)
  if (TYPE1 == 1) { cat(n, 'factor ')  }
  if (TYPE1 == 2) { cat(n, 'factors ') }
  if ((TYPE1  > 0)&(TYPE2  > 0)) { cat('and ') }
  if (TYPE2 == 1) { cat(n2, 'non-factor variable ')  }
  if (TYPE2 == 2) { cat(n2, 'non-factor variables ')  }
  if (NA.RM) { cat('(analysis ignores NA values) \n') }
  cat('\n')

  #Print results of mapping function
  if (m == 0)   { cat('There were no mappings identified \n \n') } else {
    if (m == 1) { cat('There was one mapping identified:  \n \n') } else {
      cat('There were', m, 'mappings identified: \n \n') } }

  #Print the mappings
  if (m > 0){
    ARROW.POSITION <- rep(0, m)
    for (i in 1:m) { ARROW.POSITION[i] <- which(strsplit(MM[i], '')[[1]] == '\U2192') }
    MAX.POSITION   <- max(ARROW.POSITION)
    for (i in 1:m) {
      START.SPACE <- strrep(' ', times = MAX.POSITION - ARROW.POSITION[i])
      cat('   ', START.SPACE, MM[i], '\n') }
    cat('\n') }

  #Print the redundant factors
  if (n0 > 0) {
    cat('Redundant factors: \n \n')
    for (i in 1:n0) { cat('    ', RED[i], '\n') }
    cat('\n') }

  #Print the nonredundant factors
  if (n1 > 0) {
    cat('Non-Redundant factors: \n \n')
    for (i in 1:n1) { cat('    ', NRED[i], '\n') }
    cat('\n') }

  #Print the other variables
  if (n2 > 0) {
    cat('Other variables: \n \n')
    for (i in 1:n2) { cat('    ', VARS[i], '\n') }
    cat('\n') } }

#' Plot components from data mapping
#'
#' This needs \code{ggplot2} and \code{ggdag} to function correctly.
#'
#' @param x a data mapping
#' @param node.size node size
#' @param text.size label size for a node
#' @param line.width line width
#' @param ... not used
#' @returns nothing
#'
plot.data.mappings <- function(x, node.size = 1, text.size = 1, line.width = 1, ...) {

  object <- x

  #Check inputs
  if (!('data.mappings' %in% class(object)))                           { stop('Error: This print method is for data.mappings objects') }
  if (!is.vector(node.size))                                           { stop('Error: Input node.size should be a single positive value') }
  if (!is.numeric(node.size))                                          { stop('Error: Input node.size should be a single positive value') }
  if (length(node.size) != 1)                                          { stop('Error: Input node.size should be a single positive value') }
  if (node.size <= 0)                                                  { stop('Error: Input node.size should be a single positive value') }
  if (!is.vector(text.size))                                           { stop('Error: Input text.size should be a single positive value') }
  if (!is.numeric(text.size))                                          { stop('Error: Input text.size should be a single positive value') }
  if (length(text.size) != 1)                                          { stop('Error: Input text.size should be a single positive value') }
  if (text.size <= 0)                                                  { stop('Error: Input text.size should be a single positive value') }
  if (!is.vector(line.width))                                          { stop('Error: Input line.width should be a single positive value') }
  if (!is.numeric(line.width))                                         { stop('Error: Input line.width should be a single positive value') }
  if (length(line.width) != 1)                                         { stop('Error: Input line.width should be a single positive value') }
  if (line.width <= 0)                                                 { stop('Error: Input line.width should be a single positive value') }

  #Extract information
  DATANAME <- attributes(object)$dataname
  FACTORS  <- attributes(object$uniqueness)$factor
  ALL.VARS <- attributes(object)$all.vars
  NA.RM    <- attributes(object)$na.rm
  MM       <- object$mappings
  m        <- length(MM)

  if (m == 0) {

    message('There were no mappings to print')

  } else {

    #Check installed packages and load them
    GGDAG   <- requireNamespace('ggdag',   quietly = TRUE)
    GGPLOT2 <- requireNamespace('ggplot2', quietly = TRUE)
    if (!GGDAG)   { stop('Error: Plotting a data.mappings object requires the ggdag package')   }
    if (!GGPLOT2) { stop('Error: Plotting a data.mappings object requires the ggplot2 package') }

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
    DAG$data$Type <- 'Other Variable'
    for (i in 1:length(nnn)) {
      if (nnn[i] %in% object$factor.redundant)    { DAG$data$Type[i] <- 'Redundant Factor' }
      if (nnn[i] %in% object$factor.nonredundant) { DAG$data$Type[i] <- 'Non-Redundant Factor' } }

    #Create the plot
    SUBTITLE <- ifelse(ALL.VARS, ifelse(NA.RM, paste0('(Mappings for all variables in data-frame ', DATANAME, ' --- Analysis ignores NA values) \n'),
                                               paste0('(Mappings for all variables in data-frame ', DATANAME, ') \n')),
                                 ifelse(NA.RM, paste0('(Mappings for all factor variables in data-frame ', DATANAME, ' --- Analysis ignores NA values) \n'),
                                               paste0('(Mappings for all variables in data-frame ', DATANAME, ') \n')))
    PLOT <- ggplot2::ggplot(ggplot2::aes_string(x = "x", y = "y", xend = "xend", yend = "yend", colour = "Type"), data = DAG) +
            ggdag::geom_dag_point(size = 20*sqrt(node.size)) +
            ggplot2::scale_colour_manual(breaks = c('Redundant Factor', 'Non-Redundant Factor', 'Other Variable'),
                                         values = c('green3', 'darkgreen', 'cornflowerblue')) +
            ggdag::geom_dag_edges(edge_width = 0.5*line.width) +
            ggdag::geom_dag_text(size = 4*sqrt(text.size), colour = 'white') +
            ggdag::theme_dag() +
            ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 10))) +
            ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = 'bottom') +
            ggplot2::theme(plot.title    = ggplot2::element_text(hjust = 0.5, size = 18, face = 'bold'),
                           plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, face = 'bold')) +
            ggplot2::ggtitle('Variable Mapping Plot') +
            ggplot2::labs(subtitle = SUBTITLE)

    #Print the plot
    return(PLOT) } }
