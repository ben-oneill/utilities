#' Examine mappings between variables in a data-frame
#'
#' \code{mappings} determines the mappings between variables in a data-frame
#'
#' In preliminary data analysis prior to statistical modelling, it is often useful to investigate whether there are mappings between
#' variables in a data-frame in order to see if any of the variables are redundant (i.e., fully determined by other variables).  This
#' function takes an input data-frame \code{data} and examines whether there are any mappings between the variables.  The output is a
#' list showing the mappings as character strings, the redundant and non-redundant variables, and the directed acyclic graph (DAG) of
#' the mappings (the last element requires the user to have the \code{ggdag} package installed; it is omitted if the package is not
#' installed).  If \code{plot = TRUE} the function also returns a plot of the DAG (if \code{ggdag} and \code{ggplot2} packages are installed).
#'
#' @usage \code{mappings}
#' @param data A data-frame (or an object coercible to a data-frame)
#' @param na.rm Logical value; if \code{TRUE} the function removes \code{NA} values from consideration
#' @param plot Logical value; if \code{TRUE} the function plots the DAG for the mappings (requires \code{ggplot2} and \code{ggdag} to work)
#' @return Logical matrix showing uniqueness relations between variables in data

mappings <- function(data, na.rm = TRUE, plot = TRUE) {

  #Check input data
  data <- as.data.frame(data)
  if (!is.data.frame(DATA))                                       { stop('Error: Input is not a data frame') }

  #Check input na.rm
  if (!is.logical(na.rm))                                         { stop('Error: Input na.rm should be a logical value') }
  if (length(na.rm) != 1)                                         { stop('Error: Input na.rm should be a single logical value') }

  #Check input plot
  if (!is.logical(plot))                                          { stop('Error: Input plot should be a logical value') }
  if (length(plot) != 1)                                          { stop('Error: Input plot should be a single logical value') }

  #Check installed packages and load them
  GGDAG   <- requireNamespace('ggdag',   quietly = TRUE)
  GGPLOT2 <- requireNamespace('ggplot2', quietly = TRUE)
  if (GGDAG)   { library(ggdag)  }
  if (GGPLOT2) { library(ggplot2) }

  #Extract uniqueness relations and create output list
  UU    <- uniqueness(data, na.rm = na.rm)
  NAMES <- colnames(UU)
  n     <- ncol(UU)
  for (i in 1:n) { UU[i,i] <- FALSE }
  OUT <- list()

  #Create vector of mappings
  MAP <- character(0)
  k   <- 0
  for (i in 1:n) {
  for (j in 1:n) {
    if (UU[i,j]) {
      MAP[k+1] <- paste0(NAMES[i], " \U2192 ", NAMES[j])
      k <- k+1 } } }
  OUT$Mappings <- MAP

  #Compute redundancies
  RED <- rep(FALSE, n)
  for (i in 1:n) { RED[i] <- (sum(UU[,i]) > 0) }
  OUT$Redundant     <- colnames(UU)[RED]
  OUT$Non_Redundant <- colnames(UU)[!RED]

  #Create directed acyclic graph (DAG)
  if (GGDAG) {
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
    for (i in 1:length(nnn)) { IND[i] <- (nnn[i] %in% OUT$Redundant) }
    DAG$data$Redundant      <- 'Non Redundant Variable'
    DAG$data$Redundant[IND] <- 'Redundant Variable'
    OUT$DAG <- DAG }

  #Print plot if requested
  if (all(GGDAG, GGPLOT2, plot, length(FORMULAS) > 0)) {
    PLOT <- ggplot2::ggplot(ggplot2::aes(x = x, y = y, xend = xend, yend = yend, colour = Redundant), data = DAG) +
            ggdag::geom_dag_point() +
            ggplot2::scale_colour_manual(values = c('Black', 'DarkGrey')) +
            ggdag::geom_dag_edges() +
            ggdag::geom_dag_text(colour = 'white') +
            ggdag::theme_dag() +
            ggplot2::theme(legend.position = 'none')

    print(PLOT) }

  #Return output
  OUT }


uniqueness <- function(data, na.rm = TRUE) {

  #Check input data
  DATA <- as.data.frame(data)
  if (!is.data.frame(data)) { stop('Error: Input is not a data frame') }

  #Generate output matrix
  NAMES <- colnames(data)
  n     <- length(NAMES)
  OUT   <- array(NA, dim = c(n,n))
  rownames(OUT) <- NAMES
  colnames(OUT) <- NAMES

  #Check mappings
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
      OUT[i,j] <- LU } }

  #Return output
  OUT }

