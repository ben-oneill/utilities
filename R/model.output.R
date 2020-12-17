#' Model Output
#'
#' \code{model.output} returns outputs for a statistical model
#'
#' This function produces standard output tables for statistical models, including the summary statistics table, the coefficient estimates
#' table and the ANOVA table.  The output has a custom print method that prints these tables in a user-friendly format.
#'
#'
#' @usage \code{model.output}
#' @param model An object of class \code{'lm'}, \code{'nlm'}, \code{'glm'},etc.
#' @param digits The number of digits to print
#' @return A list of useful model outputs

model.output <- function(model, digits = 4, ...) {

  #Check inputs
  if (!('lm'  %in% class(model))) {
  if (!('nlm' %in% class(model))) {
  if (!('glm' %in% class(model))) {              stop('Error: Input model should be a model of class lm, nlm or glm') } } }

  #Construct description table
  CALL <- as.character(model$call)
  if (CALL[1] == 'lm')     { CALL[1] <- 'Linear Regression' }
  if (CALL[1] == 'nlm')    { CALL[1] <- 'Nonlinear Regression' }
  if (CALL[1] == 'glm')    { CALL[1] <- 'Generalised Linear Model' }
  DESC <- data.frame(Model   = CALL[1],
                     Formula = CALL[2],
                     Data    = CALL[3])
  rownames(DESC) <- ''

  #Construct summary statistics table
  SUMMARY <- summary(model)
  STATS   <- data.frame(Observations      = as.integer(nrow(model$model)),
                        `Standard Error`  = SUMMARY$sigma,
                        `Multiple R`      = sqrt(SUMMARY$r.squared),
                        `R-Squared`       = SUMMARY$r.squared,
                        `R-Squared (Adj)` = SUMMARY$adj.r.squared,
                        check.names = FALSE)
  rownames(STATS) <- ''

  #Construct coefficient estimates table
  COEF <- as.data.frame(SUMMARY$coef)
  colnames(COEF) <- c('Estimate', 'Std Error', 'T', 'p-value')

  #Compute ANOVA quantities
  TABLE  <- anova(model, ...)
  m      <- nrow(TABLE)-1
  DF.reg <- sum(TABLE[1:m, 1])
  DF.res <- TABLE[m+1, 1]
  DF.tot <- DF.reg + DF.res
  SS.reg <- sum(TABLE[1:m, 2])
  SS.res <- TABLE[m+1, 2]
  SS.tot <- SS.reg + SS.res
  MS.reg <- SS.reg/DF.reg
  MS.res <- SS.res/DF.res
  MS.tot <- SS.tot/DF.tot
  FSTAT  <- MS.reg/MS.res
  PVALUE <- pf(FSTAT,  df1 = DF.reg, df2 = DF.res, lower.tail = FALSE)

  #Construct ANOVA table
  ANOVA <- data.frame(Df        = c(DF.reg, DF.res, DF.tot),
                      `Sum Sq`  = c(SS.reg, SS.res, SS.tot),
                      `Mean Sq` = c(MS.reg, MS.res, MS.tot),
                      `F`       = c(FSTAT, NA, NA),
                      `p-value` = c(PVALUE, NA, NA),
                      check.names = FALSE)
  rownames(ANOVA) <- c('Regression', 'Residuals', 'Total')

  #Set class
  OUT <- list(Description = DESC, Summary = STATS, Coefficients = COEF, ANOVA = ANOVA)
  class(OUT) <- 'model.output'
  attr(OUT, 'model')  <- class(model)[1]
  attr(OUT, 'digits') <- digits

  #Give the output
  OUT }


print.model.output <- function(object) {

  DIGITS <- attributes(object)$digits

  cat('\n----------------------------------------------------------------------\n Model \n \n')
  print(object$Description)
  cat('\n----------------------------------------------------------------------\n Summary Statistics Table \n \n')
  print(round(object$Summary, digits = DIGITS))
  cat('\n----------------------------------------------------------------------\n Coefficient Estimates Table \n \n')
  print.data.frame(object$Coefficients, digits = DIGITS)
  cat('\n----------------------------------------------------------------------\n ANOVA Table \n \n')
  print(as.matrix(object$ANOVA), na.print = "", quote = FALSE)
  cat('\n----------------------------------------------------------------------\n') }

