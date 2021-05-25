#' Remove attributes from an object
#'
#' \code{rm.attr} removes attributes from an object
#'
#' This function removes non-protected attributes from an \code{R} object.  If the object is a list then
#' the function will remove attributes within elements of the list by default, but it can be set to ignore
#' removal of attributes in elements of lists.
#'
#' @usage \code{rm.attr}
#' @param object An object to operate on attributes from the object
#' @param remove.in.list logical; if \code{TRUE} the function removes attributes within elements of lists
#' @param protected A character vector containing the names of protected attributes (not to be removed)
#' @return The object is returned with non-protected attributes removed

rm.attr <- function(object, remove.in.list = TRUE,
                    protected = c('class', 'dim', 'names', 'dimnames', 'rownames', 'colnames')) {

  #Remove unprotected attributes from object
  ATTR.NAMES <- names(attributes(object))
  REMOVE     <- (!(ATTR.NAMES %in% protected))
  OUT <- object
  for(i in ATTR.NAMES[REMOVE]) { attr(OUT, i) <- NULL }

  #Remove unprotected attributes
  if (remove.in.list) {
  if (is.list(OUT)) {
    for (i in 1:length(OUT)) { OUT[[i]] <- rm.attr(OUT[[i]], remove.in.list = TRUE, protected = protected) } } }

  #Return output
  OUT }
