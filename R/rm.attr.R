#' Remove (non-protected) attributes from an object
#'
#' \code{rm.attr} removes (non-protected) attributes from an object
#'
#' This function removes non-protected attributes from an \code{R} object.  If the object is a list then
#' the function will remove attributes within elements of the list down to the level specified by the
#' \code{list.levels} input.  (By default the function removes attributes from all levels of lists.)  If
#' you do not want to remove attributes from elements of a list (but still remove attributes from the outer
#' level) you can set \code{list.levels = 0} to do this..
#'
#' @param object An object to operate on attributes from the object
#' @param list.levels A non-negative integer specifying the number of levels of lists to apply the removal to
#' @param protected A character vector containing the names of protected attributes (not to be removed)
#' @return The object is returned with non-protected attributes removed
#'
#' @examples
#' a <- structure(list(structure(1, x=2, names=3),
#'                list(0, structure(3, x=4, names=5))),
#'                x=3, names = 4)
#' str(rm.attr(a, 1))
#'
rm.attr <- function(object, list.levels = Inf,
                    protected = c('class', 'dim', 'names', 'dimnames', 'rownames', 'colnames')) {

  #Check inputs
  if (!is.numeric(list.levels))             stop('Error: Input list.levels should be a single non-negative integer')
  if (length(list.levels) != 1)             stop('Error: Input list.levels should be a single non-negative integer')
  if (list.levels < 0)                      stop('Error: Input list.levels should be a single non-negative integer')
  if (!is.character(protected))             stop('Error: Input protected should be a character vector')

  #Remove unprotected attributes from object
  ATTR.NAMES <- names(attributes(object))
  OUT <- object

  attributes(OUT)[setdiff(ATTR.NAMES, protected)] <- NULL

  #Recursively remove unprotected attributes from elements
  if (list.levels > 0 && is.list(OUT)) {
      OUT[] <- lapply(OUT, rm.attr, list.levels = list.levels - 1, protected=protected) }

  #Return output
  OUT }
