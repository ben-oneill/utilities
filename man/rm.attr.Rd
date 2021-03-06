% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rm.attr.R
\name{rm.attr}
\alias{rm.attr}
\title{Remove (non-protected) attributes from an object}
\usage{
rm.attr(
  object,
  list.levels = Inf,
  protected = c("class", "dim", "names", "dimnames", "rownames", "colnames")
)
}
\arguments{
\item{object}{An object to operate on attributes from the object}

\item{list.levels}{A non-negative integer specifying the number of levels of lists to apply the removal to}

\item{protected}{A character vector containing the names of protected attributes (not to be removed)}
}
\value{
The object is returned with non-protected attributes removed
}
\description{
\code{rm.attr} removes (non-protected) attributes from an object
}
\details{
This function removes non-protected attributes from an \code{R} object.  If the object is a list then
the function will remove attributes within elements of the list down to the level specified by the
\code{list.levels} input.  (By default the function removes attributes from all levels of lists.)  If
you do not want to remove attributes from elements of a list (but still remove attributes from the outer
level) you can set \code{list.levels = 0} to do this..
}
\examples{
a <- structure(list(structure(1, x=2, names=3),
               list(0, structure(3, x=4, names=5))),
               x=3, names = 4)
str(rm.attr(a, 1))

}
