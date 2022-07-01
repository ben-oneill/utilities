#' Print a Data-Frame (allowing column/row separators)
#'
#' Custom print method for objects of type \code{data.frame}.  This function prints the data-frame in the same way
#' as the default \link[base]{print.data.frame} in the base package, except that it allows the user to add textual column/row separators in the
#' print output in specified positions.  To do this the user adds row/column values for the inputs \code{row.separator}
#' and \code{col.separator} indicating that separators should be added after those rows/columns.  The user can also
#' set \code{sep.extend} to \code{TRUE} to extend the separators into the row/column-names.
#'
#' @param x A data-frame (object of class \code{data.frame})
#' @param ... optional arguments to \code{print} or \code{plot} methods
#' @param row.separator A vector of values of rows (adds separators after those rows)
#' @param col.separator A vector of values of columns (adds separators after those columns)
#' @param sep.extend Logical value; if \code{TRUE} the separators are extended into the row/column-names
#' @param print.gap A non-negative integer specifyig the number of spaces between columns
#' @param digits the minimum number of significant digits to be used: see \code{print.default}
#' @param quote Logical value; if \code{TRUE} entries are printed with surrounding quotes
#' @param right Logical value; if \code{TRUE} strings are right-aligned
#' @param row.names Logical value or character vector; indicating whether (or what) row names are printed
#' @param max numeric or \code{NULL}, specifying the maximal number of entries to be printed. By default, when \code{NULL}, \code{getOption("max.print")} used
#' @return Prints the data frame with the specified column/row separators

.print.data.frame <- function(x, ...,
                             row.separator = NULL, col.separator = NULL,
                             sep.extend = FALSE, print.gap = 1,
                             digits = NULL, quote = FALSE, right = TRUE,
                             row.names = TRUE, max = NULL) {

  #Get row and column information
  n <- dim(x)[1]
  m <- dim(x)[2]
  ROWNAMES <- rownames(x)
  COLNAMES <- colnames(x)

  #Print output if there are no columns/rows
  if (m == 0L) {
    cat(sprintf(ngettext(n, 'data frame with 0 columns and %d row',
                         'data frame with 0 columns and %d rows'), n), '\n', sep = '') }
  if ((n == 0L)&(m > 0)) {
    print.default(names(x), quote = FALSE)
    cat(gettext('<0 rows> (or 0-length row.names)\n')) }

  #Print output in standard case
  if ((n > 0)&(m > 0)) {

    #Get maximum printing length
    if (is.null(max))       { max <- getOption('max.print', 99999L) }
    if (!is.finite(max))    { stop("invalid 'max' / getOption(\"max.print\"): ", max) }

    #Get number of printing digits
    if (is.null(digits))    { digits <- getOption('digits', 7L) }
    if (!is.finite(digits)) { stop("invalid 'max' / getOption(\"max.print\"): ", max) }

    #Construct character matrix and determine maximum characters in columns
    MAT  <- as.matrix(format.data.frame(x, digits = digits, na.encode = FALSE))
    MAXCHAR <- rep(0, m)
    for (j in 1:m) { MAXCHAR[j] <- max(nchar(MAT[,j])) }

    #Construct new character matrix with room for column separators
    r  <- length(row.separator)
    c  <- length(col.separator)
    SEPMAT <- matrix('', nrow = n, ncol = m+c)
    colnames(SEPMAT) <- rep('', m+c)
    RR <- row.separator + 1:r
    CC <- col.separator + 1:c

    #Determine maximum characters in columns
    MAXCHAR2 <- rep(1, m+c)
    if (c > 0)  { MAXCHAR2[-CC] <- MAXCHAR }
    if (c == 0) { MAXCHAR2      <- MAXCHAR }

    #Add column separators and values to new matrix
    if (c > 0) {
      SEPMAT[, -CC] <- MAT
      colnames(SEPMAT)[-CC] <- COLNAMES }
    if (c == 0) {
      SEPMAT <- MAT
      colnames(SEPMAT) <- COLNAMES }
    if (isTRUE(row.names)) { rownames(SEPMAT) <- ROWNAMES } else {
      rownames(SEPMAT) <- rep('', n) }
    for (j in 1:(m+c)) {
      if (sep.extend) { if (j %in% CC) { colnames(SEPMAT)[j] <- '|' } }
      for (i in 1:n) {
      if (j %in% CC) { SEPMAT[i,j] <- '|' } } }

    #Omit excess rows
    n0   <- max %/% length(x)
    OMIT <- (n0 < n)
    if (OMIT) { SEPMAT <- SEPMAT[seq_len(n0), , drop = FALSE] }

    #Capture the output
    OUT <- capture.output({
      print(SEPMAT, print.gap = print.gap, quote = quote, right = right, max = max);
      if (OMIT) cat(" [ reached 'max' / getOption(\"max.print\") -- omitted", n - n0, "rows ]\n") })

    #Add row separators to output
    if (r > 0) {
      if (isFALSE(row.names)) { M0 <- 0 } else { M0 <- max(nchar(ROWNAMES)) }
      ROWSEP   <- paste(c(rep(' ', M0 + print.gap),
                          rep('-', sum(MAXCHAR2 + print.gap)-1)), collapse = '')
      if (sep.extend) {
        substr(ROWSEP, start = 1, stop = M0 + print.gap) <- paste(rep('-', M0 + print.gap), collapse = '') }
      if (c > 0) {
      for (j in 1:c) {
        INDEX <- 1 + M0 + print.gap + ifelse(CC[j] == 1, 0, sum(MAXCHAR2[1:(CC[j]-1)] + print.gap))
        substr(ROWSEP, start = INDEX, stop = INDEX) <- '+' } }
      OUT2 <- rep('', length(OUT) + r)
      if (r > 0)  { OUT2[-(RR+1)] <- OUT }
      if (r == 0) { OUT2 <- OUT }
      for (i in 1:r) { OUT2[RR[i]+1] <- ROWSEP }
      OUT <- OUT2 }

    #Print the output
    for (i in 1:length(OUT)) {
      cat(OUT[i], '\n') } }

    invisible(x) }
