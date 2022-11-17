#' Convert integers to their textual rendition (in English)
#'
#' \code{integertext} returns a character vector with the textual rendition of the input integers
#'
#' This function converts a vector of integer values into their textual rendition in the English language.
#' The function can accommodate both positive and negative integers.  The function includes some options for
#' the textual rendition of numbers, including an option to add dashes in the textual rendition, an option
#' to set the appropriate descriptor for a negative number (e.g., \'negative\' or \'minus\'), and an option
#' to capitalise the first letter of the textual rendition.  The output is a character vector with the textual
#' renditions of the input numbers.
#'
#' @usage \code{integertext(x)}
#' @param x A vector of integers
#' @param negative.text A character representing the prefix used to denote negative values
#' @param dashes Logical value; if \code{TRUE} the textual rendition of numbers includes dashes
#' @param capitalise.first Logical value; if \code{TRUE} the first letter of the text is capitalised
#' @return A character vector giving the textual rendition of the integers

integertext <- function(x, negative.text = 'negative', dashes = FALSE, capitalise.first = TRUE) {

  #Check input x
  if (!is.vector(x))                       stop('Error: Input x should be a vector')
  if (!is.numeric(x))                      stop('Error: Input x should be a numeric vector')
  XX <- round(x, 0)
  NEGATIVE <- (XX < 0)
  if (any(x != XX))                        stop('Error: Input x should contain only integers')

  #Check input negative.text
  if (!is.vector(negative.text))           stop('Error: Input negative.text should be a vector')
  if (!is.character(negative.text))        stop('Error: Input negative.text should be a character string')
  if (length(negative.text) != 1)          stop('Error: Input negative.text should be a single character string')

  #Check input dashes
  if (!is.vector(dashes))                  stop('Error: Input dashes should be a vector')
  if (!is.logical(dashes))                 stop('Error: Input dashes should be a logical value')
  if (length(dashes) != 1)                 stop('Error: Input dashes should be a single logical value')

  #Check input capitalise.first
  if (!is.vector(capitalise.first))        stop('Error: Input capitalise.first should be a vector')
  if (!is.logical(capitalise.first))       stop('Error: Input capitalise.first should be a logical value')
  if (length(capitalise.first) != 1)       stop('Error: Input capitalise.first should be a single logical value')

  ###############################################################################################################
  ######################################## SET PRELIMINARY TEXTUAL VECTORS ######################################
  ###############################################################################################################

  #Set text for prefixes and large numbers
  if (dashes) {
    PREFIX <-  c('', 'one-hundred-and-', 'two-hundred-and-', 'three-hundred-and-', 'four-hundred-and-',
                 'five-hundred-and-', 'six-hundred-and-', 'seven-hundred-and-', 'eight-hundred-and-', 'nine-hundred-and-')
    NUMORDER <- c('', '-thousand', '-million', '-billion', '-trillion', '-quadrillion', '-quintillion',
                  '-sextillion', '-septillion', '-octillion', '-nonillion', '-decillion', '-undecillion',
                  '-duodecillion', '-tredecillion', '-quattuordecillion', '-quindecillion', '-sexdecillion',
                  '-septendecillion', '-octodecillion', '-novemdecillion')
  } else {
    PREFIX <-  c('', 'one hundred and ', 'two hundred and ', 'three hundred and ', 'four hundred and ',
                 'five hundred and ', 'six hundred and ', 'seven hundred and ', 'eight hundred and ', 'nine hundred and ')
    NUMORDER <- c('', ' thousand', ' million', ' billion', ' trillion', ' quadrillion', ' quintillion',
                  ' sextillion', ' septillion', ' octillion', ' nonillion', ' decillion', ' undecillion',
                  ' duodecillion', ' tredecillion', ' quattuordecillion', ' quindecillion', ' sexdecillion',
                  ' septendecillion', ' octodecillion', ' novemdecillion') }

  #Set text for small numbers
  NUM.99  <- c('',        'one',         'two',         'three',         'four',         'five',         'six',         'seven',         'eight',         'nine',
               'ten',     'eleven',      'twelve',      'thirteen',      'fourteen',     'fifteen',      'sixteen',     'seventeen',     'eighteen',      'nineteen',
               'twenty',  'twenty-one',  'twenty-two',  'twenty-three',  'twenty-four',  'twenty-five',  'twenty-six',  'twenty-seven',  'twenty-eight',  'twenty-nine',
               'thirty',  'thirty-one',  'thirty-two',  'thirty-three',  'thirty-four',  'thirty-five',  'thirty-six',  'thirty-seven',  'thirty-eight',  'thirty-nine',
               'forty',   'forty-one',   'forty-two',   'forty-three',   'forty-four',   'forty-five',   'forty-six',   'forty-seven',   'forty-eight',   'forty-nine',
               'fifty',   'fifty-one',   'fifty-two',   'fifty-three',   'fifty-four',   'fifty-five',   'fifty-six',   'fifty-seven',   'fifty-eight',   'fifty-nine',
               'sixty',   'sixty-one',   'sixty-two',   'sixty-three',   'sixty-four',   'sixty-five',   'sixty-six',   'sixty-seven',   'sixty-eight',   'sixty-nine',
               'seventy', 'seventy-one', 'seventy-two', 'seventy-three', 'seventy-four', 'seventy-five', 'seventy-six', 'seventy-seven', 'seventy-eight', 'seventy-nine',
               'eighty',  'eighty-one',  'eighty-two',  'eighty-three',  'eighty-four',  'eighty-five',  'eighty-six',  'eighty-seven',  'eighty-eight',  'eighty-nine',
               'ninety',  'ninety-one',  'ninety-two',  'ninety-three',  'ninety-four',  'ninety-five',  'ninety-six',  'ninety-seven',  'ninety-eight',  'ninety-nine')
  NUM.999      <- paste(rep(PREFIX, each = 100), rep(NUM.99, times = 10), sep = '')
  if (dashes) {
    NUM.999[101] <- 'one-hundred'
    NUM.999[201] <- 'two-hundred'
    NUM.999[301] <- 'three-hundred'
    NUM.999[401] <- 'four-hundred'
    NUM.999[501] <- 'five-hundred'
    NUM.999[601] <- 'six-hundred'
    NUM.999[701] <- 'seven-hundred'
    NUM.999[801] <- 'eight-hundred'
    NUM.999[901] <- 'nine-hundred'
  } else {
    NUM.999[101] <- 'one hundred'
    NUM.999[201] <- 'two hundred'
    NUM.999[301] <- 'three hundred'
    NUM.999[401] <- 'four hundred'
    NUM.999[501] <- 'five hundred'
    NUM.999[601] <- 'six hundred'
    NUM.999[701] <- 'seven hundred'
    NUM.999[801] <- 'eight hundred'
    NUM.999[901] <- 'nine hundred' }

  ###############################################################################################################
  ########################################## GENERATE TEXTUAL RENDITIONS ########################################
  ###############################################################################################################

  #Extract number information
  #Put numbers in TT form (trinumeral, trinumeral, ...)
  n     <- length(XX)
  XX    <- abs(XX)
  ORDER <- floor(log10(XX))
  TT    <- vector(n, mode = 'list')
  for (i in 1:n) {
    mm <- 3*floor(ORDER[i]/3)
    kk <- 0
    DD <- integer(0)
    while (mm >= 0) {
      VV <- 10^mm
      DD[kk+1] <- floor(XX[i]/VV)
      XX[i] <- XX[i] - DD[kk+1]*VV
      mm <- mm-3
      kk <- kk+1 }
    TT[[i]] <- rev(DD) }
  XX <- round(x, 0)

  #Set text
  TEXTLIST <- vector(n, mode = 'list')
  TEXT     <- rep('', n)
  for (i in 1:n) {
    r <- length(TT[[i]])
    TEXTLIST[[i]] <- NUM.999[TT[[i]]+1]
    LAYER2   <- NUMORDER[1:r]
    LAYER2[TT[[i]] == 0] <- ''
    TEXTLIST[[i]] <- paste(TEXTLIST[[i]], LAYER2, sep = '')
    TEXTLIST[[i]] <- TEXTLIST[[i]][TT[[i]] != 0]
    r <- length(TEXTLIST[[i]])
    if (r > 1) {
      LAYER3 <- c('', rep(', ', r-1))
      if (TT[[i]][1] < 100) { LAYER3[2] <- ' and ' }
      TEXTLIST[[i]] <- paste(TEXTLIST[[i]], LAYER3, sep = '') }
    TEXT[i] <- paste0(rev(TEXTLIST[[i]]), collapse = '')
    if (NEGATIVE[i]) { TEXT[i] <- paste(negative.text, TEXT[i], collapse = ' ') }

    #Special case of zero
    if (XX[i] == 0) { TEXT[i] <- 'zero' }

    #Capitalise (if required)
    if (capitalise.first) {
      FIRST <- substr(TEXT[i], start = 1, stop = 1)
      TEXT[i] <- paste0(toupper(FIRST), substr(TEXT[i], start = 2, stop = nchar(TEXT[i]))) } }

  #Return output
  TEXT }

