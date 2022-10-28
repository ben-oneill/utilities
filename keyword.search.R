#' Search for keywords in a character vector
#'
#' \code{keyword.search} searches for keywords in a character vector
#'
#' This function searches an input character vector for a set of keywords and gives the result of the search.
#' The user must enter a character vector to search and a vector of keywords; by default the search is case-
#' sensitive and searched only for whole words, but the user can change these search parameters if preferred.
#' The output of the search is a list showing the search that was used and the results of the search.  Results
#' are presented in a data-frame and matrix form showing the keywords found in each element of the character
#' vector that was searched.
#'
#' @usage \code{keyword.search(x, keywords, case.sensitive = TRUE, whole.word = TRUE)}
#' @param text A character vector over which to search
#' @param keywords A character vector containing keywords to search
#' @param case.sensitive Logical value or vector; if \code{TRUE} the search is case-sensitive
#' @param whole.word Logical value or vector; if \code{TRUE} the search is for whole words only
#' @return A list of class \code{'keyword.search'} showing the search parameters and results of the keyword search

keyword.search <- function(text, keywords, case.sensitive = TRUE, whole.word = TRUE) {

  #Check input text
  if (!is.vector(text))                   stop('Error: Input text should be a vector')
  if (!is.character(text))                stop('Error: Input text should be a character vector')
  n <- length(text)

  #Check input keywords
  if (!is.vector(keywords))               stop('Error: Input keywords should be a vector')
  if (!is.character(keywords))            stop('Error: Input keywords should be a character vector')
  m <- length(keywords)
  for (i in 1:m) { if (nchar(keywords)[i] == 0) stop(paste('Error: Keyword', i, 'is blank')) }

  #Check input case.sensitive
  if (!is.vector(case.sensitive))         stop('Error: Input case.sensitive should be a vector')
  if (!is.logical(case.sensitive))        stop('Error: Input case.sensitive should be a logical vector')
  if (!(length(case.sensitive) %in% c(1, m))) {
    stop('Error: Input case.sensitive should either be a single logical value or should have one logical value for each keyword') }
  if (length(case.sensitive) == 1) { CASE <- rep(case.sensitive, m) }
  if (length(case.sensitive) == m) { CASE <- case.sensitive }

  #Check input whole.word
  if (!is.vector(whole.word))         stop('Error: Input whole.word should be a vector')
  if (!is.logical(whole.word))        stop('Error: Input whole.word should be a logical vector')
  if (!(length(whole.word) %in% c(1, m))) {
    stop('Error: Input whole.word should either be a single logical value or should have one logical value for each keyword') }
  if (length(whole.word) == 1) { WHOLE <- rep(whole.word, m) }
  if (length(whole.word) == m) { WHOLE <- whole.word }

  #######################################################################################################

  #Set output parts
  NOTATION   <- sprintf('K[%s]', 1:m)
  OUT.SEARCH <- data.frame(Keyword = keywords, Found = 0, case.sensitive = CASE, whole.word = WHOLE)
  OUT.RESULT <- data.frame(Text = text, Found = 0, Keywords = '')
  OUT.FLAG   <- matrix(FALSE, nrow = n, ncol = m)
  rownames(OUT.FLAG) <- rownames(OUT.RESULT) <- sprintf('Text[%s]', 1:n)
  colnames(OUT.FLAG) <- rownames(OUT.SEARCH) <- NOTATION

  #Populate flagging matrix
  #Value in element (i, k) is set to TRUE if keyword k is in string i
  for (k in 1:m) {

    #Search for whole words
    if (WHOLE[k]) {
      for (i in 1:n) {

        #Remove punctuation and split into distinct words
        STRING <- gsub(pattern = "([[:punct:]])", replacement = " ", text[i])
        WORDS  <- strsplit(STRING, split = ' ')[[1]]

        #Flag keywords
        if (CASE[k])  { OUT.FLAG[i, k] <- (keywords[k] %in% WORDS) }
        if (!CASE[k]) { OUT.FLAG[i, k] <- (tolower(keywords[k]) %in% tolower(WORDS)) } } }

    #Search for partial words
    if (!WHOLE[k]) {
      for (i in 1:n) {
        OUT.FLAG[i, k] <- grepl(x = text[i], pattern = keywords[k], ignore.case = !CASE[k]) } } }

  #Populate search and results tables
  OUT.SEARCH$Found <- colSums(OUT.FLAG)
  OUT.RESULT$Found <- rowSums(OUT.FLAG)
  for (i in 1:n) {
    if (OUT.RESULT$Found[i] > 0) {
      OUT.RESULT$Keywords[i] <- paste0(NOTATION[OUT.FLAG[i, ]], collapse = ', ') } }

  #Create output
  OUT <- list(search = OUT.SEARCH, result = OUT.RESULT, flag.matrix = OUT.FLAG)
  class(OUT) <- 'keyword.search'

  #Return output
  OUT }


print.keyword.search <- function(object, filter = TRUE) {

  #Check object class
  if (!('keyword.search' %in% class(object)))          stop('Error: This print method is for objects of class \'keyword.search\'')

  #Extract parameters
  OUT.SEARCH <- object$search
  OUT.RESULT <- object$result
  OUT.FLAG   <- object$flag.matrix
  m <- nrow(OUT.SEARCH)
  n <- nrow(OUT.RESULT)

  #Print title
  cat('\n  Keyword Search\n\n')
  cat('Search of', m, 'keywords\n')
  cat('Keywords were found in', sum(OUT.RESULT$Found > 0), 'out of', n, 'strings\n')
  cat('\n')

  #Print search table
  PRINT.DF1 <- cbind(OUT.SEARCH, col1 = '', col2 = '')
  PRINT.DF1 <- PRINT.DF1[, c(1, 5, 2, 6, 3, 4)]
  PRINT.DF1[, 5] <- ifelse(OUT.SEARCH$case.sensitive, '', '*')
  PRINT.DF1[, 6] <- ifelse(OUT.SEARCH$whole.word,     '^', '')
  names(PRINT.DF1) <- c('Keyword', '', 'Found', '', '', '')
  cat('--------------------------------------------------------\n')
  cat('  Keyword Table\n\n')
  print(PRINT.DF1, right = FALSE)
  cat('\n')
  if (sum(!OUT.SEARCH$case.sensitive) > 0) { cat('* Keyword is not case-sensitive\n') }
  if (sum(OUT.SEARCH$whole.word) > 0)      { cat('^ Keyword must be matched to whole word\n') }
  if ((sum(!OUT.SEARCH$case.sensitive) > 0)|(sum(OUT.SEARCH$whole.word) > 0)) { cat('\n') }

  #Print results table
  PRINT.DF2 <- cbind(OUT.RESULT, col1 = '', col2 = '')
  PRINT.DF2 <- PRINT.DF2[, c(1, 4, 2, 5, 3)]
  names(PRINT.DF2) <- c('', ' ', 'Found', ' ', 'Keywords')
  cat('--------------------------------------------------------\n')
  if (filter) {
    cat('  Results Table\n\n')
    FF <- which(PRINT.DF2$Found > 0)
    if (length(FF) > 0)
    print(PRINT.DF2[FF, ], right = FALSE)
  } else {
    cat('  Results Table\n\n')
    print(PRINT.DF2, right = FALSE) }
  cat('\n') }

