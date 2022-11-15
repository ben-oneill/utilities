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
#' @usage \code{keyword.search(text, keywords, case.sensitive = TRUE, whole.word = TRUE)}
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


plot.keyword.search <- function(object, subtitle = TRUE,
                                show.text = TRUE, text.max = 20,
                                show.keywords = TRUE, keyword.max = 14,
                                point.size = 3, point.shape = 15,
                                point.color = 'darkblue', point.colour = point.color) {

  #Check object class
  if (!('keyword.search' %in% class(object)))          stop('Error: This print method is for objects of class \'keyword.search\'')

  #Check if required packages are installed
  GGPLOT2   <- requireNamespace('ggplot2',  quietly = TRUE)
  GGTHEMES  <- requireNamespace('ggthemes', quietly = TRUE)
  if (!GGPLOT2)  { stop('Error: This plot method requires the ggplot2 package') }
  if (!GGTHEMES) { stop('Error: This plot method requires the ggthemes package') }

  #Check input subtitle
  if (!is.vector(subtitle))              stop('Error: Input subtitle should be a vector')
  if (!is.logical(subtitle))             stop('Error: Input subtitle should be a logical value')
  if (length(subtitle) != 1)             stop('Error: Input subtitle should be a single logical value')

  #Check input show.text
  if (!is.vector(show.text))             stop('Error: Input show.text should be a vector')
  if (!is.logical(show.text))            stop('Error: Input show.text should be a logical value')
  if (length(show.text) != 1)            stop('Error: Input show.text should be a single logical value')

  #Check input text.max
  if (!is.vector(text.max))              stop('Error: Input text.max should be a vector')
  if (!is.numeric(text.max))             stop('Error: Input text.max should be numeric')
  if (length(text.max) != 1)             stop('Error: Input text.max should be a single integer')
  if (text.max != as.integer(text.max))  stop('Error: Input text.max should be an integer')
  if (min(text.max) < 10)                stop('Error: Input text.max should be at least 10')

  #Check input show.keywords
  if (!is.vector(show.keywords))         stop('Error: Input show.keywords should be a vector')
  if (!is.logical(show.keywords))        stop('Error: Input show.keywords should be a logical value')
  if (length(show.keywords) != 1)        stop('Error: Input show.keywords should be a single logical value')

  #Check input keyword.max
  if (!is.vector(keyword.max))           stop('Error: Input keyword.max should be a vector')
  if (!is.numeric(keyword.max))          stop('Error: Input keyword.max should be numeric')
  if (length(keyword.max) != 1)          stop('Error: Input keyword.max should be a single integer')
  if (keyword.max != as.integer(keyword.max))  stop('Error: Input keyword.max should be an integer')
  if (min(keyword.max) <  4)             stop('Error: Input keyword.max should be at least 4')

  #Check input point.size
  if (!is.vector(point.size))            stop('Error: Input point.size should be a vector')
  if (!is.numeric(point.size))           stop('Error: Input point.size should be numeric')
  if (length(point.size) != 1)           stop('Error: Input point.size should be a single number')
  if (min(point.size) <= 0)              stop('Error: Input point.size should be greater than zero')

  #Check input point.shape
  if (!is.vector(point.shape))           stop('Error: Input point.shape should be a vector')
  if (!is.numeric(point.shape))          stop('Error: Input point.shape should be numeric')
  if (length(point.shape) != 1)          stop('Error: Input point.shape should be a single number')
  if (!(point.shape %in% 0:25))          stop('Error: Input point.shape should be greater than zero')

  #Check input point.colour
  if (!is.vector(point.colour))          stop('Error: Input point.colour should be a vector')
  if (!is.character(point.colour))       stop('Error: Input point.colour should be numeric')
  if (length(point.colour) != 1)         stop('Error: Input point.colour should be a single integer')
  if (!(point.colour %in% colours()))    stop('Error: Input point.colour should be in \'colours()\'')

  ##############################################################################################################

  #Extract parameters
  OUT.SEARCH <- object$search
  OUT.RESULT <- object$result
  OUT.FLAG   <- object$flag.matrix
  KEYWORDS   <- OUT.SEARCH$Keyword
  CASE.SENS  <- OUT.SEARCH$case.sensitive
  WHOLE      <- OUT.SEARCH$whole.word
  TEXT       <- OUT.RESULT$Text
  m <- nrow(OUT.SEARCH)
  n <- nrow(OUT.RESULT)

  #Set keyword labels
  if (show.keywords) {
    for (i in 1:length(KEYWORDS)) {
      if (nchar(KEYWORDS[i]) > keyword.max) { KEYWORDS[i] <- paste0(substr(KEYWORDS[i], start = 1, stop = keyword.max), '...') } }
      KEYWORD_LABELS <- KEYWORDS
      for (i in 1:length(KEYWORDS)) {
        if (!WHOLE[i])     { KEYWORD_LABELS[i] <- paste0('*', KEYWORD_LABELS[i], '*') }
        if (!CASE.SENS[i]) { KEYWORD_LABELS[i] <- paste0('[', KEYWORD_LABELS[i], ']') } }
    } else {
    KEYWORD_LABELS <- rownames(OUT.SEARCH) }

  #Set text labels
  if (show.text) {
    for (i in 1:length(TEXT)) {
      if (nchar(TEXT[i]) > text.max) { TEXT[i] <- paste0(substr(TEXT[i], start = 1, stop = text.max), '...') } }
    TEXT_LABELS <- TEXT
  } else {
    TEXT_LABELS <- rownames(OUT.RESULT) }

  #Set subtitle and caption
  if (m == 1) {
    if (n == 1) {
      SUBTITLE <- paste0('Search of 1 keyword in 1 text')
    } else {
      SUBTITLE <- paste0('Search of 1 keyword in ', n, ' texts') }
  } else {
    if (n == 1) {
      SUBTITLE <- paste0('Search of ', m, ' keywords in 1 text')
    } else {
      SUBTITLE <- paste0('Search of ', m, ' keywords in ', n, ' texts') } }
  if (sum(!CASE.SENS) > 0) {
    CAPTION <- '[Bracketed keywords are not case sensitive]'
  } else {
    CAPTION <- NULL }

  #Set theme
  THEME <- ggplot2::theme(plot.title    = ggplot2::element_text(hjust = 0.5, size = 14, face = 'bold',
                                                                margin = ggplot2::margin(t = 0, r = 0, b = 10, l = 0)),
                          plot.subtitle = ggplot2::element_text(hjust = 0.5,
                                                                margin = ggplot2::margin(t = 0, r = 0, b = 6, l = 0)),
                          axis.title.x  = ggplot2::element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                          axis.title.y  = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
                          axis.text.x   = ggplot2::element_text(face = 'bold'),
                          plot.caption  = ggplot2::element_text(size = 8, hjust = 1,
                                                                margin = ggplot2::margin(t = 20, r = 0, b = 0, l = 0)),
                          axis.ticks    = ggplot2::element_blank(),
                          panel.background = ggplot2::element_rect(fill = 'white'),
                          panel.grid.major   = ggplot2::element_line(color = 'gray50', linetype = 'dotted', size = 0.5),
                          panel.grid.major.x = ggplot2::element_blank())

  #Generate plot data
  PLOTDATA <- data.frame(Text     = rep(TEXT_LABELS, each = m),
                         Keyword  = rep(KEYWORD_LABELS, times = n),
                         Found    = c(t(OUT.FLAG)))

  #Generate plot
  FIGURE  <- ggplot2::ggplot(ggplot2::aes(x = Keyword, y = Text, alpha = Found), data = PLOTDATA) +
    ggplot2::geom_point(colour = point.colour, size = point.size, shape = point.shape) +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::scale_alpha_manual(values = c(0, 1), guide = 'none') +
    THEME +
    ggplot2::labs(title = 'Keyword Search') +
    { if (!is.null(CAPTION)) { ggplot2::labs(caption = CAPTION) } } +
    { if (subtitle) { ggplot2::labs(subtitle = SUBTITLE) } } +
    ggplot2::xlab(NULL) + ggplot2::ylab(NULL)

  #Return plot
  FIGURE }
