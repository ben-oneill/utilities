#' Structure of Available Datasets
#'
#' \code{datasets.str} returns the structure of available datasets
#'
#' Datasets are often available in packages loaded into \code{R} and it is useful to know the structure of these datasets.  This function
#' shows the user the strucure of all available datasets in a specified package or packages.  (If the user does not specify a \code{package})
#' then the function searches over all available packages.
#'
#' @usage \code{datasets.str}
#' @param package The package/packages containing the datasets of interest
#' @return Print output showing the structure of all datasets

datasets.str <- function(package = NULL) {

  if (missing(package)) { package <- .packages(all.available = TRUE) }

  #Create data frame listing available datasets
  ENV           <- new.env()
  DATASETS      <- as.data.frame(data(package = package, envir = ENV)$results, stringsAsFactors = FALSE)[, c(1, 3, 4)]
  DATASETS$Item <- sapply(strsplit(DATASETS$Item, split = ' ', fixed = TRUE), '[', 1)
  DATASETS      <- DATASETS[order(DATASETS$Package, tolower(DATASETS$Item)),]

  #Print structure
  for (i in 1:nrow(DATASETS)) {
    DESC <- DATASETS[i, ]
    data(list = DESC$Item, package = DESC$Package, envir = ENV)
    DATA <- get(DESC$Item, envir = ENV)
    message(DESC$Package[i], "  |  ", DESC$Item, "  |  ", DESC$Title)
    message(str(DATA)) } }
