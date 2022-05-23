#' Structure of Available Datasets
#'
#' \code{datasets.str} returns the structure of available datasets
#'
#' Datasets are often available in packages loaded into \code{R} and it is useful to know the structure of these datasets.  This function
#' shows the user the strucure of all available datasets in a specified package or packages.  (If the user does not specify a \code{package})
#' then the function searches over all available packages.
#'
#' @param package The package/packages containing the datasets of interest
#' @return A data frame listing available data sets, invisibly
#'
#' @examples
#' datasets.str("datasets")

datasets.str <- function(package = NULL) {

  if (missing(package)) { package <- .packages(all.available = TRUE) }

  #Create data frame listing available datasets
  ENV                <- new.env()
  DATASETS           <- as.data.frame(utils::data(package = package, envir = ENV)$results, stringsAsFactors = FALSE)[, c(1, 3, 4)]
  DATASETS$dataset   <- sub(".*[(]", "", sub("[)]$", "", DATASETS$Item))
  DATASETS$Item      <- sub(" .*", "", DATASETS$Item)
  DATASETS           <- DATASETS[order(DATASETS$Package, tolower(DATASETS$Item)),]
  DATASETS$class     <- ''

  #Print structure
  for (i in 1:nrow(DATASETS)) {
    DESC <- DATASETS[i, 1:4]
    data(list = DESC$dataset, package = DESC$Package, envir = ENV)
    DATA <- get(DESC$Item, envir = ENV)
    DATASETS$class[i] <- paste0(class(DATA), collapse = ', ')
    message("\n", DESC$Package, "  |  ", DESC$Item, "  |  ", DESC$Title)
    str(DATA, max.level = 1) }

  invisible(DATASETS) }
