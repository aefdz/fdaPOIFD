#' exampleData
#'
#' An illustrative Functional Gaussian processes with different partially observed patterns
#' with outliers and without outliers.
#'
#' @format A list with three data sets (functions by columns):
#' \describe{
#' \item{PoFDintervals}{Partially observed functional data in intervals}
#' \item{PoFDextremes}{Partially Observed functional data with missing intervals at the extremes}
#' \item{PoFDextremes_outliers}{Same as above but including two magnitude and shape outliers}
#' \item{PoFDintervals_list}{Partially observed data without a common domain. Each function is one element of the list, containing the evaluation points (x) and the evaluated function (y).}
#' }
#'
#' @references Elías, Antonio, Jiménez, Raúl, Paganoni, Anna M. and Sangalli, Laura M. (2020). Integrated Depths for Partially Observed Functional Data.
#'
#' @examples
#' data(exampleData)
#' plotPOFD(exampleData$PoFDintervals)
#'
"exampleData"
