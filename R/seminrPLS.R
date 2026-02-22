#' PLS Structural Equation Modeling
#'
#' Performs Partial Least Squares Structural Equation Modeling (PLS-SEM)
#' using the seminr package.
#'
#' @param data a data frame
#' @param constructs a list of construct definitions. Each element should be a
#'   list with fields: \code{name} (character), \code{type} (one of
#'   \code{"reflective"}, \code{"composite_a"}, \code{"composite_b"}), and
#'   \code{indicators} (character vector of variable names).
#' @param paths a list of path definitions. Each element should be a list with
#'   fields: \code{from} (character, construct name) and \code{to} (character,
#'   construct name).
#' @param innerWeights inner weighting scheme: \code{"path_weighting"}
#'   (default) or \code{"factor_weighting"}.
#' @param missing missing data handling: \code{"mean_replacement"} (default)
#'   or \code{"listwise_deletion"}.
#' @param bootstrap logical; if \code{TRUE}, bootstrap the model for inference.
#' @param nboot integer; number of bootstrap samples (default 1000).
#' @param seed integer; random seed for reproducibility.
#'
#' @return A results object of class \code{SeminrPLS}.
#'
#' @examples
#' \dontrun{
#' # Using the mobi dataset from seminr
#' library(seminr)
#' data(mobi)
#'
#' result <- SeminrPLS(
#'   data = mobi,
#'   constructs = list(
#'     list(name = "Image",
#'          type = "reflective",
#'          indicators = c("IMAG1", "IMAG2", "IMAG3", "IMAG4", "IMAG5")),
#'     list(name = "Expectation",
#'          type = "reflective",
#'          indicators = c("CUEX1", "CUEX2", "CUEX3")),
#'     list(name = "Quality",
#'          type = "reflective",
#'          indicators = c("PERQ1", "PERQ2", "PERQ3", "PERQ4", "PERQ5", "PERQ6", "PERQ7")),
#'     list(name = "Value",
#'          type = "reflective",
#'          indicators = c("PERV1", "PERV2")),
#'     list(name = "Satisfaction",
#'          type = "reflective",
#'          indicators = c("CUSA1", "CUSA2", "CUSA3")),
#'     list(name = "Complaints",
#'          type = "reflective",
#'          indicators = c("CUSCO")),
#'     list(name = "Loyalty",
#'          type = "reflective",
#'          indicators = c("CUSL1", "CUSL2", "CUSL3"))
#'   ),
#'   paths = list(
#'     list(from = "Image",        to = "Expectation"),
#'     list(from = "Image",        to = "Satisfaction"),
#'     list(from = "Image",        to = "Loyalty"),
#'     list(from = "Expectation",  to = "Quality"),
#'     list(from = "Expectation",  to = "Value"),
#'     list(from = "Expectation",  to = "Satisfaction"),
#'     list(from = "Quality",      to = "Value"),
#'     list(from = "Quality",      to = "Satisfaction"),
#'     list(from = "Value",        to = "Satisfaction"),
#'     list(from = "Satisfaction", to = "Complaints"),
#'     list(from = "Satisfaction", to = "Loyalty"),
#'     list(from = "Complaints",   to = "Loyalty")
#'   ),
#'   bootstrap = TRUE,
#'   nboot = 500
#' )
#' }
#'
#' @export
SeminrPLS <- function(
        data,
        constructs = list(),
        paths = list(),
        innerWeights = "path_weighting",
        missing = "mean_replacement",
        bootstrap = FALSE,
        nboot = 1000L,
        seed = 123L) {

    if (base::missing(data))
        data <- jmvcore::marshalData(parent.frame())

    opts <- SeminrPLSOptions$new(
        constructs   = constructs,
        paths        = paths,
        innerWeights = innerWeights,
        missing      = missing,
        bootstrap    = bootstrap,
        nboot        = nboot,
        seed         = seed)

    analysis <- SeminrPLSClass$new(
        options    = opts,
        data       = data,
        datasetId  = "",
        analysisId = "",
        revision   = 0)

    analysis$run()
    analysis$results
}
