# =============================================================================
# Documentation for example datasets
# =============================================================================

#' Lead concentration in the blood of New York children, 1974
#'
#' Interval-censored measurements of lead concentration (in
#' \eqn{\mu}g/dl) in the blood of young Puerto Ricans aged 1-12 years
#' living in New York in 1974, as reported by Hasselblad et al. (1980)
#' and re-used as Example 6.1 by Lambert and Eilers (2009).
#'
#' The original instrument's resolution was limited; observations were
#' recorded only as falling into one of seven wide intervals. The last
#' interval, \emph{65+}, was bounded above by 80 in the analysis of
#' Lambert and Eilers (2009). Total \eqn{n = 139}.
#'
#' @format A data frame with 7 rows and 3 columns:
#' \describe{
#'   \item{lower}{Lower limit of the wide bin (\eqn{\mu}g/dl).}
#'   \item{upper}{Upper limit of the wide bin (\eqn{\mu}g/dl).}
#'   \item{count}{Observed frequency \eqn{m_j} in the bin.}
#' }
#'
#' @references
#' Hasselblad, V., Stead, A. G. and Galke, W. (1980). Analysis of
#' coarsely grouped data from the lognormal distribution.
#' \emph{Journal of the American Statistical Association}, 75,
#' 771--778.
#'
#' Lambert, P. and Eilers, P. H. C. (2009). Bayesian density
#' estimation from grouped continuous data. \emph{Computational
#' Statistics and Data Analysis}, 53(4), 1388--1399.
#'
#' @examples
#' data(bloodlead)
#' bloodlead
"bloodlead"

#' Tuberculosis deaths by age, The Netherlands, 1907
#'
#' Counts of deaths attributed to tuberculosis in The Netherlands in
#' 1907, classified by wide age bands. Used as Example 6.2 by Lambert
#' and Eilers (2009) to illustrate density estimation from heavily
#' grouped mortality data. Total deaths in the dataset: 9440.
#'
#' Yearly population numbers in single-year age intervals were
#' available, but TB deaths were only recorded in wide, irregular age
#' bands (a common feature of historical mortality data). An extra
#' interval [100, 120) with zero count is added at the upper tail to
#' force the estimated density to taper smoothly to zero, as suggested
#' by Lambert and Eilers (2009, Section 4).
#'
#' @format A data frame with 12 rows and 3 columns:
#' \describe{
#'   \item{lower}{Lower limit of the age band (years).}
#'   \item{upper}{Upper limit of the age band (years).}
#'   \item{count}{Number of deaths attributed to tuberculosis in the
#'     band.}
#' }
#'
#' @references
#' Lambert, P. and Eilers, P. H. C. (2009). Bayesian density
#' estimation from grouped continuous data. \emph{Computational
#' Statistics and Data Analysis}, 53(4), 1388--1399.
#'
#' @examples
#' data(tbdeaths1907)
#' tbdeaths1907
"tbdeaths1907"
