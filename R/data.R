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


# =============================================================================
# ONS datasets used by the paper reproduction scripts in inst/paper/
# -----------------------------------------------------------------------------
# These four data frames share a common schema:
#
#   age_lower  integer  lower limit of the age band (years).
#   age_upper  integer  upper limit of the age band (years). For the
#                       single-year-of-age rows this is age_lower + 1.
#                       For the open terminal bin (e.g. "105+") the
#                       upper limit is a finite cap (110 for the deaths
#                       and COVID datasets; 100 for the population
#                       dataset).
#   age_label  character display label ("0", "1", ..., "104", "105+").
#   count      integer  observed count in the band.
#
# They are extracted from public ONS releases (see each item's @source)
# and are kept deliberately small (a few rows times four columns) so
# that the paper's reproduction scripts run without needing the original
# spreadsheets. The R script that performs the extraction lives in
# data-raw/make-paper-data.R; users only need to run it when they want
# to refresh the datasets against a newer ONS release.
# =============================================================================

#' Deaths registered in England and Wales, 2022, by single year of age
#'
#' Counts of deaths registered in England and Wales during the calendar
#' year 2022, classified by single year of age (0, 1, ..., 104) plus an
#' open "105+" bin. Persons (males + females summed). Sourced from sheets
#' 5 and 6 of the ONS Deaths Registered in England and Wales (Series DR)
#' 2022 release.
#'
#' @format A data frame with 106 rows and 4 columns
#'   (\code{age_lower}, \code{age_upper}, \code{age_label}, \code{count}).
#'   Total counts equal \eqn{577{,}160}.
#' @source Office for National Statistics. \emph{Deaths registered in
#'   England and Wales (series DR), 2022}. Sheets 5 (males) and 6
#'   (females), row for 2022.
#'   \url{https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/deathsregisteredinenglandandwalesseriesdrreferencetables}
#' @examples
#' data(ons_deaths_ew_2022)
#' head(ons_deaths_ew_2022)
#' sum(ons_deaths_ew_2022$count)   # 577160
"ons_deaths_ew_2022"

#' UK mid-2024 population estimates by single year of age
#'
#' Mid-2024 population estimates for the United Kingdom (persons,
#' i.e. males and females combined), classified by single year of age
#' from 0 to 89 plus an open "90+" bin. Extracted from the
#' \code{MYE2 - Persons} sheet of the ONS mid-year population
#' estimates, row for the United Kingdom (code \code{K02000001}).
#'
#' \strong{Cross-release rounding caveat.} The "90+" count in this
#' dataset (\eqn{625{,}236}) is the published aggregate from the
#' \code{MYE2} table. The sum of the single-year and 105+ counts in
#' the companion dataset \code{\link{ons_centenarians_uk_2024}} for
#' the same year is \eqn{625{,}250}, differing by 14 persons (about
#' 0.002\%). This is not a bug: the two ONS releases are rounded
#' independently and use slightly different bases. Reproductions of
#' the open-tail analysis (\code{inst/paper/03_open_tail.R}) should
#' mention the discrepancy in a footnote where the "truth" for ages
#' 90+ is described, because anyone re-running the figure will notice
#' the small mismatch.
#'
#' @format A data frame with 91 rows and 4 columns
#'   (\code{age_lower}, \code{age_upper}, \code{age_label}, \code{count}).
#'   Total counts equal \eqn{69{,}281{,}437}.
#' @source Office for National Statistics. \emph{Estimates of the
#'   population for the UK, England, Wales, Scotland, and Northern
#'   Ireland}, mid-2024 edition. Sheet \code{MYE2 - Persons}.
#'   Released under the Open Government Licence v3.0.
#'   \url{https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland}
#' @seealso \code{\link{ons_centenarians_uk_2024}} for the
#'   complementary single-year breakdown at ages 90 and over.
#' @examples
#' data(ons_pop_uk_2024)
#' head(ons_pop_uk_2024)
#' tail(ons_pop_uk_2024)    # the "90+" open bin in the last row
"ons_pop_uk_2024"

#' UK 2024 population at very old ages
#'
#' Mid-2024 population estimates for the United Kingdom (persons),
#' classified by single year of age from 90 to 104 plus an open
#' "105+" bin. Provided by ONS as a separate release for the very old
#' because the main \code{MYE2} table closes the age distribution at
#' "90+".
#'
#' The 105+ count (\eqn{610}) is the input the open-tail analysis
#' relies on (see \code{inst/paper/03_open_tail.R}).
#'
#' \strong{Cross-release rounding caveat.} The total of this dataset
#' (single-year counts at ages 90 to 104 plus the 105+ count) is
#' \eqn{625{,}250}, which differs from the \code{MYE2} "90+"
#' aggregate of \eqn{625{,}236} in \code{\link{ons_pop_uk_2024}} by
#' 14 persons (about 0.002\%). This is not a bug: the two ONS
#' releases are rounded independently. Manuscripts reproducing the
#' open-tail figure should mention the discrepancy in a footnote
#' where the "truth" at ages 90+ is described.
#'
#' @format A data frame with 16 rows and 4 columns
#'   (\code{age_lower}, \code{age_upper}, \code{age_label}, \code{count}).
#' @source Office for National Statistics. \emph{Mid-year population
#'   estimates of the very old, including centenarians: UK}, edition
#'   covering 2002 to 2024. Row for Persons, 2024. Released under the
#'   Open Government Licence v3.0.
#'   \url{https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/ageing/datasets/midyearpopulationestimatesoftheveryoldincludingcentenariansunitedkingdom}
#' @seealso \code{\link{ons_pop_uk_2024}} for the complementary
#'   single-year breakdown at ages 0 to 89.
#' @examples
#' data(ons_centenarians_uk_2024)
#' ons_centenarians_uk_2024
"ons_centenarians_uk_2024"

#' Deaths due to COVID-19, England and Wales, March 2020 to March 2023
#'
#' Counts of deaths where COVID-19 was the underlying cause ("due to
#' COVID-19"), aggregated across all 37 monthly cohorts from March 2020
#' to March 2023 inclusive, England and Wales, by single year of age
#' from 0 to 104 plus an open "105+" bin. Persons (males and females
#' combined).
#'
#' @format A data frame with 106 rows and 4 columns
#'   (\code{age_lower}, \code{age_upper}, \code{age_label}, \code{count}).
#'   Total counts equal \eqn{168{,}418}.
#' @source Office for National Statistics. \emph{Single year of age and
#'   average age of death of people whose death was due to or involved
#'   coronavirus (COVID-19), England and Wales}. Sheet 1, rows whose
#'   Cause is "Due to COVID-19", aggregated across the 37 cohorts from
#'   March 2020 to March 2023 inclusive. Released under the Open
#'   Government Licence v3.0.
#'   \url{https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/singleyearofageandaverageageofdeathofpeoplewhosedeathwasduetoorinvolvedcovid19}
#' @examples
#' data(ons_covid_due_ew)
#' head(ons_covid_due_ew)
#' sum(ons_covid_due_ew$count)   # 168418
"ons_covid_due_ew"
