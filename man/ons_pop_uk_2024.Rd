\name{ons_pop_uk_2024}
\alias{ons_pop_uk_2024}
\docType{data}
\title{UK mid-2024 population estimates by single year of age}
\description{
Mid-2024 population estimates for the United Kingdom (persons),
classified by single year of age from 0 to 89 plus an open "90+"
bin. Extracted from the \code{MYE2 - Persons} sheet of the ONS
mid-year population estimates, row for the United Kingdom (code
\code{K02000001}).

\strong{Cross-release rounding caveat.} The "90+" count in this
dataset (625,236) is the published aggregate from the \code{MYE2}
table. The sum of the single-year and 105+ counts in the companion
dataset \code{\link{ons_centenarians_uk_2024}} for the same year is
625,250, differing by 14 persons (about 0.002\%). This is not a
bug: the two ONS releases are rounded independently and use slightly
different bases. Reproductions of the open-tail analysis
(\code{inst/paper/03_open_tail.R}) should mention the discrepancy in
a footnote where the "truth" for ages 90+ is described, because
anyone re-running the figure will notice the small mismatch.
}
\usage{data(ons_pop_uk_2024)}
\format{
A data frame with 91 rows and 4 columns:
\describe{
  \item{age_lower}{Lower limit of the age band (years), integer.}
  \item{age_upper}{Upper limit of the age band (years), integer. For
    the open "90+" bin this is the finite cap 100.}
  \item{age_label}{Display label ("0", "1", ..., "89", "90+").}
  \item{count}{Observed count in the band. Total counts equal
    69,281,437.}
}
}
\source{
Office for National Statistics. \emph{Estimates of the population
for the UK, England, Wales, Scotland, and Northern Ireland},
mid-2024 edition. Sheet \code{MYE2 - Persons}. Released under the
Open Government Licence v3.0.
\url{https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland}
}
\seealso{
\code{\link{ons_centenarians_uk_2024}} for the complementary
single-year breakdown at ages 90 and over.
}
\examples{
data(ons_pop_uk_2024)
head(ons_pop_uk_2024)
tail(ons_pop_uk_2024)    # the "90+" open bin in the last row
}
\keyword{datasets}
