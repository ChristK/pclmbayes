\name{ons_centenarians_uk_2024}
\alias{ons_centenarians_uk_2024}
\docType{data}
\title{UK 2024 population at very old ages}
\description{
Mid-2024 population estimates for the United Kingdom (persons),
classified by single year of age from 90 to 104 plus an open "105+"
bin. ONS publishes this separately from the main \code{MYE2} table
because \code{MYE2} closes the age distribution at "90+". The 105+
count (610) is the input the open-tail analysis relies on (see
\code{inst/paper/03_open_tail.R}).

\strong{Cross-release rounding caveat.} The total of this dataset
(single-year counts at ages 90 to 104 plus the 105+ count) is
625,250, which differs from the \code{MYE2} "90+" aggregate of
625,236 in \code{\link{ons_pop_uk_2024}} by 14 persons (about
0.002\%). This is not a bug: the two ONS releases are rounded
independently. Manuscripts reproducing the open-tail figure should
mention the discrepancy in a footnote where the "truth" at ages 90+
is described.
}
\usage{data(ons_centenarians_uk_2024)}
\format{
A data frame with 16 rows and 4 columns:
\describe{
  \item{age_lower}{Lower limit of the age band (years), integer.}
  \item{age_upper}{Upper limit of the age band (years), integer. For
    the open "105+" bin this is the finite cap 110.}
  \item{age_label}{Display label ("90", "91", ..., "104", "105+").}
  \item{count}{Observed count in the band.}
}
}
\source{
Office for National Statistics. \emph{Mid-year population estimates
of the very old, including centenarians: UK}, edition covering 2002
to 2024. Row for Persons, 2024. Released under the Open Government
Licence v3.0.
\url{https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/ageing/datasets/midyearpopulationestimatesoftheveryoldincludingcentenariansunitedkingdom}
}
\seealso{
\code{\link{ons_pop_uk_2024}} for the complementary single-year
breakdown at ages 0 to 89.
}
\examples{
data(ons_centenarians_uk_2024)
ons_centenarians_uk_2024
}
\keyword{datasets}
