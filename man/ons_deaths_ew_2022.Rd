\name{ons_deaths_ew_2022}
\alias{ons_deaths_ew_2022}
\docType{data}
\title{Deaths registered in England and Wales, 2022, by single year of age}
\description{
Counts of deaths registered in England and Wales during the calendar
year 2022, classified by single year of age (0, 1, ..., 104) plus an
open "105+" bin. Persons (males and females summed). Used by the
paper's ONS benchmark in \code{inst/paper/02_ons_benchmark.R}.
}
\usage{data(ons_deaths_ew_2022)}
\format{
A data frame with 106 rows and 4 columns:
\describe{
  \item{age_lower}{Lower limit of the age band (years), integer.}
  \item{age_upper}{Upper limit of the age band (years), integer. For
    the open "105+" bin this is the finite cap 110.}
  \item{age_label}{Display label ("0", "1", ..., "104", "105+").}
  \item{count}{Observed count in the band. Total counts equal 577,160.}
}
}
\source{
Office for National Statistics. \emph{Deaths registered in England and
Wales (series DR), 2022}. Sheets 5 (males) and 6 (females), row for
2022.
\url{https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/deathsregisteredinenglandandwalesseriesdrreferencetables}
}
\examples{
data(ons_deaths_ew_2022)
head(ons_deaths_ew_2022)
sum(ons_deaths_ew_2022$count)   # 577160
}
\keyword{datasets}
