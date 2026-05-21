\name{ons_covid_due_ew}
\alias{ons_covid_due_ew}
\docType{data}
\title{Deaths due to COVID-19, England and Wales, March 2020 to March 2023}
\description{
Counts of deaths where COVID-19 was the underlying cause ("due to
COVID-19"), aggregated across all 37 monthly cohorts from March 2020
to March 2023 inclusive, England and Wales, by single year of age
from 0 to 104 plus an open "105+" bin. Persons (males and females
combined).
}
\usage{data(ons_covid_due_ew)}
\format{
A data frame with 106 rows and 4 columns:
\describe{
  \item{age_lower}{Lower limit of the age band (years), integer.}
  \item{age_upper}{Upper limit of the age band (years), integer. For
    the open "105+" bin this is the finite cap 110.}
  \item{age_label}{Display label ("0", "1", ..., "104", "105+").}
  \item{count}{Observed count in the band. Total counts equal 168,418.}
}
}
\source{
Office for National Statistics. \emph{Single year of age and average
age of death of people whose death was due to or involved
coronavirus (COVID-19), England and Wales}. Sheet 1, rows whose
Cause is "Due to COVID-19", aggregated across the 37 cohorts from
March 2020 to March 2023 inclusive. Released under the Open
Government Licence v3.0.
\url{https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/datasets/singleyearofageandaverageageofdeathofpeoplewhosedeathwasduetoorinvolvedcovid19}
}
\examples{
data(ons_covid_due_ew)
head(ons_covid_due_ew)
sum(ons_covid_due_ew$count)   # 168418
}
\keyword{datasets}
