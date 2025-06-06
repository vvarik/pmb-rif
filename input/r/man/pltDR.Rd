% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pltDR.R
\name{pltDR}
\alias{pltDR}
\title{Plot dose-response for monotherapies}
\usage{
pltDR(
  formula,
  dat,
  col = "red3",
  by = cond,
  type = average,
  ci95 = T,
  zero = 10,
  ...
)
}
\arguments{
\item{formula}{A formula for fit}

\item{dat}{Data in data.table}

\item{col}{The color of lines, points}

\item{by}{A grouping variable for fit}

\item{type}{A type of drc plot, one of 'average', 'all', 'none',
'confidence', 'bars', 'obs'.}
}
\description{
Plot dose-response for monotherapies
}
