% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{.set_repr}
\alias{.set_repr}
\title{Set Representation Format}
\usage{
.set_repr(x, repr, .keep = "unused")
}
\arguments{
\item{x}{An R object.}

\item{repr}{The representation type (\code{"cplx"}, \code{"rect"}, \code{"polr"}).}

\item{.keep}{Controls preservation of data.}
}
\value{
The object in the specified representation format.
}
\description{
Converts an object into a specific representation format.
}
\keyword{internal}
