% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{.is_angular}
\alias{.is_angular}
\title{Check Whether an Object is Angular}
\usage{
.is_angular(x)
}
\description{
An internal helper function that checks whether an object \code{x} is marked as
using angular frequency. It looks for the \code{.is_angular} attribute.
}
\details{
This function is not exported and is only used internally by \code{to_angf()}
and \code{to_cycf()}.
}
\keyword{internal}
