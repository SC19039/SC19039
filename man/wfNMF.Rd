% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NMF.R
\name{wfNMF}
\alias{wfNMF}
\title{Standard non-negative matrix factorization algorithm.}
\usage{
wfNMF(V, r, maxiter)
}
\arguments{
\item{V}{the data}

\item{r}{The number of columns in W, or the number of rows in H.}

\item{maxiter}{Number of cycles}
}
\value{
Generated H matrix
}
\description{
Apply Non-negative matrix factorization to a picture.
}
\examples{
\dontrun{
data(V)
H1<-wfNMF(V,30,500)
}
}
