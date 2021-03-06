% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NMF.R
\name{SNMF}
\alias{SNMF}
\title{Sparse non-negative matrix factorization using BeSS.}
\usage{
SNMF(V, r, maxiter, kong)
}
\arguments{
\item{V}{the data}

\item{r}{the number of columns in W, or the number of rows in H.}

\item{maxiter}{the number of cycles}

\item{kong}{t he number of non-zero elements of the matrix H}
}
\value{
Generated H matrix
}
\description{
Apply Sparse nonnegative matrix factorization to a picture.
}
\examples{
\dontrun{
data(V)
H2<-SNMF(V,30,500,800)
}
}
