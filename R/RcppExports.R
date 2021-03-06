# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title A markov sampler using Rcpp
#' @descriptionImplement a random walk Metropolis sampler for generating the standard Laplace distribution
#' @param sigma the variance
#' @param xx the Initial value
#' @param N Number of random numbers
#' @return Returns random sequence and rejection probability
#' @examples
#' \dontrun{
#' rp <-(2,6,1000)
#' }
#' @export
rp <- function(sigma, xx, N) {
    .Call(`_SC19039_rp`, sigma, xx, N)
}

