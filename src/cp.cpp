#include <Rcpp.h>
using namespace Rcpp;

//' @title A markov sampler using Rcpp
//' @descriptionImplement a random walk Metropolis sampler for generating the standard Laplace distribution
//' @param sigma the variance
//' @param xx the Initial value
//' @param N Number of random numbers
//' @return Returns random sequence and rejection probability
//' @examples
//' \dontrun{
//' rp <-(2,6,1000)
//' }
//' @export
// [[Rcpp::export]]
NumericVector rp(double sigma,double xx,int N){
#include<Rcpp.h>
  NumericVector x(N+1);
  x[0]=xx;
  int k=0;
  for(int i=1;i<N;++i){
    srand(i);
    double u=rand()/(RAND_MAX+1.0);
    double y=as<double>(rnorm(1, x[i-1], sigma));
    if(u<=(0.5*exp(-abs(y)))/(0.5*exp(-abs(x[i-1]))))
    {x[i]=y;}
    else{
      x[i]=x[i-1];
      k=k+1;
    }
  }
  x[N]=k;
  return x;
}