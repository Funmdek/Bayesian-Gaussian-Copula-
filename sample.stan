functions{
 real my_custom_likelihood_lpdf(vector y, real sigma,  vector beta, matrix x, int N){
  real loglike = 0;
  for (i in 1:N) {
  loglike = loglike+(-log(sigma) - 0.5 * (y[i] - (x[i,]*beta) )  / (sigma * sigma));
  }
  return loglike;
}
}
data {
  int <lower=0> N; // number of data items
  vector[N] y;       // Observed data (dependent variable)
  matrix[N,2] x;
  //int<lower=0,upper=1> y[N];
}
parameters { 
  vector[2] beta;
  real<lower=0> sigma;
  }

model {
  // Prior for the regression coefficients
  beta ~ multi_normal([0,0],[[100,0],[0,100]]);// mu is a vector and sigma is a matrix 
  sigma ~ gamma(0.1, 0.6);
  // Custom likelihood using the 'my_custom_likelihood' distribution
  y ~ my_custom_likelihood(sigma, beta, x, N);
}

