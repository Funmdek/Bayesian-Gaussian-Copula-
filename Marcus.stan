functions{
  real binormal_cdf(real z1, real z2, real rho) {
  if (z1 != 0 || z2 != 0) {
    real denom = fabs(rho) < 1.0 ? sqrt((1 + rho) * (1 - rho)) : not_a_number();
    real a1 = (z2 / z1 - rho) / denom;
    real a2 = (z1 / z2 - rho) / denom;
    real product = z1 * z2;
    real delta = product < 0 || (product == 0 && (z1 + z2) < 0);
    return 0.5 * (Phi(z1) + Phi(z2) - delta) - owens_t(z1, a1) - owens_t(z2, a2);
  }
  return log((0.25 + asin(rho) / (2 * pi())));
}
  real my_custom_likelihood_lpdf(vector y, matrix x, matrix z, real alpha, real beta10, real beta20, int N, vector rho, vector beta1, vector beta2){
    int t = size(y);
    // int f = t/2;
    int f = t/2;
    vector[f] y1;
    vector[f] y2;
    for (i in 1:f){
      y1[i] = y[i];
      int ff = f+i;
      y2[i] = y[ff]; 
    }
    real loglike = 0;
    real loglike1 = 0;
    real loglike2 = 0;
    real loglike3 = 0;
    real loglike4 = 0 ;
      for (i in 1:N) {
        loglike1 = loglike+((1-y1[i])*(1-y2[i]))* (binormal_cdf(Phi_approx(-x[i]*beta1)-beta10, Phi_approx(-z[i]*beta2)-beta20, rho[i]));
        loglike2 = loglike+ ((y1[i])*(1-y2[i]))* ((binormal_cdf(1,(Phi_approx(-x[i]*beta1)-beta10), rho[i])) - (binormal_cdf((Phi_approx(-x[i]*beta1)-beta10), (Phi_approx(-z[i]*beta2)-beta20), rho[i])));
        loglike3 = loglike+((1-y1[i])*(y2[i]))* ((binormal_cdf(((Phi_approx(-x[i]*beta1)-beta10 - alpha)), 1, rho[i]) - (binormal_cdf((Phi_approx(-x[i]*beta1)-beta10), (Phi_approx(-z[i]*beta2)-beta20), rho[i]))));
        loglike4 = loglike+((y1[i])*(y2[i]))* (1-binormal_cdf(Phi_approx((-x[i]*beta1)-beta10-alpha), 1, rho[i]) - binormal_cdf(1,(Phi_approx(-x[i]*beta1)-beta10), rho[i]) + (binormal_cdf((Phi_approx(-x[i]*beta1)-beta10), Phi_approx((-z[i]*beta2)-beta20),rho[i])));
      }
      loglike = loglike1 + loglike2 + loglike3 + loglike4;
      return loglike;
    }
}
data {
  int <lower=0> N; // number of data items
  //int t = 2*N;
  vector[2 * N] y;       // Observed data (dependent variable)
  matrix[N,3] x;
  matrix[N,3] z;
  //int<lower=0,upper=1> y[N];
}
parameters { 
  vector[3] beta1;
  vector[3] beta2;
  real beta10;
  real beta20;
  vector[3] beth;
  real betc;
  //real<lower=0> sigma1;
  //real<lower=0> sigma2;
  real alpha;
  //real rho;
  //matrix[2,2] Sigma;
  //real nu;
  }
transformed parameters{
  vector[N] rho;
  vector[N] h;
  for (i in 1:N){
    h[i] = x[i]*beth+betc;
    rho[i] = tanh(h[i]);
  }
}
model {
  // Prior for the regression coefficients
  beta1 ~ normal (0, 100);
  beta2 ~ normal (0, 100); 
  beth ~ normal (0, 100);
  betc ~ normal(0, 100);
  //beta1 ~ multi_normal([0,0],[[100,0],[0,100]]);// mu is a vector and sigma is a matrix 
  //beta2 ~ multi_normal([0,0.5],[[100,4],[0,100]]);
  alpha ~ gamma(0.3, 0.1);
  //sigma1 ~ gamma(0.1, 0.6);
  //sigma2 ~ gamma(0.2, 0.1);
  beta10 ~ normal(0, 100);
  beta20 ~ normal(0, 100);
  //h ~ normal(0.3, 0.8);
  //nu ~ normal(0,1);
  //w ~ wishart_rng(nu, Sigma)
  // Custom likelihood using the 'my_custom_likelihood' distribution'
  //target += target + (y2, beta1, beta2, beta10, beta20,  z,  x, rho, alpha, N);
  y ~ my_custom_likelihood(x, z, alpha, beta10, beta20, N, rho, beta1, beta2);
  }
