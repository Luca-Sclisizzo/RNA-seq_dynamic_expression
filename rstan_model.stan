//
// This Stan program defines a hiherarchical model, with a
// vector of values 'y' modeled as Negative-Binomial distributed
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=1> N;            // number of observations
  int<lower=1> K;            // spline basis dimension
  int<lower=1> S;            // number of subjects

  vector[N] age;

  array[N] int<lower=1> subject;  // subject ID (1..S)

  matrix[N, K] B;            // spline basis

  array[N] int<lower=0> y;   // counts (expression)
}
// The parameters accepted by the model
parameters {
  real alpha; // intercept of the spline
  vector[K] beta_spline;
  vector[S] z_subject; // standard normal (non-centered) - to solve the funnel structure

  real<lower=0> sigma_u;
  real<lower=1e-6> phi;  // NB dispersion, a little offset to prevent exact zero
}
// Non centered parametrization to solve funnel structure
transformed parameters {
  vector[S] u_subject = sigma_u * z_subject;  // building the subject-wise distribution using z_subject
}
// The model to be estimated. We model the output
// 'y' to be NB distributed
model {
  // priors
  beta_spline ~ normal(0, 1);
  z_subject ~ normal(0, 1);     // prior on raw parameter
  sigma_u ~ normal(0, 1); // subject-wise dispersion parameter (the shrinkage)
  phi ~ exponential(1); // NB dispersion parameter
  alpha ~ normal(0, 2); // spline's intercept
  
  vector[N] eta = alpha + B * beta_spline + u_subject[subject]; // pointwise eta to use for the generated_quantities{} chunk
  y ~ neg_binomial_2_log(eta, phi);
}
// Generate quantities: point likelihood for loo estimate after the fit
generated quantities {
  vector[N] log_lik;
  vector[N] eta = alpha + B * beta_spline + u_subject[subject];
  
  for (n in 1:N)
    log_lik[n] = neg_binomial_2_log_lpmf(y[n] | eta[n], phi);
}
