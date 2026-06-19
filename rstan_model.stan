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
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> S;
  int<lower=1> M;
  int<lower=1> G;                          // numero di geni

  array[N] int<lower=1, upper=M> module;
  array[N] int<lower=1, upper=S> subject;
  array[N] real log_offset;
  array[N] int<lower=1, upper=G> gene;     // indice gene per ogni obs
  matrix[N, K] B;
  array[N] int<lower=0> y;
  
  // mappa gene -> modulo (per il prior gerarchico)
  array[G] int<lower=1, upper=M> gene_module;
}

parameters {
  real alpha;

  // Trend del modulo
  matrix[M, K] beta_module;
  vector[K] mu_beta;
  vector<lower=0>[K] sigma_beta;

  // Magnitudine gene-specifica (non-centered)
  vector[G] z_gene;
  vector<lower=0>[M] sigma_gene;   // varianza di magnitudine per modulo

  // Soggetti
  vector[S] z_subject;
  real<lower=0> sigma_u;

  real<lower=1e-6> phi;
}

transformed parameters {
  vector[S] u_subject = sigma_u * z_subject;

  // Effetto gene: non-centered, varianza dipende dal modulo
  vector[G] gamma_gene;
  for (g in 1:G)
    gamma_gene[g] = sigma_gene[gene_module[g]] * z_gene[g];
}

model {
  // Iperpriori modulo
  mu_beta    ~ normal(0, 1);
  sigma_beta ~ normal(0, 0.5);
  for (m in 1:M)
    beta_module[m] ~ normal(mu_beta, sigma_beta);

  // Prior magnitudine gene
  z_gene     ~ normal(0, 1);
  sigma_gene ~ normal(0, 0.5);    // shrinkage: geni dello stesso modulo

  // Soggetti
  z_subject ~ normal(0, 1);
  sigma_u   ~ normal(0, 1);

  // Globali
  alpha ~ normal(0, 2);
  phi   ~ exponential(1);
  
  // Identificabilità
  sum(z_gene)    ~ normal(0, 0.1 * sqrt(G));
  sum(z_subject) ~ normal(0, 0.1 * sqrt(S));
  
  // Likelihood
  vector[N] eta;
  for (n in 1:N)
    eta[n] = alpha
             + log_offset[n]
             + dot_product(B[n], beta_module[module[n]])
             + gamma_gene[gene[n]]       // within-gene variability in the module
             + u_subject[subject[n]];

  y ~ neg_binomial_2_log(eta, phi);
}

generated quantities {
  vector[N] log_lik;
  array[N] int y_rep;

  for (n in 1:N) {
    real eta_n = alpha
                 + log_offset[n]
                 + dot_product(B[n], beta_module[module[n]])
                 + gamma_gene[gene[n]]
                 + u_subject[subject[n]];
    
    log_lik[n] = neg_binomial_2_log_lpmf(y[n] | eta_n, phi);
    y_rep[n]   = neg_binomial_2_log_rng(eta_n, phi);          // y_rep is commented and will be used later to calculate 100 y_rep values for PPC
  }
}
