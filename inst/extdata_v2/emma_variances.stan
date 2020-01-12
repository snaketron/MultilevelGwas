data {
  int<lower=1> S; // number of SNPs
  int<lower=1> K; // number of strains
  int<lower=1> N; // number of observations
  matrix[N,S] X;  // design matrix
  matrix[N,K] Z;  // random effects design matrix
  vector[N] Y;    // traits
  matrix[K,K] A;  // kinship (relationship) matrix
}


transformed data {
  matrix[K,K] LA;
  LA = cholesky_decompose(A);
}


parameters {
  vector [K] a_decompose; // breeding values
  real b0;               // intercept
  real <lower=0> sigma_G; // genetic standard deviation
  real <lower=0> sigma_R; // residual standard deviation
}


model {
  vector[N] mu;
  vector[K] a;
  a_decompose ~ normal(0, 1);
  a = sigma_G * (LA * a_decompose);
  mu = b0 + Z * a;
  Y ~ normal(mu, sigma_R);
  b0 ~ normal(0, 100);
  sigma_G ~ student_t(4, 0, 1);
  sigma_R ~ student_t(4, 0, 1);
}


generated quantities{
  real sigma_U;
  real sigma_E;
  real h2;
  matrix[K,K] CD;
  sigma_U = sigma_G * sigma_G; // genetic variance
  sigma_E = sigma_R * sigma_R; // residual variance
  h2 = sigma_U / (sigma_E + sigma_U);
  CD = cholesky_decompose(A);
}

