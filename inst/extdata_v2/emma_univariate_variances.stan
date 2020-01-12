data {
  int<lower=1> S; // number of SNPs
  int<lower=1> K; // number of strains
  int<lower=1> N; // number of observations
  matrix[N,S] X;  // design matrix
  matrix[N,K] Z;  // random effects design matrix
  vector[N] Y;    // traits
  matrix[K,K] A;  // kinship (relationship) matrix
  int Ks [N];     // strain ids of different mice
}


transformed data {
  matrix[K,K] LA;
  LA = cholesky_decompose(A);
}


parameters {
  vector[K] a_decompose; // breeding values
  real b0;               // intercept
  real<lower=0> sigma_G; // genetic standard deviation
  real<lower=0> sigma_R; // residual standard deviation

  vector [S] b; // snp effects
  real b_mu;
  real <lower=0> b_sigma;
  real <lower=0> b_nu;
  vector [S] b_z;
}


// transformed parameters {
//   vector [S] b; // snp effects
//   b = b_mu + b_sigma * b_z;
// }

model {
    // vector[N] mu;
    matrix [N,S] mus;
    vector[K] a;
    a_decompose ~ normal(0, 1);
    a = sigma_G * (LA * a_decompose);

    // for(s in 1:S) {
    //   mu = b0 + X[, s] * b[s] + Z * a;
    //   Y ~ normal(mu, sigma_R);
    // }

    for(n in 1:N) {
      mus[n,] = to_row_vector(b0+to_vector(X[n,]).*b+a[Ks[n]]);
    }
    for(s in 1:S) {
      Y ~ normal(mus[, s], sigma_R);
    }
    // to_vector(b) ~ normal(0, 5);
    b0 ~ normal(0, 100);
    sigma_G ~ cauchy(0, 1);
    sigma_R ~ cauchy(0, 1);
    b_mu ~ normal(0, 5);
    b_sigma ~ cauchy(0, 1);
    b_nu ~ gamma(2, 0.1);
    b_z ~ normal(0, 5);

    b ~ student_t(b_nu, b_mu, b_sigma);
}


generated quantities{
  real sigma_U;
  real sigma_E;
  real h2;
  sigma_U = sigma_G * sigma_G; // genetic variance
  sigma_E = sigma_R * sigma_R; // residual variance
  h2 = sigma_U / (sigma_E + sigma_U);
}


