data {
  int<lower=1> S; // number of SNPs
  int<lower=1> K; // number of strains
  int<lower=1> N; // number of observations
  matrix[N,S] X;  // design matrix
  matrix[N,K] Z;  // random effects design matrix
  vector[N] Y;    // traits
  matrix[K,K] A;  // kinship (relationship) matrix
}



parameters {
  real b0;               // intercept
  vector[S] b;           // SNP effects
  real<lower=0> sigma_R; // residual standard deviation
}


model {
    for(n in 1:N) {
      Y[n] ~ normal(b0 + b .* to_vector(X[n, ]), sigma_R);
    }
    b0 ~ normal(0, 100);
    to_vector(b) ~ normal(0, 1);
    sigma_R ~ student_t(4, 0, 1);
}
