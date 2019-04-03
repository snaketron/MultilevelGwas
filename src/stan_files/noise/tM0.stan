data {
  int N; // number of all entries
  int Ntq; // number of continuous traits
  int Ntd; // number of dichotomous traits
  int Ns; // number of all SNPs
  real Yq[N, Ntq] ; // number of hits response
  int Yd[N, Ntd] ; // number of hits response
  vector [Ns] X [N]; // index of all individuals
}

parameters {
  vector [Ntq+Ntd] alpha;
  vector <lower = 0> [Ntq] sigma;
  matrix [Ntq+Ntd, 1] mu_beta;
  vector <lower = 0> [Ntq+Ntd] sigma_beta;
  vector<lower=0> [Ntq+Ntd] u; // helper
  cholesky_factor_corr [Ntq+Ntd] L_rho;
  matrix [Ntq+Ntd, Ns] z;
}

transformed parameters {
  matrix [Ntq+Ntd, Ns] beta;
  for(s in 1:Ns) {
     // multi student
     beta[, s] = mu_beta[, 1] + sqrt(sigma_beta ./ u) .* z[, s];
  }
}

model {
  for(i in 1:N) {
    if(Ntq > 0) {
      for(t in 1:Ntq) {
        Yq[i,t] ~ normal(alpha[t] + rows_dot_product(X[i], to_vector(beta[t, ])), sigma[t]);
      }
    }
    if(Ntd > 0) {
      for(d in 1:Ntd) {
        Yd[i,d] ~ bernoulli_logit(alpha[d+Ntq] + rows_dot_product(X[i], to_vector(beta[d+Ntq, ])));
      }
    }
  }

  alpha ~ student_t(1, 0, 100);
  mu_beta[, 1] ~ student_t(1, 0, 10);
  sigma ~ cauchy(0, 5);
  sigma_beta ~ cauchy(0, 5);
  u ~ chi_square(sigma_beta);

  for(s in 1:Ns) {
    z[, s] ~ normal(0, 1);
  }
}
