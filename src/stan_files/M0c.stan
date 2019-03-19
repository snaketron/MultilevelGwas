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
  vector <lower = 0> [Ntq+Ntd] sd_beta;
  cholesky_factor_corr [Ntq+Ntd] L_rho;
  matrix [Ntq+Ntd, Ns] z;
}

transformed parameters {
  matrix [Ntq+Ntd, Ns] beta;
  for(s in 1:Ns) {
     // multi_normal
     beta[, s] = mu_beta[, 1] + diag_pre_multiply(sd_beta, L_rho)*z[, s];
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
  sd_beta ~ cauchy(0, 5);
  L_rho ~ lkj_corr_cholesky(2);

  for(s in 1:Ns) {
    z[, s] ~ normal(0, 1);
  }
}

generated quantities {
  corr_matrix[Ntq+Ntd] rho;
  matrix[Ntq+Ntd, Ntq+Ntd] SD_cov;
  rho = multiply_lower_tri_self_transpose(L_rho);
  SD_cov = quad_form_diag(rho, sd_beta);
}
