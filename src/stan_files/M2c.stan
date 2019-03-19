data {
  int N; // number of all entries
  int Ntc; // number of traits
  int Ntd; // number of traits
  int Ns; // number of all SNPs
  int Nk; // number of all strains
  real Yc[N, Ntc] ; // number of hits response
  int Yd[N, Ntd] ; // number of hits response
  vector [Ns] X [N]; // index of all individuals
  int K[N]; //index to all strains
}

parameters {
  vector [Ntc+Ntd] alpha;
  matrix [Ntc+Ntd, 1] grand_mu_beta;
  vector <lower = 0> [Nk] sigma;
  real <lower = 0> sigma_trait;
  real <lower = 0> mean_trait;
  vector <lower = 0> [Ntc+Ntd] sigma_beta;
  vector <lower = 0> [Ntc+Ntd] grand_sigma_beta;
  matrix [Ns, Nk] z [Ntc+Ntd];
  matrix [Ntc+Ntd, Ns] grand_z;
  cholesky_factor_corr [Ntc+Ntd] L_rho;
}

transformed parameters {
  matrix [Ns, Nk] beta [Ntc+Ntd];
  matrix [Ntc+Ntd, Ns] mu_beta;

  for(s in 1:Ns) {
    // multi_normal
    mu_beta[, s] = grand_mu_beta[, 1] + diag_pre_multiply(grand_sigma_beta, L_rho)*grand_z[, s];
  }


  for(t in 1:(Ntc+Ntd)) {
    for(s in 1:Ns) {
      for(k in 1:Nk) {
        // normal
        beta[t][s, k] = mu_beta[t, s] + z[t][s, k]*sigma_beta[t];
      }
    }
  }
}

model {
  for(i in 1:N) {
    if(Ntc > 0) {
      for(t in 1:Ntc) {
        Yc[i,t] ~ normal(alpha[t] + rows_dot_product(X[i], to_vector(beta[t][, K[i]])), sigma[t]);
      }
    }
    if(Ntd > 0) {
      for(d in 1:Ntd) {
        Yd[i,d] ~ bernoulli_logit(alpha[d+Ntc] + rows_dot_product(X[i], to_vector(beta[d+Ntc][, K[i]])));
      }
    }
  }

  alpha ~ student_t(1, 0, 100);
  to_vector(grand_mu_beta) ~ student_t(1, 0, 10);

  sigma_beta ~ cauchy(0, 5);
  grand_sigma_beta ~ cauchy(0, 5);
  // log_normal solution
  to_vector(sigma) ~ lognormal(mean_trait, sigma_trait);
  mean_trait ~ student_t(1, 0, 10);
  sigma_trait ~ cauchy(0, 5);

  // priors on the dummy z's
  for(t in 1:(Ntc+Ntd)) {
    to_vector(z[t]) ~ normal(0, 1);
    grand_z[t, ] ~ normal(0, 1);
  }

  L_rho ~ lkj_corr_cholesky(2);
}


generated quantities {
  corr_matrix[Ntc+Ntd] rho;
  matrix[Ntc+Ntd, Ntc+Ntd] SD_cov;
  rho = multiply_lower_tri_self_transpose(L_rho);
  SD_cov = quad_form_diag(rho, grand_sigma_beta);
}
