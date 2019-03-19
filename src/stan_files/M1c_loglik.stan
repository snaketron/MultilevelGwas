data {
  int N; // number of all entries
  int Ntq; // number of traits
  int Ntd; // number of traits
  int Ns; // number of all SNPs
  int Nk; // number of all strains
  real Yq[N, Ntq] ; // number of hits response
  int Yd[N, Ntd] ; // number of hits response
  vector [Ns] X [N]; // index of all individuals
  int K[N]; //index to all strains
}

parameters {
  vector [Ntq+Ntd] alpha;
  matrix [Ntq+Ntd, 1] grand_mu_beta;
  vector <lower = 0> [Ntq] sigma;
  vector <lower = 0> [Ntq+Ntd] sigma_beta;
  vector <lower = 0> [Ntq+Ntd] grand_sigma_beta;
  matrix [Ns, Nk] z [Ntq+Ntd];
  matrix [Ntq+Ntd, Ns] grand_z;
  cholesky_factor_corr [Ntq+Ntd] L_rho;
}

transformed parameters {
  matrix [Ns, Nk] beta [Ntq+Ntd];
  matrix [Ntq+Ntd, Ns] mu_beta;


  for(s in 1:Ns) {
    // multi_normal
      mu_beta[, s] = grand_mu_beta[, 1] + diag_pre_multiply(grand_sigma_beta, L_rho)*grand_z[, s];
    }


  for(t in 1:(Ntq+Ntd)) {
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
    if(Ntq > 0) {
      for(t in 1:Ntq) {
        Yq[i,t] ~ normal(alpha[t] + rows_dot_product(X[i], to_vector(beta[t][, K[i]])), sigma[t]);
      }
    }
    if(Ntd > 0) {
      for(d in 1:Ntd) {
        Yd[i,d] ~ bernoulli_logit(alpha[d+Ntq] + rows_dot_product(X[i], to_vector(beta[d+Ntq][, K[i]])));
      }
    }
  }

  alpha ~ student_t(1, 0, 100);
  to_vector(grand_mu_beta) ~ student_t(1, 0, 10);

  sigma ~ cauchy(0, 5);
  sigma_beta ~ cauchy(0, 5);
  grand_sigma_beta ~ cauchy(0, 5);

  for(t in 1:(Ntq+Ntd)) {
    to_vector(z[t]) ~ normal(0, 1);
    grand_z[t, ] ~ normal(0, 1);
  }

  L_rho ~ lkj_corr_cholesky(2);
}


generated quantities {
  matrix [N, Ns] log_lik [Ntq+Ntd];
  corr_matrix[(Ntq+Ntd)] rho;
  matrix[(Ntq+Ntd), (Ntq+Ntd)] SD_cov;

  rho = multiply_lower_tri_self_transpose(L_rho);
  SD_cov = quad_form_diag(rho, grand_sigma_beta);

  for(s in 1:Ns) {
    for(i in 1:N) {
      if(Ntq > 0) {
        for(t in 1:Ntq) {
          log_lik[t][i,s] = normal_lpdf(Yq[i,t] | alpha[t]+X[i][s]*beta[t][s,K[i]], sigma[t]);
        }
      }
      if(Ntd > 0) {
        for(d in 1:Ntd) {
          log_lik[d+Ntq][i,s] = bernoulli_logit_lpmf(Yd[i,d] | alpha[d+Ntq]+X[i][s]*beta[d+Ntq][s,K[i]]);
        }
      }
    }
  }
}
