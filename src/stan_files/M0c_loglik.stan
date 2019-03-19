data {
  int N; // number of all entries
  int Ntc; // number of continuous traits
  int Ntd; // number of dichotomous traits
  int Ns; // number of all SNPs
  real Yc[N, Ntc] ; // number of hits response
  int Yd[N, Ntd] ; // number of hits response
  vector [Ns] X [N]; // index of all individuals
}

parameters {
  vector [Ntc+Ntd] alpha;
  vector <lower = 0> [Ntc] sigma;
  matrix [Ntc+Ntd, 1] mu_beta;
  vector <lower = 0> [Ntc+Ntd] sd_beta;
  cholesky_factor_corr [Ntc+Ntd] L_rho;
  matrix [Ntc+Ntd, Ns] z;
}

transformed parameters {
  matrix [Ntc+Ntd, Ns] beta;
  for(s in 1:Ns) {
     beta[, s] = mu_beta[, 1] + diag_pre_multiply(sd_beta, L_rho)*z[, s]; // multi_normal
  }
}

model {
  for(i in 1:N) {
    if(Ntc > 0) {
      for(t in 1:Ntc) {
        Yc[i,t] ~ normal(alpha[t] + rows_dot_product(X[i], to_vector(beta[t, ])), sigma[t]);
      }
    }
    if(Ntd > 0) {
      for(d in 1:Ntd) {
        Yd[i,d] ~ bernoulli_logit(alpha[d+Ntc] + rows_dot_product(X[i], to_vector(beta[d+Ntc, ])));
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
  matrix [N, Ns] log_lik [Ntc+Ntd];
  corr_matrix[Ntc+Ntd] rho;
  matrix[Ntc+Ntd, Ntc+Ntd] SD_cov;

  rho = multiply_lower_tri_self_transpose(L_rho);
  SD_cov = quad_form_diag(rho, sd_beta);

  for(s in 1:Ns) {
    for(i in 1:N) {
      for(t in 1:Ntc) {
        log_lik[t][i,s] = normal_lpdf(Yc[i,t] | alpha[t]+X[i][s]*beta[t, s], sigma[t]);
      }

      for(d in 1:Ntd) {
        log_lik[d+Ntc][i,s] = bernoulli_logit_lpmf(Yd[i, d] | alpha[Ntc+d]+X[i][s]*beta[Ntc+d, s]);
      }
    }
  }
}
