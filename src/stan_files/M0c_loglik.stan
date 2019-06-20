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
  vector [Ntq+Ntd] mu_beta;
  vector <lower = 0> [Ntq+Ntd] nu;
  vector <lower = 2> [Ntq+Ntd] nu_help;
  vector [Ntq+Ntd] z [Ns];
  cholesky_factor_corr [Ntq+Ntd] L_rho;
}

transformed parameters {
  matrix [Ntq+Ntd, Ns] beta;
  for(s in 1:Ns) {
    beta[, s] = mu_beta + sqrt(nu_help ./ nu) .* (L_rho * z[s]);
  }
}

model {
  for(i in 1:N) {
    if(Ntq > 0) {
      for(t in 1:Ntq) {
        Yq[i,t] ~ normal(alpha[t] + X[i] .* to_vector(beta[t, ]), sigma[t]);
      }
    }
    if(Ntd > 0) {
      for(d in 1:Ntd) {
        Yd[i,d] ~ bernoulli_logit(alpha[d+Ntq] + X[i] .* to_vector(beta[d+Ntq, ]));
      }
    }
  }

  alpha ~ student_t(1, 0, 100);
  mu_beta ~ student_t(1, 0, 10);
  sigma ~ cauchy(0, 5);
  (nu_help-2) ~ exponential(0.5);
  nu ~ chi_square(nu_help);

  for(s in 1:Ns) {
    z[s] ~ normal(0, 1);
  }

  L_rho ~ lkj_corr_cholesky(2);
}


generated quantities {
  matrix [N, Ns] log_lik2 [Ntq+Ntd];
  matrix [N, Ntq+Ntd] log_lik;
  corr_matrix[Ntq+Ntd] rho;
  matrix [N, Ns] Yhat [Ntq+Ntd];
  rho = multiply_lower_tri_self_transpose(L_rho);

  for(i in 1:N) {
    for(s in 1:Ns) {
      if(Ntq > 0) {
        for(t in 1:Ntq) {
          log_lik2[t][i,s] = normal_lpdf(Yq[i,t] | alpha[t]+X[i][s]*beta[t, s], sigma[t]);
          Yhat[t][i,s] = normal_rng(alpha[t]+X[i][s]*beta[t][s], sigma[t]);
        }
      }
      if(Ntd > 0) {
        for(d in 1:Ntd) {
          log_lik2[Ntq+d][i,s] = bernoulli_logit_lpmf(Yd[i, d] | alpha[Ntq+d]+X[i][s]*beta[Ntq+d, s]);
          Yhat[Ntq+d][i,s] = bernoulli_rng(inv_logit(alpha[Ntq+d]+X[i][s]*beta[Ntq+d][s]));
        }
      }
    }

    for(t in 1:(Ntq+Ntd)) {
      log_lik[i, t] = mean(log_lik2[t][i, ]);
    }
  }
}
