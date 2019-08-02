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
  vector [Ns] z [Ntq+Ntd];
}

transformed parameters {
  vector [Ns] beta [Ntq+Ntd];

  for(t in 1:(Ntq+Ntd)) {
    beta[t] = mu_beta[t] + sqrt(nu_help[t]/nu[t]) * z[t];
  }
}

model {
  for(i in 1:N) {
    if(Ntq > 0) {
      for(t in 1:Ntq) {
        Yq[i,t] ~ normal(alpha[t] + X[i] .* beta[t], sigma[t]);
      }
    }
    if(Ntd > 0) {
      for(d in 1:Ntd) {
        Yd[i,d] ~ bernoulli_logit(alpha[d+Ntq] + X[i] .* beta[d+Ntq]);
      }
    }
  }

  alpha ~ student_t(1, 0, 100);
  mu_beta ~ student_t(1, 0, 10);
  sigma ~ cauchy(0, 5);
  (nu_help-2) ~ exponential(0.5);
  nu ~ chi_square(nu_help);

  for(t in 1:(Ntq+Ntd)) {
    z[t] ~ normal(0, 1);
  }
}

generated quantities {
  matrix [N, Ns] log_lik2 [Ntq+Ntd];
  matrix [N, Ntq+Ntd] log_lik;
  matrix [N, Ns] Yhat [Ntq+Ntd];
  matrix [2, Ns] Yhat_snp [Ntq+Ntd];

  for(i in 1:N) {
    for(s in 1:Ns) {
      if(Ntq > 0) {
        for(t in 1:Ntq) {
          log_lik2[t][i,s] = normal_lpdf(Yq[i,t] | alpha[t]+X[i][s]*beta[t][s], sigma[t]);
          Yhat[t][i,s] = normal_rng(alpha[t]+X[i][s]*beta[t][s], sigma[t]);
        }
      }
      if(Ntd > 0) {
        for(d in 1:Ntd) {
          log_lik2[Ntq+d][i,s] = bernoulli_logit_lpmf(Yd[i, d] | alpha[Ntq+d]+X[i][s]*beta[Ntq+d][s]);
          Yhat[Ntq+d][i,s] = bernoulli_rng(inv_logit(alpha[Ntq+d]+X[i][s]*beta[Ntq+d][s]));
        }
      }
    }

    for(t in 1:(Ntq+Ntd)) {
      log_lik[i, t] = mean(log_lik2[t][i, ]);
    }
  }

  // make SNP-level predictions
  for(s in 1:Ns) {
    if(Ntq > 0) {
      for(t in 1:Ntq) {
        Yhat_snp[t][1, s] = alpha[t]+(+1)*beta[t][s];
        Yhat_snp[t][2, s] = alpha[t]+(-1)*beta[t][s];
      }
    }
    if(Ntd > 0) {
      for(d in 1:Ntd) {
        Yhat_snp[Ntq+d][1, s] = inv_logit(alpha[Ntq+d]+(+1)*beta[Ntq+d][s]);
        Yhat_snp[Ntq+d][2, s] = inv_logit(alpha[Ntq+d]+(-1)*beta[Ntq+d][s]);
      }
    }
  }
}
