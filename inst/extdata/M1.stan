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
  int Xk [Nk, Ns]; // design matrix (K x S)
}

parameters {
  vector [Ntq+Ntd] alpha;
  vector [Ntq+Ntd] grand_mu_beta;
  vector <lower = 0> [Ntq] sigma;
  vector <lower = 0> [Ntq+Ntd] sigma_beta;
  vector <lower = 0> [Ntq+Ntd] nu;
  vector <lower = 2> [Ntq+Ntd] nu_help;
  matrix [Ns, Nk] z [Ntq+Ntd];
  matrix [Ntq+Ntd, Ns] grand_z;

}

transformed parameters {
  matrix [Ns, Nk] beta [Ntq+Ntd];
  matrix [Ntq+Ntd, Ns] mu_beta;

  for(t in 1:(Ntq+Ntd)) {
    for(s in 1:Ns) {
       mu_beta[t][s] = grand_mu_beta[t] + sqrt(nu_help[t]/nu[t])*grand_z[t,s];
    }
  }


  for(t in 1:(Ntq+Ntd)) {
    for(s in 1:Ns) {
      for(k in 1:Nk) {
        beta[t][s, k] = mu_beta[t, s] + z[t][s, k]*sigma_beta[t];
      }
    }
  }
}

model {
  for(i in 1:N) {
    if(Ntq > 0) {
      for(t in 1:Ntq) {
        Yq[i,t] ~ normal(alpha[t] + X[i] .* beta[t][, K[i]], sigma[t]);
      }
    }
    if(Ntd > 0) {
      for(d in 1:Ntd) {
        Yd[i,d] ~ bernoulli_logit(alpha[d+Ntq] + X[i] .* beta[d+Ntq][, K[i]]);
      }
    }
  }

  alpha ~ student_t(1, 0, 100);
  to_vector(grand_mu_beta) ~ student_t(1, 0, 10);

  sigma ~ cauchy(0, 5);
  sigma_beta ~ cauchy(0, 5);
  (nu_help-2) ~ exponential(0.5);
  nu ~ chi_square(nu_help);

  // priors on the dummy z's
  for(t in 1:(Ntq+Ntd)) {
    to_vector(z[t]) ~ normal(0, 1);
    grand_z[t, ] ~ normal(0, 1);
  }
}

generated quantities {
  matrix [N, Ns] log_lik2 [Ntq+Ntd];
  matrix [N, Ntq+Ntd] log_lik;
  matrix [N, Ns] Yhat [Ntq+Ntd];
  matrix [Ns, Nk] Yhat_strain [Ntq+Ntd];
  matrix [Ns, 2] Yhat_snp [Ntq+Ntd];

  for(i in 1:N) {
    for(s in 1:Ns) {
      if(Ntq > 0) {
        for(t in 1:Ntq) {
          log_lik2[t][i,s] = normal_lpdf(Yq[i,t] | alpha[t]+X[i][s]*beta[t][s,K[i]], sigma[t]);
          Yhat[t][i,s] = normal_rng(alpha[t]+X[i][s]*beta[t][s,K[i]], sigma[t]);
        }
      }
      if(Ntd > 0) {
        for(d in 1:Ntd) {
          log_lik2[d+Ntq][i,s] = bernoulli_logit_lpmf(Yd[i,d] | alpha[d+Ntq]+X[i][s]*beta[d+Ntq][s,K[i]]);
          Yhat[Ntq+d][i,s] = bernoulli_rng(inv_logit(alpha[d+Ntq]+X[i][s]*beta[d+Ntq][s,K[i]]));
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
        Yhat_snp[t][s, 1] = alpha[t]+(+1)*mu_beta[t, s];
        Yhat_snp[t][s, 2] = alpha[t]+(-1)*mu_beta[t, s];

        // make strain-level predictions
        for(k in 1:Nk) {
           Yhat_strain[t][s, k] = alpha[t]+Xk[s, k]*beta[t][s, k];
        }
      }
    }
    if(Ntd > 0) {
      for(d in 1:Ntd) {
        Yhat_snp[Ntq+d][s, 1] = inv_logit(alpha[Ntq+d]+(+1)*mu_beta[Ntq+d, s]);
        Yhat_snp[Ntq+d][s, 2] = inv_logit(alpha[Ntq+d]+(-1)*mu_beta[Ntq+d, s]);

        // make strain-level predictions
        for(k in 1:Nk) {
           Yhat_strain[Ntq+d][s, k] = inv_logit(alpha[Ntq+d]+Xk[s, k]*beta[Ntq+d][s, k]);
        }
      }
    }
  }
}
