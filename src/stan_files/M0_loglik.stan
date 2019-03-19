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
  matrix [Ntq+Ntd, Ns] z;
}

transformed parameters {
  matrix [Ntq+Ntd, Ns] beta;

  for(t in 1:(Ntq+Ntd)) {
    for(s in 1:Ns) {
      beta[t, s] = mu_beta[t, 1] + z[t, s]*sd_beta[t];
    }
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

  for(s in 1:Ns) {
    z[, s] ~ normal(0, 1);
  }
}

generated quantities {
  matrix [N, Ns] log_lik [Ntq+Ntd];

  for(i in 1:N) {
    for(s in 1:Ns) {
      if(Ntq > 0) {
        for(t in 1:Ntq) {
          log_lik[t][i,s] = normal_lpdf(Yq[i,t] | alpha[t]+X[i][s]*beta[t, s], sigma[t]);
        }
      }
      if(Ntd > 0) {
        for(d in 1:Ntd) {
          log_lik[Ntq+d][i,s] = bernoulli_logit_lpmf(Yd[i, d] | alpha[Ntq+d]+X[i][s]*beta[Ntq+d, s]);
        }
      }
    }
  }
}
