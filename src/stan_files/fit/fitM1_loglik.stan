data {
  int N; // number of all entries
  int Ntq; // number of traits
  int Ntd; // number of dichotomous traits
  int Nk; // number of all strains
  real Yq[N, Ntq] ; // number of hits response
  int Yd[N, Ntd] ; // number of hits response
  int K[N]; //index to all strains
}

parameters {
  vector <lower = 0> [Nk] sigma;
  real <lower = 0> sigma_mu;
  real <lower = 0> sigma_sigma;
  vector <lower = 0> [Ntq+Ntd] grand_sigma;
  vector [Ntq+Ntd] grand_mu;
  matrix [Ntq+Ntd, Nk] z;
}

transformed parameters {
  matrix [Ntq+Ntd, Nk] mu;

  for(t in 1:(Ntq+Ntd)) {
    for(k in 1:Nk) {
      mu[t, k] = grand_mu[t] + z[t, k]*grand_sigma[t];
    }
  }
}

model {
  for(i in 1:N) {
    if(Ntq > 0) {
      for(t in 1:Ntq) {
        Yq[i,t] ~ normal(mu[t, K[i]], sigma[K[i]]);
      }
    }

    if(Ntd > 0) {
      for(d in 1:Ntd) {
        Yd[i,d] ~ bernoulli_logit(mu[d+Ntq, K[i]]);
      }
    }
  }

  grand_mu ~ student_t(1, 0, 10);

  to_vector(sigma) ~ lognormal(sigma_mu, sigma_sigma);
  sigma_mu ~ student_t(1, 0, 10);
  sigma_sigma ~ cauchy(0, 5);

  grand_sigma ~ cauchy(0, 5);
  to_vector(z) ~ normal(0, 1);
}

generated quantities {
  real log_lik [N, Ntq+Ntd];
  for(i in 1:N) {

    if(Ntq > 0) {
      for(t in 1:Ntq) {
        log_lik[i, t] = normal_lpdf(Yq[i,t] | mu[t, K[i]], sigma[K[i]]);
      }
    }

    if(Ntd > 0) {
      for(d in 1:Ntd) {
        log_lik[i, (Ntq+d)] = bernoulli_logit_lpmf(Yd[i,d] | mu[(Ntq+d), K[i]]);
      }
    }
  }
}
