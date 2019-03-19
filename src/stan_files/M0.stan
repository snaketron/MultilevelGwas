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
  matrix [Ntc+Ntd, Ns] z;
}

transformed parameters {
  matrix [Ntc+Ntd, Ns] beta;

  for(t in 1:(Ntc+Ntd)) {
    for(s in 1:Ns) {
      beta[t, s] = mu_beta[t, 1] + z[t, s]*sd_beta[t];
    }
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

  for(s in 1:Ns) {
    z[, s] ~ normal(0, 1);
  }
}
