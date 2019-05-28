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
    mu_beta[t, ] = grand_mu_beta[t] + sqrt(nu_help[t]/nu[t])*grand_z[t, ];
  }

<<<<<<< HEAD
  for(t in 1:(Ntq+Ntd)) {
    for(k in 1:Nk) {
      beta[t][,k] = to_vector(mu_beta[t, ]) + z[t][,k]*sigma_beta[t];
=======

  for(t in 1:(Ntq+Ntd)) {
    for(s in 1:Ns) {
      beta[t][s,] = mu_beta[t, s] + z[t][s,]*sigma_beta[t];
>>>>>>> a58ac4fc2c528fee7484196f0d57850ec6e2f01a
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
  }
  to_vector(grand_z) ~ normal(0, 1);
}
