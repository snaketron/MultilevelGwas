data {
  int N;                // number of all entries
  int Ns;               // number of all SNPs
  int Nk;               // number of all strains
  real Y [N];           // number of hits response
  vector [Ns] X [N];    // index of all individuals
  int K[N];             //index to all strains
}


parameters {
  real alpha_trait;
  real beta_trait;
  real <lower = 0> sigma;
  real <lower = 0> sigma_snp;
  real <lower = 0> sigma_trait;
  real <lower = 0> nu_trait;
  vector <offset=beta_trait, multiplier=sigma_trait> [Ns] beta_snp;
  vector [Ns] z [Nk];
}


transformed parameters {
  vector [Ns] beta_strain [Nk];

  for(k in 1:Nk) {
    beta_strain[k] = beta_snp + z[k]*sigma_snp;
  }
}

model {
  for(i in 1:N) {
    Y[i] ~ normal(alpha_trait + X[i] .* beta_strain[K[i]], sigma);
  }

  alpha_trait ~ normal(0, 20);
  beta_trait ~ normal(0, 5);

  beta_snp ~ student_t(nu_trait, beta_trait, sigma_trait);

  sigma ~ cauchy(0, 5);
  sigma_snp ~ cauchy(0, 5);
  sigma_trait ~ cauchy(0, 5);
  nu_trait ~ gamma(2, 0.1);

  for(k in 1:Nk) {
    z[k] ~ std_normal();
  }
}
