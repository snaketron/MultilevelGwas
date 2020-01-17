data {
  int N;                    // number of all entries
  int Ns;                   // number of all SNPs
  int Nk;                   // number of all strains
  real Y [N];               // number of hits response
  vector [Ns] X [N];        // index of all individuals
  int K[N];                 //index to all strains
  int <lower=0, upper=1> P; // trait type: 1=cont, 0=dich
}


parameters {
  real alpha_trait;
  real <lower = 0> sigma [P];
  real <lower = 0> sigma_snp;
  real <lower = 0> sigma_trait;
  real <lower = 1> nu_trait;
  vector <offset=0, multiplier=sigma_trait> [Ns] beta_snp;
  vector [Ns] z [Nk];
}


transformed parameters {
  vector [Ns] beta_strain [Nk];

  for(k in 1:Nk) {
    beta_strain[k] = beta_snp + z[k]*sigma_snp;
  }
}

model {
  int Yint;
  if(P == 1) {
    for(i in 1:N) {
      Y[i] ~ normal(alpha_trait + X[i] .* beta_strain[K[i]], sigma[1]);
    }
    sigma[1] ~ cauchy(0, 5);
  }
  if(P == 0) {
    for(i in 1:N) {
      Y[i]>0?1:0 ~ bernoulli_logit(alpha_trait + X[i] .* beta_strain[K[i]]);
    }
  }

  alpha_trait ~ normal(0, 20);
  beta_snp ~ student_t(nu_trait, 0, sigma_trait);
  sigma_snp ~ cauchy(0, 5);
  sigma_trait ~ cauchy(0, 5);
  nu_trait ~ gamma(2, 0.1);

  for(k in 1:Nk) {
    z[k] ~ std_normal();
  }
}
