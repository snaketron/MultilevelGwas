data {
  int N;                    // number of all entries
  int Ns;                   // number of all SNPs
  real Y [N];               // number of hits response
  vector [Ns] X [N];        // index of all individuals
  int <lower=0, upper=1> P; // trait type: 1=cont, 0=dich
}


parameters {
  real alpha_trait;
  real beta_trait;
  real <lower = 0> sigma [P];
  real <lower = 0> sigma_trait;
  real <lower = 0> nu_trait;
  vector <offset = beta_trait, multiplier = sigma_trait> [Ns] beta_snp;
}


model {
  int Yint;
  if(P == 1) {
    for(i in 1:N) {
      Y[i] ~ normal(alpha_trait + X[i] .* beta_snp, sigma);
    }
    sigma ~ cauchy(0, 5);
  }
  if(P == 0) {
    for(i in 1:N) {
      Y[i]>0?1:0 ~ bernoulli_logit(alpha_trait + X[i] .* beta_snp);
    }
  }

  alpha_trait ~ normal(0, 20);
  beta_trait ~ normal(0, 5);

  beta_snp ~ student_t(nu_trait, beta_trait, sigma_trait);

  sigma_trait ~ cauchy(0, 5);
  nu_trait ~ gamma(2, 0.1);
}

