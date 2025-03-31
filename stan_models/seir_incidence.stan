functions {
 real[] seir(real t, real[] y, real[] theta,
            real[] x_r,  int[] x_i) {
    real S = y[1];
    real E = y[2];
    real I = y[3];
    real R = y[4];
    real N = x_i[1];

    real beta = theta[1];
    real gamma = theta[2];
    real sigma = theta[3];

    real dS_dt = -beta * I * S / N;
    real dE_dt = beta * I * S / N - gamma * E;
    real dI_dt =  gamma * E - sigma * I;
    real dR_dt =  sigma * I;
    real dC_dt = gamma * E;

    return {dS_dt, dE_dt, dI_dt, dR_dt, dC_dt};
  }
}

data {
  int<lower=1> n_days;
  real y0[5];
  real t0;
  real ts[n_days];
  int N;
  int cases[n_days];
}

transformed data {
  real x_r[0];
  int x_i[1] = {N};
}

parameters {
  real<lower=0> beta;
  real<lower=0> gamma;
  real<lower=0> sigma;
  real<lower=0> phi_inv;
}

transformed parameters{
  real y[n_days, 5];
  real phi = 1. / phi_inv;
  real incidence[n_days];
  real theta[3]; // model parameters

  theta[1] = beta;
  theta[2] = gamma;
  theta[3] = sigma;

  y = integrate_ode_rk45(seir, y0, t0, ts, theta, x_r, x_i);

  incidence[1] = y[1, 5];
  for (i in 2:n_days)
    incidence[i] = y[i, 5] - y[i-1, 5];

}

model {
  //priors
  beta ~ normal(4, 2); //truncated at 0 2,1
  gamma ~ normal(0.9, 0.3); //truncated at 0
  sigma ~ normal(0.5, 0.25); //truncated at 0
  phi_inv ~ exponential(5);

  //likelihood
  cases ~ neg_binomial_2(incidence, phi);
}

generated quantities {
  real recovery_time = 1 / sigma;
  real incubation_period = 1 / gamma;
}
