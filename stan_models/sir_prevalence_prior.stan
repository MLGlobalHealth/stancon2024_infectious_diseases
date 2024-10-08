functions {
 real[] sir(real t, real[] y, real[] theta,
            real[] x_r,  int[] x_i) {
    real S = y[1];
    real I = y[2];
    real R = y[3];
    real N = x_i[1];

    real beta = theta[1];
    real sigma = theta[2];

    real dS_dt = -beta * I * S / N;
    real dR_dt =  sigma * I;
    real dI_dt =  -(dS_dt + dR_dt);

    return {dS_dt, dI_dt, dR_dt};
  }
}

data {
  int<lower=1> n_days;
  real y0[3];
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
  real<lower=0> sigma;
  real<lower=0> beta;
  real<lower=0> phi_inv;
}

transformed parameters{
  real y[n_days, 3];
  real phi = 1. / phi_inv;
  {
    real theta[2]; // model parameters
    theta[1] = beta;
    theta[2] = sigma;

    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
  }
}

model {
  //priors
  beta ~ normal(2, 1); //truncated at 0
  sigma ~ normal(0.4, 0.25); //truncated at 0
  phi_inv ~ exponential(5);

}

generated quantities {
  real recovery_time = 1 / sigma;
  real pred_cases[n_days];
  pred_cases = neg_binomial_2_rng(col(to_matrix(y), 2) + 1e-5, phi);
}
