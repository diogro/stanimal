// Univariate animal model parametrized by the genetic and residual variances.
// Breeding values are marginalized: with A = Q diag(lambda) Q', the rotated
// records are independent, y*_i ~ N(X*_i beta, lambda_i sigma_G^2 + sigma_E^2).
data {
  int<lower=1> J; // number of fixed effects
  int<lower=2> N; // number of individuals
  matrix[N, J] X; // fixed effects design matrix
  vector[N] Y; // response variable, with missing entries set to zero
  array[N] int<lower=0, upper=1> miss; // 1 where the record is missing
  int<lower=0> n_missing; // number of gaps; each is one parameter
  array[n_missing] int<lower=1, upper=N> mis_row; // row of each gap
  matrix[N, N] A; // known relationship matrix (symmetric, PSD)
  vector[J] x_mean_weight; // X * x_mean_weight is a column of ones
}
transformed data {
  vector[N] lambda;
  matrix[N, N] Q = eigenvectors_sym(A);
  real y_sd;
  real y_mean;
  vector[J] x_center;
  vector[J] x_scale;
  vector[J] x_scale_inv;
  vector[J] beta_prior_sd;
  array[J] int x_is_constant;
  vector[J] x_constant_value;
  vector[N] Ys_observed; // rotated standardized response, gaps at zero
  matrix[N, J] Xs;

  {
    vector[N] ev = eigenvalues_sym(A);
    for (n in 1:N)
      lambda[n] = fmax(ev[n], 0);
  }

  {
    matrix[N, J] X_std;
    for (j in 1:J) {
      int is_indicator = 1;
      int is_constant = 1;
      for (n in 1:N) {
        if (X[n, j] != X[1, j])
          is_constant = 0;
        if (X[n, j] != 0 && X[n, j] != 1)
          is_indicator = 0;
      }
      x_is_constant[j] = is_constant;
      x_constant_value[j] = X[1, j];
      x_center[j] = 0;
      x_scale[j] = 1;
      if (is_indicator == 0 && is_constant == 0) {
        x_center[j] = mean(X[, j]);
        x_scale[j] = sd(X[, j]);
      }
      x_scale_inv[j] = 1 / x_scale[j];
      beta_prior_sd[j] = 1;
      if (is_constant == 1 && x_constant_value[j] != 0)
        beta_prior_sd[j] = 2.5 / abs(x_constant_value[j]);
      X_std[, j] = (X[, j] - x_center[j]) * x_scale_inv[j];
    }
    // Centre and scale on the observed records only, then rotate. A gap is a
    // parameter: the rotation is linear, so its contribution is added back per
    // draw in the model block.
    {
      int present = 0;
      real total = 0;
      for (n in 1:N)
        if (miss[n] == 0) {
          present += 1;
          total += Y[n];
        }
      y_mean = total / present;
      real sum_sq = 0;
      for (n in 1:N)
        if (miss[n] == 0)
          sum_sq += square(Y[n] - y_mean);
      y_sd = sqrt(sum_sq / (present - 1));
      vector[N] Y_std;
      for (n in 1:N)
        Y_std[n] = miss[n] == 1 ? 0 : (Y[n] - y_mean) / y_sd;
      Ys_observed = Q' * Y_std;
    }
    Xs = Q' * X_std;
  }
}
parameters {
  vector[n_missing] Y_mis; // the gaps, on the standardized response scale
  vector[J] beta_scaled; // fixed effects on scaled predictors
  real<lower=0> s_G; // genetic sd, standardized response scale
  real<lower=0> s_E; // residual sd, standardized response scale
}
transformed parameters {
  vector[J] beta; // fixed effects on original predictor and response scales
  for (j in 1:J) {
    beta[j] = y_sd * beta_scaled[j] * x_scale_inv[j];
  }
  // The response is centred and scaled, so the model predicts deviations from
  // its mean. Put the mean back where the design can hold it: X x_mean_weight is
  // a column of ones, so the shift lands on every prediction, whether the design
  // has an intercept column or the indicators of a factor.
  {
    real mean_shift = y_mean - dot_product(x_center, beta);
    for (j in 1:J)
      beta[j] += mean_shift * x_mean_weight[j];
  }
}
model {
  // A gap in individual n shifts every rotated record by Q[n, i] * value, so
  // the rotated response is assembled here rather than in transformed data.
  vector[N] Ys = Ys_observed;
  for (j in 1:n_missing)
    Ys += Q[mis_row[j]]' * Y_mis[j];

  Ys ~ normal(Xs * beta_scaled, sqrt(lambda * square(s_G) + square(s_E)));

  for (j in 1:J)
    beta_scaled[j] ~ normal(0, beta_prior_sd[j]);
  s_G ~ normal(0, 1);
  s_E ~ normal(0, 1);
}
generated quantities {
  real sigma_G = y_sd * s_G;
  real sigma_E = y_sd * s_E;
  real h2 = square(s_G) / (square(s_G) + square(s_E));
  // Pointwise log likelihood: the rotated records are independent, so the terms
  // sum to the joint log density. Term n belongs to rotated record n, which is
  // individual n only when A has a trivial eigenbasis.
  vector[N] log_lik;
  {
    vector[N] Ys = Ys_observed;
    for (j in 1:n_missing)
      Ys += Q[mis_row[j]]' * Y_mis[j];
    vector[N] mu = Xs * beta_scaled;
    vector[N] sd_rotated = sqrt(lambda * square(s_G) + square(s_E));
    for (n in 1:N)
      log_lik[n] = normal_lpdf(Ys[n] | mu[n], sd_rotated[n]);
  }
}
