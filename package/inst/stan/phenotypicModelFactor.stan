// Multivariate linear model for the phenotypic covariance matrix P with a
// sparse factor structure on the correlation scale.
// P = diag(y_sd) (Lambda Lambda' + Psi) diag(y_sd)
data {
  int<lower=1> K; // number of traits
  int<lower=1> J; // number of fixed effects
  int<lower=2> N; // number of individuals
  int<lower=1, upper=K - 1> M; // number of factors
  matrix[N, J] X; // fixed effects design matrix
  matrix[N, K] Y; // response variables
  vector[J] x_mean_weight; // X * x_mean_weight is a column of ones
  real<lower=0> hs_global_scale; // half-Cauchy scale for the global shrinkage
  real<lower=0> hs_slab_scale; // slab scale of the regularized horseshoe
  real<lower=0, upper=1> column_decay; // sd multiplier for successive factors
}
transformed data {
  vector[K] y_sd;
  vector[K] y_mean;
  vector[J] x_center;
  vector[J] x_scale;
  vector[J] x_scale_inv;
  vector[J] beta_prior_sd;
  array[J] int x_is_constant;
  vector[J] x_constant_value;
  matrix[N, K] Ys;
  matrix[N, J] Xs;
  vector[M] column_scale;

  for (m in 1:M)
    column_scale[m] = pow(column_decay, m - 1);

  for (k in 1:K) {
    y_sd[k] = sd(Y[, k]);
    y_mean[k] = mean(Y[, k]);
    Ys[, k] = (Y[, k] - y_mean[k]) / y_sd[k];
  }
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
    Xs[, j] = (X[, j] - x_center[j]) * x_scale_inv[j];
  }
}
parameters {
  matrix[K, J] beta_scaled; // fixed effects on scaled predictors
  matrix[K, M] loadings_raw; // non-centred loadings
  matrix<lower=0>[K, M] local_scale; // local horseshoe scales
  real<lower=0> global_scale; // global horseshoe scale
}
transformed parameters {
  matrix[K, J] beta; // fixed effects on original predictor and response scales
  matrix[K, M] loadings_tp;
  vector<lower=0>[K] uniqueness_tp;
  vector<lower=0>[K] communality_tp;
  for (k in 1:K) {
    for (j in 1:J) {
      beta[k, j] = y_sd[k] * beta_scaled[k, j] * x_scale_inv[j];
    }
  }
  for (k in 1:K) {
    real mean_shift = y_mean[k] - dot_product(x_center, beta[k]');
    for (j in 1:J)
      beta[k, j] += mean_shift * x_mean_weight[j];
  }
  // Regularized horseshoe with deterministic column decay, followed by the
  // row normalization that makes the factor part a correlation matrix.
  for (k in 1:K) {
    real raw_sq = 0;
    for (m in 1:M) {
      real local_variance = fmin(
        square(global_scale) * square(local_scale[k, m]) *
        square(column_scale[m]), 1e150);
      real slab_variance = square(hs_slab_scale);
      real shrunk_variance = slab_variance * local_variance /
        (slab_variance + local_variance);
      raw_sq += shrunk_variance * square(loadings_raw[k, m]);
    }
    real row_scale = inv_sqrt(1 + raw_sq);
    uniqueness_tp[k] = square(row_scale);
    communality_tp[k] = 1 - uniqueness_tp[k];
    for (m in 1:M) {
      real local_variance = fmin(
        square(global_scale) * square(local_scale[k, m]) *
        square(column_scale[m]), 1e150);
      real slab_variance = square(hs_slab_scale);
      real shrunk_variance = slab_variance * local_variance /
        (slab_variance + local_variance);
      loadings_tp[k, m] = row_scale * loadings_raw[k, m] * sqrt(shrunk_variance);
    }
  }
}
model {
  matrix[K, K] P_std = loadings_tp * loadings_tp' + diag_matrix(uniqueness_tp);
  matrix[K, K] L_P = cholesky_decompose(P_std);
  vector[K] log_diag_P = log(diagonal(L_P));
  for (k in 1:K)
    L_P[k, k] = fmax(L_P[k, k], 1e-12);

  matrix[N, K] resid = Ys - Xs * beta_scaled';
  matrix[K, N] W = mdivide_left_tri_low(L_P, resid');
  target += -0.5 * (N * K * log(2 * pi())
                    + 2 * N * sum(log_diag_P)
                    + sum(square(W)));

  for (k in 1:K)
    for (j in 1:J)
      beta_scaled[k, j] ~ normal(0, beta_prior_sd[j]);
  to_vector(loadings_raw) ~ normal(0, 1);
  to_vector(local_scale) ~ cauchy(0, 1);
  global_scale ~ cauchy(0, hs_global_scale);
}
generated quantities {
  matrix[K, M] loadings = loadings_tp;
  vector[K] uniqueness = uniqueness_tp;
  vector[K] communality = communality_tp;
  vector[K] sigma_P = y_sd;
  matrix[K, K] P_std = loadings_tp * loadings_tp' + diag_matrix(uniqueness_tp);
  matrix[K, K] P = quad_form_diag(P_std, y_sd);
  matrix[K, K] corrP = P_std;
  vector[M] factor_magnitude;
  for (m in 1:M)
    factor_magnitude[m] = dot_product(loadings_tp[, m], loadings_tp[, m]);
}
