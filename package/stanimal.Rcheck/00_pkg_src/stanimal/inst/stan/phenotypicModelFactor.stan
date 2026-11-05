// Multivariate linear model for the phenotypic covariance matrix P with a
// sparse factor structure on the correlation scale.
// P = diag(sigma_P) (Lambda Lambda' + Psi) diag(sigma_P), sigma_P = y_sd .* s_P
// The scale s_P is a parameter: P describes the residuals, and removing the
// fixed effects lowers their variance below the marginal variance of Y. Fixing
// the diagonal at the marginal variance would leave the correlations as the
// only free part and inflate them to compensate. (ask me how I know...)
//
// Incomplete traits are handled by pattern, as in phenotypicModel.stan: each
// pattern evaluates the density on its observed subvector of P, so the missing
// values are integrated out and cost no parameters.
data {
  int<lower=1> K; // number of traits
  int<lower=1> J; // number of fixed effects
  int<lower=2> N; // number of individuals
  int<lower=1, upper=K - 1> M; // number of factors
  matrix[N, J] X; // fixed effects design matrix
  matrix[N, K] Y; // response variables, with missing entries set to zero
  array[N, K] int<lower=0, upper=1> miss; // 1 where that trait is missing
  int<lower=1, upper=N> n_pattern; // number of distinct missingness patterns
  array[n_pattern, K] int<lower=0, upper=1> pattern_obs; // observed traits per pattern
  array[N] int<lower=1, upper=N> row_order; // rows grouped by pattern
  array[n_pattern] int<lower=0> n_rows; // rows in each pattern
  array[n_pattern] int<lower=1> first_row; // first reordered row of each pattern
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
  matrix[N, K] Ys; // standardized responses, rows grouped by pattern
  matrix[N, J] Xs; // transformed design, rows grouped by pattern
  vector[M] column_scale;

  for (m in 1:M)
    column_scale[m] = pow(column_decay, m - 1);

  {
    matrix[N, K] Yp; // the response in pattern-contiguous order
    matrix[N, J] Xp;
    for (n in 1:N) {
      Yp[n] = Y[row_order[n]];
      Xp[n] = X[row_order[n]];
    }
    for (k in 1:K) {
      int present = 0;
      real total = 0;
      for (n in 1:N)
        if (miss[row_order[n], k] == 0) {
          present += 1;
          total += Yp[n, k];
        }
      y_mean[k] = total / present;
      real sum_sq = 0;
      for (n in 1:N)
        if (miss[row_order[n], k] == 0)
          sum_sq += square(Yp[n, k] - y_mean[k]);
      y_sd[k] = sqrt(sum_sq / (present - 1));
      for (n in 1:N)
        Ys[n, k] = miss[row_order[n], k] == 1 ? 0 : (Yp[n, k] - y_mean[k]) / y_sd[k];
    }
    for (j in 1:J) {
      int is_indicator = 1;
      int is_constant = 1;
      for (n in 1:N) {
        if (Xp[n, j] != Xp[1, j])
          is_constant = 0;
        if (Xp[n, j] != 0 && Xp[n, j] != 1)
          is_indicator = 0;
      }
      x_is_constant[j] = is_constant;
      x_constant_value[j] = Xp[1, j];
      x_center[j] = 0;
      x_scale[j] = 1;
      if (is_indicator == 0 && is_constant == 0) {
        x_center[j] = mean(Xp[, j]);
        x_scale[j] = sd(Xp[, j]);
      }
      x_scale_inv[j] = 1 / x_scale[j];
      beta_prior_sd[j] = 1;
      if (is_constant == 1 && x_constant_value[j] != 0)
        beta_prior_sd[j] = 2.5 / abs(x_constant_value[j]);
      Xs[, j] = (Xp[, j] - x_center[j]) * x_scale_inv[j];
    }
  }
}
parameters {
  matrix[K, J] beta_scaled; // fixed effects on scaled predictors
  vector<lower=0>[K] s_P; // phenotypic sds, standardized response scale
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

  for (p in 1:n_pattern) {
    // Traits this pattern observes.
    int m = 0;
    array[K] int obs;
    for (k in 1:K)
      if (pattern_obs[p, k] == 1) {
        m += 1;
        obs[m] = k;
      }
    // Principal submatrix of P on that pattern's observed traits, on the
    // standardized response scale. P_std is a correlation matrix, so the scale
    // s_P enters as the two-sided multiplication below. The diagonal is floored
    // so that an extreme proposal gives a finite, hugely negative log density
    // rather than an exception out of the Cholesky.
    matrix[m, m] P_sub = P_std[obs[1:m], obs[1:m]];
    for (i in 1:m)
      for (j in 1:m)
        P_sub[i, j] *= s_P[obs[i]] * s_P[obs[j]];
    for (i in 1:m)
      P_sub[i, i] = fmax(P_sub[i, i], 1e-12);
    matrix[m, m] L_sub = cholesky_decompose(P_sub);

    int first = first_row[p];
    int last = first + n_rows[p] - 1;
    matrix[n_rows[p], m] resid = Ys[first:last, obs[1:m]]
                                 - Xs[first:last, ] * beta_scaled[obs[1:m]]';
    matrix[m, n_rows[p]] W = mdivide_left_tri_low(L_sub, resid');
    target += -0.5 * (n_rows[p]
                      * (m * log(2 * pi()) + 2 * sum(log(diagonal(L_sub))))
                      + sum(square(W)));
  }

  for (k in 1:K)
    for (j in 1:J)
      beta_scaled[k, j] ~ normal(0, beta_prior_sd[j]);
  s_P ~ normal(0, 1);
  to_vector(loadings_raw) ~ normal(0, 1);
  to_vector(local_scale) ~ cauchy(0, 1);
  global_scale ~ cauchy(0, hs_global_scale);
}
generated quantities {
  matrix[K, M] loadings = loadings_tp;
  vector[K] uniqueness = uniqueness_tp;
  vector[K] communality = communality_tp;
  vector[K] sigma_P = y_sd .* s_P;
  matrix[K, K] P_std = loadings_tp * loadings_tp' + diag_matrix(uniqueness_tp);
  matrix[K, K] P = quad_form_diag(P_std, y_sd .* s_P);
  matrix[K, K] corrP = P_std;
  vector[M] factor_magnitude;
  for (m in 1:M)
    factor_magnitude[m] = dot_product(loadings_tp[, m], loadings_tp[, m]);

  // Pointwise log likelihood, one term per row of the response as the caller
  // passed it: the reordered blocks are mapped back through row_order.
  vector[N] log_lik;
  {
    for (p in 1:n_pattern) {
      int m = 0;
      array[K] int obs;
      for (k in 1:K)
        if (pattern_obs[p, k] == 1) {
          m += 1;
          obs[m] = k;
        }
      matrix[m, m] P_sub = P_std[obs[1:m], obs[1:m]];
      for (i in 1:m)
        for (j in 1:m)
          P_sub[i, j] *= s_P[obs[i]] * s_P[obs[j]];
      for (i in 1:m)
        P_sub[i, i] = fmax(P_sub[i, i], 1e-12);
      matrix[m, m] L_sub = cholesky_decompose(P_sub);
      int first = first_row[p];
      int last = first + n_rows[p] - 1;
      matrix[n_rows[p], m] resid = Ys[first:last, obs[1:m]]
                                   - Xs[first:last, ] * beta_scaled[obs[1:m]]';
      matrix[m, n_rows[p]] W = mdivide_left_tri_low(L_sub, resid');
      real constant = -0.5 * (m * log(2 * pi())
                              + 2 * sum(log(diagonal(L_sub))));
      for (n in 1:n_rows[p])
        log_lik[row_order[first + n - 1]] = constant - 0.5 * sum(square(col(W, n)));
    }
  }
}
