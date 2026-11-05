// Multivariate linear model for the phenotypic (residual) covariance matrix P,
// for data with no relatedness or genomic information. Parametrized like the
// animal model: scales times an LKJ-distributed correlation matrix.
//
// Incomplete traits: without relatedness the records are independent, so an
// incomplete record's marginal density is the normal on its observed subvector,
// with the matching principal submatrix of P. Rows are grouped by which traits
// they observe (a pattern), and each pattern contributes one Cholesky of its own
// submatrix. The missing values are integrated out analytically, so they cost no
// parameters, and a complete response is the single-pattern case of the same
// code. The rows are reordered inside the model so that each pattern is a
// contiguous block; every returned quantity is mapped back to the caller's order.
data {
  int<lower=1> K; // number of traits
  int<lower=1> J; // number of fixed effects
  int<lower=2> N; // number of individuals
  matrix[N, J] X; // fixed effects design matrix
  matrix[N, K] Y; // response variables, with missing entries set to zero
  array[N, K] int<lower=0, upper=1> miss; // 1 where that trait is missing
  int<lower=1, upper=N> n_pattern; // number of distinct missingness patterns
  array[n_pattern, K] int<lower=0, upper=1> pattern_obs; // observed traits per pattern
  array[N] int<lower=1, upper=N> row_order; // rows grouped by pattern
  array[n_pattern] int<lower=0> n_rows; // rows in each pattern
  array[n_pattern] int<lower=1> first_row; // first reordered row of each pattern
  real<lower=0> lkj_prior; // LKJ shape for the correlation matrix
  vector[J] x_mean_weight; // X * x_mean_weight is a column of ones
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
      // Centre and scale on the observed values only.
      y_mean[k] = total / present;
      real sum_sq = 0;
      for (n in 1:N)
        if (miss[row_order[n], k] == 0)
          sum_sq += square(Yp[n, k] - y_mean[k]);
      y_sd[k] = sqrt(sum_sq / (present - 1));
      // A missing entry stays at zero. The likelihood reads only the observed
      // entries of each pattern, so the placeholder is never used as a value.
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
  cholesky_factor_corr[K] L_Omega_P;
}
transformed parameters {
  matrix[K, J] beta; // fixed effects on original predictor and response scales
  for (k in 1:K) {
    for (j in 1:J) {
      beta[k, j] = y_sd[k] * beta_scaled[k, j] * x_scale_inv[j];
    }
  }
  // The responses are centred and scaled, so the model predicts deviations from
  // their mean. Put the mean back where the design can hold it: X
  // x_mean_weight is a column of ones, so the shift lands on every prediction,
  // whether the design has an intercept column or the indicators of a factor.
  for (k in 1:K) {
    real mean_shift = y_mean[k] - dot_product(x_center, beta[k]');
    for (j in 1:J)
      beta[k, j] += mean_shift * x_mean_weight[j];
  }
}
model {
  // Cholesky factor of P on the standardized response scale.
  matrix[K, K] L_P = diag_pre_multiply(s_P, L_Omega_P);
  matrix[K, K] P_std = L_P * L_P';

  for (p in 1:n_pattern) {
    // Traits this pattern observes.
    int m = 0;
    array[K] int obs;
    for (k in 1:K)
      if (pattern_obs[p, k] == 1) {
        m += 1;
        obs[m] = k;
      }
    // Principal submatrix of P over those traits. The diagonal is floored so
    // that an extreme proposal (a scale underflowing to zero) gives a finite,
    // hugely negative log density instead of an exception out of the Cholesky.
    // A floored diagonal keeps the submatrix positive semidefinite.
    matrix[m, m] P_sub = P_std[obs[1:m], obs[1:m]];
    for (i in 1:m)
      P_sub[i, i] = fmax(P_sub[i, i], 1e-12);
    matrix[m, m] L_sub = cholesky_decompose(P_sub);

    int first = first_row[p];
    int last = first + n_rows[p] - 1;
    // Whitened residuals: column n of W is L_sub^{-1} r_n over the observed
    // traits, so sum(square(W)) is the quadratic form of the density.
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
  L_Omega_P ~ lkj_corr_cholesky(lkj_prior);
}
generated quantities {
  vector[K] sigma_P = y_sd .* s_P; // phenotypic sds, original response scale
  matrix[K, K] P = quad_form_diag(
    multiply_lower_tri_self_transpose(diag_pre_multiply(s_P, L_Omega_P)), y_sd);
  corr_matrix[K] corrP = multiply_lower_tri_self_transpose(L_Omega_P);

  // Pointwise log likelihood, one term per row of the response as the caller
  // passed it: the rotated order is mapped back through row_order.
  vector[N] log_lik;
  {
    matrix[K, K] P_std = multiply_lower_tri_self_transpose(
      diag_pre_multiply(s_P, L_Omega_P));
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
