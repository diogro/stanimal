// Multivariate animal model with shared factors for G and E.
// G = Lambda diag(h2_factor) Lambda' + diag(psi_G)
// E = Lambda diag(1 - h2_factor) Lambda' + diag(psi_E)
// The factor variances sum to one, so h2_factor is each factor's heritability.
// E is therefore low-rank plus diagonal, and the same factors explain genetic
// and residual covariance.
functions {
  // The rotated response with the gaps filled in. See animalModelG.stan for why
  // incomplete traits need parameters here instead of per-pattern marginals: a
  // missing value in individual n shifts every rotated record by Q[n, i].
  matrix rotated_response(matrix Ys_observed, matrix Q, vector Y_mis,
                          array[] int mis_row, array[] int mis_col) {
    matrix[rows(Ys_observed), cols(Ys_observed)] Ys = Ys_observed;
    for (j in 1:size(Y_mis))
      Ys[, mis_col[j]] += Q[mis_row[j]]' * Y_mis[j];
    return Ys;
  }
}
data {
  int<lower=1> K; // number of traits
  int<lower=1> J; // number of fixed effects
  int<lower=2> N; // number of individuals
  int<lower=1, upper=K - 1> M; // number of shared factors
  matrix[N, J] X; // fixed effects design matrix
  matrix[N, K] Y; // response variables, with missing entries set to zero
  array[N, K] int<lower=0, upper=1> miss; // 1 where that trait is missing
  int<lower=0> n_missing; // number of gaps; each is one parameter
  array[n_missing] int<lower=1, upper=N> mis_row; // row of each gap
  array[n_missing] int<lower=1, upper=K> mis_col; // trait of each gap
  matrix[N, N] Q; // eigenvectors of the relationship matrix A, as columns
  vector<lower=0>[N] lambda; // matching eigenvalues of A, floored at 0
  vector[J] x_mean_weight; // X * x_mean_weight is a column of ones
  real<lower=0> hs_global_scale; // half-Cauchy scale for the global shrinkage
  real<lower=0> hs_slab_scale; // slab scale of the regularized horseshoe
  real<lower=0, upper=1> column_decay; // sd multiplier for successive factors
  real<lower=0> h2_prior_a; // beta prior shape for factor heritabilities
  real<lower=0> h2_prior_b; // beta prior shape for factor heritabilities
}
transformed data {
  vector[M] column_scale;
  vector[K] y_sd;
  vector[K] y_mean;
  vector[J] x_center;
  vector[J] x_scale;
  vector[J] x_scale_inv;
  vector[J] beta_prior_sd;
  array[J] int x_is_constant;
  vector[J] x_constant_value;
  matrix[N, K] Ys_observed; // rotated standardized responses, gaps at zero
  matrix[N, J] Xs;

  for (m in 1:M)
    column_scale[m] = pow(column_decay, m - 1);

  {
    matrix[N, K] Y_std;
    matrix[N, J] X_std;
    for (k in 1:K) {
      int present = 0;
      real total = 0;
      for (n in 1:N)
        if (miss[n, k] == 0) {
          present += 1;
          total += Y[n, k];
        }
      // Centre and scale on the observed values only.
      y_mean[k] = total / present;
      real sum_sq = 0;
      for (n in 1:N)
        if (miss[n, k] == 0)
          sum_sq += square(Y[n, k] - y_mean[k]);
      y_sd[k] = sqrt(sum_sq / (present - 1));
      // A gap is a parameter. Here it only needs a placeholder, because the
      // rotation below is linear and its contribution is added back per draw.
      for (n in 1:N)
        Y_std[n, k] = miss[n, k] == 1 ? 0 : (Y[n, k] - y_mean[k]) / y_sd[k];
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
      X_std[, j] = (X[, j] - x_center[j]) * x_scale_inv[j];
    }
    Ys_observed = Q' * Y_std;
    Xs = Q' * X_std;
  }
}
parameters {
  matrix[K, J] beta_scaled; // fixed effects on scaled predictors
  vector<lower=0>[K] sqrt_psi_G; // trait-specific genetic sds
  vector[n_missing] Y_mis; // the gaps, on the standardized response scale
  vector<lower=0>[K] sqrt_psi_E; // trait-specific residual sds
  vector<lower=0, upper=1>[M] h2_factor; // heritability of each factor
  matrix[K, M] loadings_raw; // non-centred loadings
  matrix<lower=0>[K, M] local_scale; // local horseshoe scales
  real<lower=0> global_scale; // global horseshoe scale
}
transformed parameters {
  matrix[K, J] beta; // fixed effects on original predictor and response scales
  matrix[K, M] loadings_tp; // shared factor loadings
  vector<lower=0>[K] psi_G; // trait-specific genetic variances
  vector<lower=0>[K] psi_E; // trait-specific residual variances
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
  for (k in 1:K) {
    psi_G[k] = square(sqrt_psi_G[k]);
    psi_E[k] = square(sqrt_psi_E[k]);
  }
  // Regularized horseshoe with a deterministic column decay. The same
  // loadings enter G and E, so only one loading matrix is estimated.
  for (k in 1:K) {
    for (m in 1:M) {
      real local_variance = fmin(
        square(global_scale) * square(local_scale[k, m]) *
        square(column_scale[m]), 1e150);
      real slab_variance = square(hs_slab_scale);
      real shrunk_variance = slab_variance * local_variance /
        (slab_variance + local_variance);
      loadings_tp[k, m] = loadings_raw[k, m] * sqrt(shrunk_variance);
    }
  }
}
model {
  matrix[N, K] Ys = rotated_response(Ys_observed, Q, Y_mis, mis_row, mis_col);
  vector[M] h2_f = h2_factor;
  vector[M] h2_c = 1 - h2_factor;
  matrix[K, K] G_std = loadings_tp * diag_matrix(h2_f) * loadings_tp' +
    diag_matrix(psi_G);
  matrix[K, K] E_std = loadings_tp * diag_matrix(h2_c) * loadings_tp' +
    diag_matrix(psi_E);
  matrix[K, K] L_G = cholesky_decompose(G_std);
  matrix[K, K] L_E = cholesky_decompose(E_std);

  // log det E uses the true diagonal; the floor only keeps an extreme HMC
  // proposal from producing NaN/Inf in the whitening.
  vector[K] log_diag_E = log(diagonal(L_E));
  for (k in 1:K)
    L_E[k, k] = fmax(L_E[k, k], 1e-12);

  matrix[K, K] C = mdivide_left_tri_low(L_E, L_G);
  matrix[K, K] Mmat = C * C';
  Mmat = (Mmat + Mmat') / 2;
  for (i in 1:K)
    for (j in 1:K)
      Mmat[i, j] = fmin(fmax(Mmat[i, j], -1e150), 1e150);
  vector[K] d = eigenvalues_sym(Mmat);
  for (k in 1:K)
    d[k] = fmax(d[k], 0);

  real gap = positive_infinity();
  for (k in 2:K)
    gap = fmin(gap, d[k] - d[k - 1]);

  if (gap > 1e-4 * fmax(d[K], 1.0)) {
    matrix[K, K] U = eigenvectors_sym(Mmat);
    matrix[N, K] resid = Ys - Xs * beta_scaled';
    matrix[K, N] W = mdivide_left_tri_low(L_E, resid');
    matrix[N, K] Z = W' * U;
    matrix[N, K] denom = rep_matrix(1, N, K) + lambda * d';
    target += -0.5 * (N * K * log(2 * pi())
                      + 2 * N * sum(log_diag_E)
                      + sum(log(denom))
                      + sum(square(Z) ./ denom));
  } else {
    matrix[N, K] mu = Xs * beta_scaled';
    for (n in 1:N)
      target += multi_normal_cholesky_lpdf(Ys[n] | mu[n],
                                           cholesky_decompose(lambda[n] * G_std + E_std));
  }

  for (k in 1:K)
    for (j in 1:J)
      beta_scaled[k, j] ~ normal(0, beta_prior_sd[j]);
  sqrt_psi_G ~ normal(0, 1);
  sqrt_psi_E ~ normal(0, 1);
  h2_factor ~ beta(h2_prior_a, h2_prior_b);
  to_vector(loadings_raw) ~ normal(0, 1);
  to_vector(local_scale) ~ cauchy(0, 1);
  global_scale ~ cauchy(0, hs_global_scale);
}
generated quantities {
  matrix[K, M] loadings = loadings_tp;
  vector[K] psi_G_out = psi_G;
  vector[K] psi_E_out = psi_E;
  vector[K] genetic_uniqueness = psi_G;
  vector[K] residual_uniqueness = psi_E;
  vector[K] uniqueness = psi_G + psi_E;
  vector[K] communality;
  vector[M] factor_magnitude;
  vector[M] factor_heritability = h2_factor;
  vector[M] h2_f = h2_factor;
  vector[M] h2_c = 1 - h2_factor;
  matrix[K, K] G_std = loadings_tp * diag_matrix(h2_f) * loadings_tp' +
    diag_matrix(psi_G);
  matrix[K, K] E_std = loadings_tp * diag_matrix(h2_c) * loadings_tp' +
    diag_matrix(psi_E);
  matrix[K, K] G = quad_form_diag(G_std, y_sd);
  matrix[K, K] E = quad_form_diag(E_std, y_sd);
  matrix[K, K] P = G + E;
  vector[K] sigma_G = sqrt(diagonal(G));
  vector[K] sigma_E = sqrt(diagonal(E));
  vector[K] h2;
  vector[K] inv_sd_G = inv_sqrt(diagonal(G_std));
  matrix[K, K] corrG = diag_post_multiply(
    diag_pre_multiply(inv_sd_G, G_std), inv_sd_G);
  vector[N] log_lik;

  for (k in 1:K) {
    h2[k] = G[k, k] / P[k, k];
    real row_sq = 0;
    for (m in 1:M)
      row_sq += square(loadings[k, m]);
    communality[k] = row_sq;
  }
  for (m in 1:M)
    factor_magnitude[m] = dot_product(loadings[, m], loadings[, m]);

  {
    // The rotated response belongs here, not to the saved quantities: it is a
    // derived matrix the size of the data, and saving it would bloat every draw.
    matrix[N, K] Ys = rotated_response(Ys_observed, Q, Y_mis, mis_row, mis_col);
    matrix[K, K] L_E = cholesky_decompose(E_std);
    matrix[K, K] L_G = cholesky_decompose(G_std);
    vector[K] log_diag_E = log(diagonal(L_E));
    for (k in 1:K)
      L_E[k, k] = fmax(L_E[k, k], 1e-12);
    matrix[K, K] C = mdivide_left_tri_low(L_E, L_G);
    matrix[K, K] Mmat = C * C';
    Mmat = (Mmat + Mmat') / 2;
    for (i in 1:K)
      for (j in 1:K)
        Mmat[i, j] = fmin(fmax(Mmat[i, j], -1e150), 1e150);
    vector[K] d = eigenvalues_sym(Mmat);
    for (k in 1:K)
      d[k] = fmax(d[k], 0);
    real gap = positive_infinity();
    for (k in 2:K)
      gap = fmin(gap, d[k] - d[k - 1]);

    if (gap > 1e-4 * fmax(d[K], 1.0)) {
      matrix[K, K] U = eigenvectors_sym(Mmat);
      matrix[N, K] resid = Ys - Xs * beta_scaled';
      matrix[K, N] W = mdivide_left_tri_low(L_E, resid');
      matrix[N, K] Z = W' * U;
      for (n in 1:N) {
        real quad = 0;
        real log_denom = 0;
        for (k in 1:K) {
          real denom = 1 + lambda[n] * d[k];
          quad += square(Z[n, k]) / denom;
          log_denom += log(denom);
        }
        log_lik[n] = -0.5 * (K * log(2 * pi())
                             + 2 * sum(log_diag_E)
                             + log_denom + quad);
      }
    } else {
      matrix[N, K] mu = Xs * beta_scaled';
      for (n in 1:N)
        log_lik[n] = multi_normal_cholesky_lpdf(
          Ys[n] | mu[n],
          cholesky_decompose(lambda[n] * G_std + E_std));
    }
  }
}
