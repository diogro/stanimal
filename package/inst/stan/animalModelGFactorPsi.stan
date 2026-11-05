// Multivariate animal model with a sparse factor structure for G, with the
// trait-specific genetic variance as a free parameter.
//
// G = Lambda Lambda' + diag(psi) on the standardized response scale, with
// regularized-horseshoe priors on the loadings and half-normal priors on the
// specific genetic standard deviations psi and on the residual standard
// deviations s_E. E remains unstructured or diagonal as in animalModelG.stan.
//
// This differs from animalModelGFactor.stan in where the heritability comes
// from. There the loadings were normalized to a correlation matrix, so the
// trait-specific genetic variance was 1 / (1 + ||raw||^2), a deterministic
// function of the loadings, and the heritability had to be a parameter of its
// own. That made h2 and the communality inseparable: the data only fix the split
// of the phenotypic covariance between G and E, so a prior on the loadings moved
// the heritability. Here psi carries the specific genetic variance, so the
// loadings only decide how much of G is shared and h2 is reported as
// diag(G) / diag(P) rather than imposed. The parameterization of G follows
// Runcie and Mukherjee (2013), with unstructured residuals instead of their
// diagonal ones.
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
  int<lower=1, upper=K - 1> M; // number of genetic factors
  matrix[N, J] X; // fixed effects design matrix
  matrix[N, K] Y; // response variables, with missing entries set to zero
  array[N, K] int<lower=0, upper=1> miss; // 1 where that trait is missing
  int<lower=0> n_missing; // number of gaps; each is one parameter
  array[n_missing] int<lower=1, upper=N> mis_row; // row of each gap
  array[n_missing] int<lower=1, upper=K> mis_col; // trait of each gap
  matrix[N, N] A; // known relationship matrix (symmetric, PSD)
  real<lower=0> lkj_prior; // LKJ shape for the residual correlation
  vector[J] x_mean_weight; // X * x_mean_weight is a column of ones
  int<lower=0, upper=1> unstructured_E; // 1: full residual covariance, 0: diagonal
  real<lower=0> hs_global_scale; // half-Cauchy scale for the global shrinkage
  real<lower=0> hs_slab_scale; // slab scale of the regularized horseshoe
  real<lower=0, upper=1> column_decay; // sd multiplier for successive factors
  vector<lower=0>[2] beta_prior; // unused here, kept so the data list matches
}
transformed data {
  int E_dim = unstructured_E ? K : 0;
  vector[N] lambda;
  matrix[N, N] Q = eigenvectors_sym(A);
  vector[M] column_scale;
  for (m in 1:M)
    column_scale[m] = pow(column_decay, m - 1);
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

  {
    vector[N] ev = eigenvalues_sym(A);
    for (n in 1:N)
      lambda[n] = fmax(ev[n], 0);
  }

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
  vector<lower=0>[K] psi_sd; // specific genetic sds, standardized scale
  vector<lower=0>[K] s_E; // residual sds, standardized response scale
  vector[n_missing] Y_mis; // the gaps, on the standardized response scale
  cholesky_factor_corr[E_dim] L_Omega_E;
  matrix[K, M] loadings_raw; // non-centred loadings
  matrix<lower=0>[K, M] local_scale; // local horseshoe scales
  real<lower=0> global_scale; // global horseshoe scale
}
transformed parameters {
  vector[K] psi = square(psi_sd); // specific genetic variances
  matrix[K, J] beta; // fixed effects on original predictor and response scales
  matrix[K, M] loadings; // loadings, standardized response scale
  vector[K] communality; // share of G explained by the factors
  vector[K] uniqueness; // the remainder, psi / (shared + psi)
  matrix[K, K] G_std; // genetic covariance on the standardized scale
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
  // Regularized horseshoe (Piironen & Vehtari 2017) with an optional
  // deterministic column decay: factor m has standard-deviation multiplier
  // column_decay^(m - 1). There is no row normalization here: the loadings are
  // on the response scale and psi absorbs what the factors do not explain.
  for (k in 1:K) {
    for (m in 1:M) {
      // Clamp only to keep an overflowing Cauchy proposal from producing
      // inf / inf; the horseshoe shape is unchanged for finite values.
      real local_variance = fmin(
        square(global_scale) * square(local_scale[k, m]) *
        square(column_scale[m]), 1e150);
      real slab_variance = square(hs_slab_scale);
      real shrunk_variance = slab_variance * local_variance /
        (slab_variance + local_variance);
      loadings[k, m] = loadings_raw[k, m] * sqrt(shrunk_variance);
    }
  }
  G_std = loadings * loadings' + diag_matrix(psi);
  {
    for (k in 1:K) {
      real shared = 0;
      for (m in 1:M)
        shared += square(loadings[k, m]);
      communality[k] = shared / (shared + psi[k]);
      uniqueness[k] = psi[k] / (shared + psi[k]);
    }
  }
}
model {
  matrix[N, K] Ys = rotated_response(Ys_observed, Q, Y_mis, mis_row, mis_col);

  // Cholesky factor of the residual covariance.
  matrix[K, K] L_E;
  if (unstructured_E)
    L_E = diag_pre_multiply(s_E, L_Omega_E);
  else
    L_E = diag_matrix(s_E);

  // The genetic covariance on the standardized scale, built in the transformed
  // parameters block from the loadings and the specific genetic variances.
  matrix[K, K] G = G_std;
  matrix[K, K] L_G = cholesky_decompose(G);

  // log det E uses the true diagonal. The diagonal is floored only for the
  // whitening below, so that an extreme HMC proposal yields a finite, hugely
  // negative log density instead of NaN/Inf.
  vector[K] log_diag_E = log(diagonal(L_E));
  for (k in 1:K)
    L_E[k, k] = fmax(L_E[k, k], 1e-12);

  matrix[K, K] E = L_E * L_E';

  // C = L_E^{-1} L_G is lower triangular, so C C' = L_E^{-1} G L_E^{-T}.
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
                                           cholesky_decompose(lambda[n] * G + E));
  }

  for (k in 1:K)
    for (j in 1:J)
      beta_scaled[k, j] ~ normal(0, beta_prior_sd[j]);
  psi_sd ~ normal(0, 1);
  s_E ~ normal(0, 1);
  if (unstructured_E)
    L_Omega_E ~ lkj_corr_cholesky(lkj_prior);
  to_vector(loadings_raw) ~ normal(0, 1);
  to_vector(local_scale) ~ cauchy(0, 1);
  global_scale ~ cauchy(0, hs_global_scale);
}
generated quantities {
  // loadings, communality and uniqueness are transformed parameters and are
  // saved already; what follows converts everything to the data scale and
  // reports the heritability, which is derived here rather than imposed.
  vector[M] factor_magnitude;
  for (m in 1:M)
    factor_magnitude[m] = dot_product(loadings[, m], loadings[, m]);
  matrix[K, K] E_std;
  if (unstructured_E)
    E_std = multiply_lower_tri_self_transpose(
      diag_pre_multiply(s_E, L_Omega_E));
  else
    E_std = diag_matrix(square(s_E));
  matrix[K, K] G = quad_form_diag(G_std, y_sd);
  matrix[K, K] E = quad_form_diag(E_std, y_sd);
  matrix[K, K] P = G + E;
  vector[K] sigma_G = y_sd .* sqrt(diagonal(G_std));
  vector[K] sigma_E = sqrt(diagonal(E));
  vector[K] h2;
  matrix[K, K] corrG;
  for (k in 1:K)
    h2[k] = G[k, k] / P[k, k];
  for (i in 1:K)
    for (j in 1:K)
      corrG[i, j] = G_std[i, j] / sqrt(G_std[i, i] * G_std[j, j]);
  vector[N] log_lik;

  {
    // The rotated response belongs here, not to the saved quantities: it is a
    // derived matrix the size of the data, and saving it would bloat every draw.
    matrix[N, K] Ys = rotated_response(Ys_observed, Q, Y_mis, mis_row, mis_col);
    matrix[K, K] L_E;
    if (unstructured_E)
      L_E = diag_pre_multiply(s_E, L_Omega_E);
    else
      L_E = diag_matrix(s_E);
    vector[K] log_diag_E = log(diagonal(L_E));
    for (k in 1:K)
      L_E[k, k] = fmax(L_E[k, k], 1e-12);
    matrix[K, K] L_G = cholesky_decompose(G_std);
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
