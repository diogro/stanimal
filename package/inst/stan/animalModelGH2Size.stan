// Multivariate animal model with a size-factor prior on the correlations
// (lmm_animal(correlation_prior = "size" or "shared_size")).
//
// This is animalModelGH2.stan with one change: each correlation matrix is one
// dense factor plus an LKJ remainder, normalised by row,
//   Omega = D (l l' + R) D,   D = diag(1 / sqrt(1 + l_k^2)),   R ~ LKJ(lkj_prior),
//   l_k = mu + sigma z_k,   z_k ~ N(0, 1),   sigma ~ half-normal(0, 1),
// with the mean loading written through rho = mu^2 / (1 + mu^2) ~ beta(1, 1),
// the weight of the size axis when every loading equals mu. rho = 0 with
// sigma = 0 is the LKJ model of animalModelGH2.stan. Unlike the horseshoe of
// the factor model, the loadings are pooled around a common mean, because a
// size factor is expected to load on every trait. R is a full correlation
// matrix, so Omega can be any correlation matrix: the factor moves where the
// prior shrinks to, it does not lower the rank.
//
// Why: LKJ shrinks every correlation toward zero, and for traits that share a
// size axis the truth is far from zero. In simulations of size-dominated
// traits the LKJ fits halved the genetic correlations and the heritabilities;
// this prior lets the data set where the correlations are pulled to.
//
// G and E get their own rho, sigma and LKJ remainder. With shared_pattern = 1
// they also share z, so the genetic and residual size axes deviate from equal
// loadings in the same direction, with their own strength: a weak
// shared-development assumption about the direction of size only. On the
// opossum skull data it improved leave-one-out prediction over independent
// patterns. It needs an unstructured E; with a diagonal E there is no residual
// size axis and the flag is refused in R.
//
// Anything that is not the correlation prior is shared with animalModelGH2.stan
// (and animalModelG.stan), so a change to the likelihood, the missing-data
// handling or the log_lik block has to be made in all of them.
//
functions {
  // Cholesky factor of D (l l' + L L') D, D = diag(1 / sqrt(1 + l^2)), from a
  // thin QR of B' with B = D [l, L]: Omega = B B' = R' R.
  matrix mixed_cholesky(vector l, matrix L) {
    int K = rows(L);
    matrix[K, K + 1] B = diag_pre_multiply(inv_sqrt(1 + square(l)),
                                           append_col(l, L));
    matrix[K, K] R = qr_thin_R(B');
    for (k in 1:K)
      if (R[k, k] < 0)
        R[k] = -R[k];
    return R';
  }
  // The rotated response with the gaps filled in.
  //
  // A missing value for individual n in trait k is not confined to record n of
  // the rotation: the rotated responses are Q' y, so it shifts every record, by
  // Q[n, i] * value in column k. That is also why incomplete traits cannot be
  // integrated out pattern by pattern here the way the phenotypic model does:
  // with relatedness a missing value is informed by every relative's observed
  // traits, so the marginal does not factorize over rows. Declaring the gaps as
  // parameters with no prior, and adding their contribution here, gives exactly
  // the observed-data marginal likelihood.
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
  matrix[N, J] X; // fixed effects design matrix
  matrix[N, K] Y; // response variables, with missing entries set to zero
  array[N, K] int<lower=0, upper=1> miss; // 1 where that trait is missing
  int<lower=0> n_missing; // number of gaps; each is one parameter
  array[n_missing] int<lower=1, upper=N> mis_row; // row of each gap
  array[n_missing] int<lower=1, upper=K> mis_col; // trait of each gap
  matrix[N, N] Q; // eigenvectors of the relationship matrix A, as columns
  vector<lower=0>[N] lambda; // matching eigenvalues of A, floored at 0
  real<lower=0> lkj_prior; // LKJ shape for the correlation matrices
  vector[J] x_mean_weight; // X * x_mean_weight is a column of ones
  int<lower=0, upper=1> unstructured_E; // 1: full residual covariance, 0: diagonal
  vector<lower=0>[2] beta_prior; // beta prior on the trait heritabilities
  int<lower=0, upper=1> shared_pattern; // 1: G and E share the size pattern z
}
transformed data {
  int E_dim = unstructured_E ? K : 0;
  int z_E_dim = shared_pattern ? 0 : E_dim; // E's own pattern, if not shared
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
  vector<lower=0, upper=1>[K] h2; // trait heritability, standardized scale
  vector<lower=0>[K] s_P; // phenotypic sds, standardized response scale
  vector[n_missing] Y_mis; // the gaps, on the standardized response scale
  cholesky_factor_corr[K] L_R_G; // LKJ part of the genetic correlations
  cholesky_factor_corr[E_dim] L_R_E; // LKJ part of the residual correlations
  real<lower=0, upper=1> rho_G; // mean loading of the size axis, as a weight
  array[unstructured_E] real<lower=0, upper=1> rho_E; // same for Omega_E
  real<lower=0> sigma_size_G; // spread of the size loadings around their mean
  array[unstructured_E] real<lower=0> sigma_size_E;
  vector[K] z_size_G; // standardised size loadings (also E's when shared)
  vector[z_E_dim] z_size_E;
}
transformed parameters {
  // The two scale vectors, rebuilt from the parameters so that every later
  // expression keeps its meaning. They are non-negative by construction, so
  // they carry no lower bound of their own.
  vector[K] s_G = sqrt(h2) .* s_P; // genetic sds
  vector[K] s_E = sqrt(1 - h2) .* s_P; // residual sds
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
  // Cholesky factors of the mixed correlation matrices, without forming them.
  // Omega = B B' with B = D [l, L_R], a K x (K + 1) matrix, so with B' = Q R
  // (thin QR) Omega = R' R and R' is a Cholesky factor once its rows are
  // signed to give a positive diagonal. Factorising the formed matrix instead
  // fails on rounding whenever Omega is close to singular, and every failure is
  // a rejected proposal.
  vector[K] l_G = sqrt(rho_G / (1 - rho_G)) + sigma_size_G * z_size_G;
  vector[E_dim] l_E;
  matrix[K, K] L_Omega_G = mixed_cholesky(l_G, L_R_G);
  matrix[E_dim, E_dim] L_Omega_E;
  if (unstructured_E) {
    l_E = sqrt(rho_E[1] / (1 - rho_E[1]))
          + sigma_size_E[1] * (shared_pattern ? z_size_G : z_size_E);
    L_Omega_E = mixed_cholesky(l_E, L_R_E);
  }
}
model {
  matrix[N, K] Ys = rotated_response(Ys_observed, Q, Y_mis, mis_row, mis_col);

  // Cholesky factors and covariances of G and E.
  matrix[K, K] L_G = diag_pre_multiply(s_G, L_Omega_G);
  matrix[K, K] L_E;
  if (unstructured_E)
    L_E = diag_pre_multiply(s_E, L_Omega_E);
  else
    L_E = diag_matrix(s_E);

  // log det E uses the true diagonal. The diagonal is floored only for the
  // whitening below, so that an extreme HMC proposal (e.g. a residual SD that
  // underflows to zero) yields a finite, hugely negative log density instead
  // of NaN/Inf, which would abort the chain.
  vector[K] log_diag_E = log(diagonal(L_E));
  for (k in 1:K)
    L_E[k, k] = fmax(L_E[k, k], 1e-12);

  matrix[K, K] G = L_G * L_G';
  matrix[K, K] E = L_E * L_E';

  // C = L_E^{-1} L_G is lower triangular, so C C' = L_E^{-1} G L_E^{-T}.
  matrix[K, K] C = mdivide_left_tri_low(L_E, L_G);
  matrix[K, K] M = C * C';
  M = (M + M') / 2; // symmetrize against rounding for the eigendecomposition
  for (i in 1:K)
    for (j in 1:K)
      M[i, j] = fmin(fmax(M[i, j], -1e150), 1e150); // guard overflow in the tails
  vector[K] d = eigenvalues_sym(M);
  for (k in 1:K)
    d[k] = fmax(d[k], 0); // M is PSD; guard rounding before log()

  // Eigenvector gradients are only well defined when the eigenvalues are well
  // separated. Near-coincident eigenvalues (the identity starting point, or G
  // proportional to E) fall back to the per-individual Cholesky route.
  real gap = positive_infinity();
  for (k in 2:K)
    gap = fmin(gap, d[k] - d[k - 1]);

  if (gap > 1e-4 * fmax(d[K], 1.0)) {
    // Fast path: rotate once into the eigenvectors of E^{-1/2} G E^{-1/2}.
    // Column n of W is L_E^{-1} r_n, and Z[n, k] = u_k' L_E^{-1} r_n.
    matrix[K, K] U = eigenvectors_sym(M);
    matrix[N, K] resid = Ys - Xs * beta_scaled';
    matrix[K, N] W = mdivide_left_tri_low(L_E, resid');
    matrix[N, K] Z = W' * U;

    // denom[n, k] = 1 + lambda_n d_k
    matrix[N, K] denom = rep_matrix(1, N, K) + lambda * d';
    target += -0.5 * (N * K * log(2 * pi())
                      + 2 * N * sum(log_diag_E)
                      + sum(log(denom))
                      + sum(square(Z) ./ denom));
  } else {
    // Robust fallback: per-individual Cholesky of lambda_n G + E.
    matrix[N, K] mu = Xs * beta_scaled';
    for (n in 1:N)
      target += multi_normal_cholesky_lpdf(Ys[n] | mu[n],
                                           cholesky_decompose(lambda[n] * G + E));
  }

  for (k in 1:K)
    for (j in 1:J)
      beta_scaled[k, j] ~ normal(0, beta_prior_sd[j]);
  s_P ~ normal(0, 1);
  h2 ~ beta(beta_prior[1], beta_prior[2]);
  L_R_G ~ lkj_corr_cholesky(lkj_prior);
  rho_G ~ beta(1, 1);
  sigma_size_G ~ normal(0, 1);
  z_size_G ~ std_normal();
  if (unstructured_E) {
    L_R_E ~ lkj_corr_cholesky(lkj_prior);
    rho_E ~ beta(1, 1);
    sigma_size_E ~ normal(0, 1);
    z_size_E ~ std_normal();
  }
}
generated quantities {
  matrix[K, K] G; // additive genetic covariance, original response scale
  matrix[K, K] E; // residual covariance, original response scale
  matrix[K, K] P; // phenotypic covariance
  corr_matrix[K] corrG;
  vector[K] sigma_G = y_sd .* s_G;
  vector[K] sigma_E = y_sd .* s_E;
  vector[N] log_lik;
  vector[K] size_loadings_G = l_G ./ sqrt(1 + square(l_G)); // correlation scale
  vector[E_dim] size_loadings_E = l_E ./ sqrt(1 + square(l_E));
  real integration_G; // relative eigenvalue variance of corrG
  real integration_E; // the same for the residual correlations
  matrix[K, K] G_std = multiply_lower_tri_self_transpose(diag_pre_multiply(s_G, L_Omega_G));
  matrix[K, K] E_std;
  if (unstructured_E)
    E_std = multiply_lower_tri_self_transpose(diag_pre_multiply(s_E, L_Omega_E));
  else
    E_std = diag_matrix(square(s_E));
  G = quad_form_diag(G_std, y_sd);
  E = quad_form_diag(E_std, y_sd);
  P = G + E;
  corrG = multiply_lower_tri_self_transpose(L_Omega_G);
  {
    // evolqg::CalcEigenVar(rel = TRUE) on a correlation matrix: the eigenvalue
    // variance over its maximum, (K - 1) for a K x K correlation matrix.
    vector[K] ev = eigenvalues_sym(corrG);
    integration_G = K > 1 ? (sum(square(ev - 1)) / K) / (K - 1) : 0;
    matrix[K, K] corrE = unstructured_E
      ? multiply_lower_tri_self_transpose(L_Omega_E) : diag_matrix(rep_vector(1, K));
    ev = eigenvalues_sym(corrE);
    integration_E = K > 1 ? (sum(square(ev - 1)) / K) / (K - 1) : 0;
  }

  {
    // The rotated response belongs here, not to the saved quantities: it is a
    // derived matrix the size of the data, and saving it would bloat every draw.
    matrix[N, K] Ys = rotated_response(Ys_observed, Q, Y_mis, mis_row, mis_col);
    matrix[K, K] L_G = diag_pre_multiply(s_G, L_Omega_G);
    matrix[K, K] L_E;
    if (unstructured_E)
      L_E = diag_pre_multiply(s_E, L_Omega_E);
    else
      L_E = diag_matrix(s_E);
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
