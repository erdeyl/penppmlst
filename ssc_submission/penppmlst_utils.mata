*! version 0.5.0  05jan2026
*! Utility functions for penppmlst package
*! Stata implementation by Erdey, László (2026)
*!   Faculty of Economics and Business, University of Debrecen, Hungary
*! Based on R penppml by Breinlich, Corradi, Rocha, Ruta, Santos Silva, Zylkin

version 17.0
mata:
mata set matastrict on

// ============================================================================
// SOFT-THRESHOLDING OPERATOR
// Core operation for lasso: S(z, gamma) = sign(z) * max(|z| - gamma, 0)
// ============================================================================

real colvector soft_threshold(real colvector z, real colvector gamma)
{
    real colvector result

    result = sign(z) :* max((abs(z) :- gamma, J(rows(z), 1, 0)))
    return(result)
}

// Scalar version for single coefficient updates
real scalar soft_threshold_scalar(real scalar z, real scalar gamma)
{
    if (z > gamma) return(z - gamma)
    if (z < -gamma) return(z + gamma)
    return(0)
}

// ============================================================================
// WEIGHTED RIDGE REGRESSION SOLVER
// Solves: min_beta sum_i w_i (y_i - x_i'beta)^2 + lambda * ||beta||_2^2
// Analytical solution: beta = (X'WX + lambda*I)^{-1} X'Wy
// ============================================================================

real colvector ridge_solve(real matrix X, real colvector y, real colvector w,
                           real scalar lambda)
{
    real matrix XWX
    real colvector XWy
    real scalar p
    real matrix XWX_ridge
    real colvector beta

    p = cols(X)

    // Compute weighted cross-products using quad precision for stability
    XWX = quadcross(X, w, X)
    XWy = quadcross(X, w, y)

    // Add ridge penalty to diagonal
    XWX_ridge = XWX + lambda * I(p)

    // Solve using Cholesky decomposition (X'WX + lambda*I is positive definite)
    beta = cholsolve(XWX_ridge, XWy)

    // Fallback to QR if Cholesky fails (near-singular case)
    if (hasmissing(beta)) {
        beta = qrsolve(XWX_ridge, XWy)
    }

    return(beta)
}

// ============================================================================
// WEIGHTED COORDINATE DESCENT FOR LASSO / ELASTIC NET
// Solves: min_beta 0.5 * sum_i w_i (y_i - x_i'beta)^2
//                + lambda * (alpha * ||beta||_1 + (1-alpha)/2 * ||beta||_2^2)
// Uses cyclic coordinate descent with soft-thresholding updates
// ============================================================================

real colvector coordinate_descent(real matrix X, real colvector y, real colvector w,
                                  real scalar lambda, real colvector psi,
                                  real scalar alpha,
                                  | real scalar tol, real scalar maxiter,
                                    real scalar verbose)
{
    real scalar n, p, j, iter, converged
    real colvector beta, beta_old, residual
    real colvector Xj, wXj
    real scalar Xj_norm_sq, partial_resid, threshold, update
    real scalar max_change
    real rowvector XX_diag
    real matrix wX

    // Set defaults
    if (args() < 7 | tol == .) tol = 1e-7
    if (args() < 8 | maxiter == .) maxiter = 1000
    if (args() < 9 | verbose == .) verbose = 0

    n = rows(X)
    p = cols(X)

    // Initialize beta with zeros (cold start) or ridge solution (warm start)
    if (alpha < 1) {
        // Warm start with ridge for elastic net
        beta = ridge_solve(X, y, w, lambda * (1 - alpha))
    } else {
        beta = J(p, 1, 0)
    }

    // Precompute weighted X and diagonal of X'WX for efficiency
    wX = X :* sqrt(w)
    XX_diag = colsum(wX :^ 2)'  // p x 1 vector of ||w^{1/2} X_j||^2

    // Initial residual
    residual = y - X * beta

    converged = 0
    for (iter = 1; iter <= maxiter; iter++) {
        beta_old = beta

        // Cycle through all coordinates
        for (j = 1; j <= p; j++) {
            Xj = X[., j]
            Xj_norm_sq = XX_diag[j]

            // Skip if column is zero (degenerate case)
            if (Xj_norm_sq < 1e-10) continue

            // Compute partial residual: r_j = y - X * beta + X_j * beta_j
            // Equivalently: X_j' * W * r_j = X_j' * W * (residual + X_j * beta_j)
            partial_resid = quadcross(Xj, w, residual) + Xj_norm_sq * beta[j]

            // Compute threshold for this coefficient
            threshold = lambda * alpha * psi[j]

            // Soft-thresholding update
            if (alpha > 0) {
                update = soft_threshold_scalar(partial_resid, threshold) /
                         (Xj_norm_sq + lambda * (1 - alpha))
            } else {
                // Pure ridge (no thresholding)
                update = partial_resid / (Xj_norm_sq + lambda)
            }

            // Update residual if coefficient changed
            if (update != beta[j]) {
                residual = residual - Xj * (update - beta[j])
                beta[j] = update
            }
        }

        // Check convergence: max absolute change in coefficients
        max_change = max(abs(beta - beta_old))

        if (verbose) {
            printf("Iteration %g: max change = %g\n", iter, max_change)
        }

        if (max_change < tol) {
            converged = 1
            break
        }
    }

    if (!converged & verbose) {
        printf("Warning: coordinate descent did not converge in %g iterations\n", maxiter)
    }

    return(beta)
}

// ============================================================================
// ACTIVE SET COORDINATE DESCENT (FASTER FOR SPARSE SOLUTIONS)
// Only updates non-zero coefficients after initial passes
// ============================================================================

real colvector coordinate_descent_active(real matrix X, real colvector y,
                                         real colvector w, real scalar lambda,
                                         real colvector psi, real scalar alpha,
                                         | real scalar tol, real scalar maxiter,
                                           real scalar verbose)
{
    real scalar n, p, j, iter, converged, full_pass
    real colvector beta, beta_old, residual
    real colvector Xj
    real scalar Xj_norm_sq, partial_resid, threshold, update
    real scalar max_change
    real rowvector XX_diag
    real colvector active_set
    real scalar n_active, idx

    // Set defaults
    if (args() < 7 | tol == .) tol = 1e-7
    if (args() < 8 | maxiter == .) maxiter = 1000
    if (args() < 9 | verbose == .) verbose = 0

    n = rows(X)
    p = cols(X)

    // Initialize
    beta = J(p, 1, 0)
    XX_diag = colsum((X :* sqrt(w)) :^ 2)'
    residual = y  // Initial residual is just y since beta = 0
    active_set = 1::p  // Start with all variables active

    converged = 0
    full_pass = 1  // First pass is always full

    for (iter = 1; iter <= maxiter; iter++) {
        beta_old = beta
        n_active = rows(active_set)

        // Cycle through active set
        for (idx = 1; idx <= n_active; idx++) {
            j = active_set[idx]
            Xj = X[., j]
            Xj_norm_sq = XX_diag[j]

            if (Xj_norm_sq < 1e-10) continue

            partial_resid = quadcross(Xj, w, residual) + Xj_norm_sq * beta[j]
            threshold = lambda * alpha * psi[j]

            if (alpha > 0) {
                update = soft_threshold_scalar(partial_resid, threshold) /
                         (Xj_norm_sq + lambda * (1 - alpha))
            } else {
                update = partial_resid / (Xj_norm_sq + lambda)
            }

            if (update != beta[j]) {
                residual = residual - Xj * (update - beta[j])
                beta[j] = update
            }
        }

        max_change = max(abs(beta - beta_old))

        if (max_change < tol) {
            if (full_pass) {
                // Converged on full pass - done
                converged = 1
                break
            } else {
                // Converged on active set - do a full pass to check
                full_pass = 1
                active_set = 1::p
            }
        } else {
            // Update active set to non-zero coefficients
            active_set = select(1::p, beta :!= 0)
            if (rows(active_set) == 0) active_set = 1::p  // Prevent empty
            full_pass = 0
        }
    }

    if (!converged & verbose) {
        printf("Warning: active set CD did not converge in %g iterations\n", maxiter)
    }

    return(beta)
}

// ============================================================================
// LAMBDA SEQUENCE GENERATION
// Creates a decreasing sequence of lambda values for regularization path
// ============================================================================

real colvector generate_lambda_sequence(real matrix X, real colvector y,
                                        real colvector w, real colvector psi,
                                        real scalar alpha,
                                        | real scalar nlambda,
                                          real scalar lambda_min_ratio)
{
    real scalar lambda_max, lambda_min
    real colvector lambdas
    real colvector gradient
    real scalar j, p

    // Set defaults
    if (args() < 6 | nlambda == .) nlambda = 100
    if (args() < 7 | lambda_min_ratio == .) lambda_min_ratio = 0.01

    p = cols(X)

    // Compute lambda_max: smallest lambda that sets all coefficients to zero
    // At beta = 0, the gradient of the unpenalized objective is X'Wy
    // Lambda_max = max_j |X_j'Wy| / (alpha * psi_j)
    gradient = abs(quadcross(X, w, y))

    lambda_max = 0
    for (j = 1; j <= p; j++) {
        if (psi[j] > 0 & alpha > 0) {
            lambda_max = max((lambda_max, gradient[j] / (alpha * psi[j])))
        }
    }

    // Small adjustment to ensure first lambda gives empty model
    lambda_max = lambda_max * 1.01

    // Lambda_min
    lambda_min = lambda_max * lambda_min_ratio

    // Generate log-spaced sequence
    lambdas = exp(rangen(ln(lambda_max), ln(lambda_min), nlambda))

    return(lambdas)
}

// ============================================================================
// POISSON DEVIANCE COMPUTATION
// Deviance = 2 * sum[ y*log(y/mu) - (y - mu) ]
// ============================================================================

real scalar compute_deviance(real colvector y, real colvector mu)
{
    real scalar n, dev
    real colvector y_pos, mu_pos
    real colvector log_term

    n = rows(y)

    // Handle y = 0 case: 0 * log(0) = 0 by convention
    y_pos = y :+ (y :== 0)  // Replace 0 with 1 for log computation
    mu_pos = max((mu, J(n, 1, 1e-10)))  // Prevent log(0)

    log_term = y :* ln(y_pos :/ mu_pos)
    log_term = log_term :* (y :> 0)  // Zero out where y = 0

    dev = 2 * sum(log_term - (y - mu))

    return(dev)
}

// ============================================================================
// WEIGHTED COLUMN STANDARDIZATION
// Centers and scales columns of X using weights w
// Returns: standardized X, means, standard deviations
// ============================================================================

void standardize_weighted(real matrix X, real colvector w,
                          real matrix X_std, real rowvector means,
                          real rowvector sds)
{
    real scalar n, p, j, w_sum
    real colvector wX_j
    real scalar mean_j, var_j

    n = rows(X)
    p = cols(X)

    X_std = J(n, p, .)
    means = J(1, p, .)
    sds = J(1, p, .)

    w_sum = sum(w)

    for (j = 1; j <= p; j++) {
        wX_j = w :* X[., j]
        mean_j = sum(wX_j) / w_sum
        means[j] = mean_j

        var_j = sum(w :* (X[., j] :- mean_j):^2) / w_sum
        sds[j] = sqrt(var_j)

        // Prevent division by zero for constant columns
        if (sds[j] < 1e-10) {
            sds[j] = 1
            X_std[., j] = J(n, 1, 0)
        } else {
            X_std[., j] = (X[., j] :- mean_j) / sds[j]
        }
    }
}

// ============================================================================
// UNSTANDARDIZE COEFFICIENTS
// Transforms coefficients back to original scale after estimation
// ============================================================================

real colvector unstandardize_coefs(real colvector beta_std,
                                   real rowvector means, real rowvector sds)
{
    real colvector beta

    beta = beta_std :/ sds'

    return(beta)
}

// ============================================================================
// CLAMP VALUES TO PREVENT NUMERICAL OVERFLOW
// ============================================================================

real colvector clamp(real colvector x, real scalar lo, real scalar hi)
{
    return(max((min((x, J(rows(x), 1, hi))), J(rows(x), 1, lo))))
}

// More efficient element-wise clamp
real colvector clamp_vec(real colvector x, real scalar lo, real scalar hi)
{
    real colvector result
    real scalar i, n

    n = rows(x)
    result = x

    for (i = 1; i <= n; i++) {
        if (result[i] < lo) result[i] = lo
        else if (result[i] > hi) result[i] = hi
    }

    return(result)
}

// ============================================================================
// IDENTIFY SELECTED VARIABLES (non-zero coefficients)
// ============================================================================

real colvector get_selected(real colvector beta, | real scalar tol)
{
    if (args() < 2 | tol == .) tol = 1e-10

    return(abs(beta) :> tol)
}

real scalar count_selected(real colvector beta, | real scalar tol)
{
    if (args() < 2 | tol == .) tol = 1e-10

    return(sum(abs(beta) :> tol))
}

// ============================================================================
// INFORMATION CRITERIA FOR MODEL SELECTION
// ============================================================================

real scalar compute_aic(real scalar ll, real scalar k)
{
    return(-2 * ll + 2 * k)
}

real scalar compute_bic(real scalar ll, real scalar k, real scalar n)
{
    return(-2 * ll + k * ln(n))
}

real scalar compute_ebic(real scalar ll, real scalar k, real scalar n,
                         real scalar p, | real scalar gamma)
{
    // Extended BIC with gamma parameter (default gamma = 0.5)
    if (args() < 5 | gamma == .) gamma = 0.5

    return(-2 * ll + k * ln(n) + 2 * gamma * k * ln(p))
}

// ============================================================================
// END OF UTILITY FUNCTIONS
// ============================================================================



// ============================================================================
// COMPRESS FE IDENTIFIERS TO CONTIGUOUS 1..G
// Maps arbitrary FE values to sequential integers for safe array indexing
// Returns: compressed IDs and number of groups per FE dimension
// ============================================================================

void compress_fe_ids(real matrix fe_ids_in, real matrix fe_ids_out,
                     real rowvector n_groups)
{
    real scalar n, n_fe, k, i
    real colvector unique_vals, fe_col
    real scalar n_unique

    n = rows(fe_ids_in)
    n_fe = cols(fe_ids_in)

    fe_ids_out = J(n, n_fe, .)
    n_groups = J(1, n_fe, .)

    for (k = 1; k <= n_fe; k++) {
        fe_col = fe_ids_in[., k]
        unique_vals = uniqrows(fe_col)
        n_unique = rows(unique_vals)
        n_groups[k] = n_unique

        // Map each value to its position in unique_vals (1..n_unique)
        for (i = 1; i <= n; i++) {
            fe_ids_out[i, k] = sum(unique_vals :<= fe_col[i])
        }
    }
}

// Single column version for convenience
real colvector compress_fe_col(real colvector fe_col, | real scalar n_groups)
{
    real colvector unique_vals, fe_compressed
    real scalar n, i

    n = rows(fe_col)
    unique_vals = uniqrows(fe_col)
    n_groups = rows(unique_vals)

    fe_compressed = J(n, 1, .)
    for (i = 1; i <= n; i++) {
        fe_compressed[i] = sum(unique_vals :<= fe_col[i])
    }

    return(fe_compressed)
}

// ============================================================================
// CLUSTER-ROBUST VARIANCE MATRIX FOR SCORES
// Computes: sum_g (sum_{i in g} x_i * e_i)' (sum_{i in g} x_i * e_i)
// Where g indexes clusters
// ============================================================================

real matrix cluster_matrix(real colvector scores, real colvector cluster_ids,
                           real matrix X)
{
    real scalar n, p, n_clusters, g, i
    real colvector unique_clusters, cluster_compressed
    real matrix cluster_scores, V

    n = rows(X)
    p = cols(X)

    // Compress cluster IDs to 1..G
    cluster_compressed = compress_fe_col(cluster_ids, n_clusters)

    // Accumulate scores by cluster: X_g' * e_g for each cluster g
    cluster_scores = J(n_clusters, p, 0)

    for (i = 1; i <= n; i++) {
        g = cluster_compressed[i]
        cluster_scores[g, .] = cluster_scores[g, .] + scores[i] * X[i, .]
    }

    // Compute variance: sum of outer products
    V = quadcross(cluster_scores, cluster_scores)

    return(V)
}

// ============================================================================
// COMPUTE WEIGHTED GROUP MEANS EFFICIENTLY
// For FE demeaning with compressed FE IDs
// ============================================================================

real colvector compute_group_means(real colvector v, real colvector w,
                                   real colvector fe_ids, real scalar n_groups)
{
    real scalar n, i, g
    real colvector group_sums, group_weights, group_means

    n = rows(v)
    group_sums = J(n_groups, 1, 0)
    group_weights = J(n_groups, 1, 0)

    // Accumulate
    for (i = 1; i <= n; i++) {
        g = fe_ids[i]
        group_sums[g] = group_sums[g] + w[i] * v[i]
        group_weights[g] = group_weights[g] + w[i]
    }

    // Compute means (avoid division by zero)
    group_means = group_sums :/ (group_weights :+ (group_weights :== 0))

    return(group_means)
}

// ============================================================================
// RECOVER FIXED EFFECT VALUES FROM ESTIMATION
// Implements PPML FOC iteration: FE_g = sum_g(y) / sum_g(mu)
// Similar to R penppml's compute_fes()
// ============================================================================

real colvector recover_fe_values(real colvector y, real colvector mu,
                                 real colvector fe_ids, real scalar n_groups,
                                 | real scalar tol, real scalar maxiter)
{
    real scalar n, iter, converged
    real colvector fe_values, y_sums, mu_adj_sums, adj
    real scalar crit, old_dev, dev, g, i
    real colvector dev_term

    if (args() < 5 | tol == .) tol = 1e-8
    if (args() < 6 | maxiter == .) maxiter = 100

    n = rows(y)

    // Initialize FE values to 1 (multiplicative)
    fe_values = J(n_groups, 1, 1)

    // Compute sum of y by group (constant)
    y_sums = J(n_groups, 1, 0)
    for (i = 1; i <= n; i++) {
        y_sums[fe_ids[i]] = y_sums[fe_ids[i]] + y[i]
    }

    old_dev = 0
    converged = 0

    for (iter = 1; iter <= maxiter; iter++) {
        // Compute sum of adjusted mu by group
        mu_adj_sums = J(n_groups, 1, 0)
        for (i = 1; i <= n; i++) {
            g = fe_ids[i]
            mu_adj_sums[g] = mu_adj_sums[g] + mu[i] * fe_values[g]
        }

        // Update FE values using PPML FOC: FE_g *= y_sum_g / mu_adj_sum_g
        adj = y_sums :/ (mu_adj_sums :+ (mu_adj_sums :== 0))
        adj = adj :* (y_sums :> 0) + (y_sums :== 0)  // Set to 1 where y_sum = 0
        fe_values = fe_values :* adj

        // Check convergence via deviance
        dev_term = J(n, 1, 0)
        for (i = 1; i <= n; i++) {
            g = fe_ids[i]
            if (y[i] > 0) {
                dev_term[i] = y[i] * ln(y[i] / (mu[i] * fe_values[g])) -
                              (y[i] - mu[i] * fe_values[g])
            } else {
                dev_term[i] = mu[i] * fe_values[g]
            }
        }
        dev = 2 * sum(dev_term)

        crit = abs(dev - old_dev) / max((min((dev, old_dev)), 0.1))
        if (crit < tol) {
            converged = 1
            break
        }
        old_dev = dev
    }

    return(fe_values)
}


// ============================================================================
// APPLY FE VALUES TO TEST DATA FOR CV SCORING
// Recovers FE contribution for out-of-sample observations
// Uses the approach from R penppml: estimate FE values from training data
// and apply them to test observations that share FE groups
// ============================================================================

// Binary search for a value in a sorted real vector.
// Returns index in [1..rows(v)] if found, else 0.
real scalar bsearch_sorted_real(real colvector v, real scalar key)
{
    real scalar lo, hi, mid

    if (key == .) return(0)
    lo = 1
    hi = rows(v)

    while (lo <= hi) {
        mid = floor((lo + hi) / 2)
        if (v[mid] == key) return(mid)
        if (v[mid] < key) lo = mid + 1
        else hi = mid - 1
    }

    return(0)
}

real colvector apply_fe_to_test(real colvector y_train, real colvector mu_train,
                                real colvector beta, real matrix X_test,
                                real matrix fe_train, real matrix fe_test)
{
    real scalar n_test, n_train, n_fe, k, i, g
    real colvector fe_contrib_test
    real colvector fe_train_compressed, n_groups_vec
    real scalar n_groups
    real colvector fe_values_k, y_sums, mu_sums, fe_test_k
    real matrix fe_test_compressed
    real scalar fe_val
    real colvector fe_train_k, fe_train_k_sorted
    real colvector fe_group_sorted
    real colvector fe_orig_by_group
    real colvector ord
    real scalar grp

    n_test = rows(X_test)
    n_train = rows(y_train)
    n_fe = cols(fe_train)

    // Initialize test FE contribution to 0 (log scale)
    fe_contrib_test = J(n_test, 1, 0)

    // Process each FE dimension
    for (k = 1; k <= n_fe; k++) {
        fe_train_k = fe_train[., k]

        // Compress training FE IDs
        fe_train_compressed = compress_fe_col(fe_train_k, n_groups)

        // Compute FE values from training data: FE_g = sum_g(y) / sum_g(mu)
        y_sums = J(n_groups, 1, 0)
        mu_sums = J(n_groups, 1, 0)

        for (i = 1; i <= n_train; i++) {
            g = fe_train_compressed[i]
            y_sums[g] = y_sums[g] + y_train[i]
            mu_sums[g] = mu_sums[g] + mu_train[i]
        }

        // FE values (multiplicative scale) avoiding division by zero
        fe_values_k = y_sums :/ (mu_sums :+ (mu_sums :== 0))
        fe_values_k = fe_values_k :* (mu_sums :> 0) + (mu_sums :== 0)

        // Build a fast lookup for test IDs by:
        // 1) sorting training IDs,
        // 2) recording one representative original ID per compressed group,
        // 3) binary searching that sorted representative vector.
        ord = order(fe_train_k, 1)
        fe_train_k_sorted = fe_train_k[ord]
        fe_group_sorted = fe_train_compressed[ord]

        fe_orig_by_group = J(n_groups, 1, .)
        for (i = 1; i <= n_train; i++) {
            g = fe_group_sorted[i]
            if (fe_orig_by_group[g] == .) fe_orig_by_group[g] = fe_train_k_sorted[i]
        }

        // Map test FE IDs to training FE groups
        fe_test_k = fe_test[., k]

        for (i = 1; i <= n_test; i++) {
            // Find matching group for this test observation
            fe_val = 1  // Default: no FE adjustment
            grp = bsearch_sorted_real(fe_orig_by_group, fe_test_k[i])
            if (grp > 0) fe_val = fe_values_k[grp]
            fe_contrib_test[i] = fe_contrib_test[i] + ln(max((fe_val, 1e-10)))
        }
    }

    return(fe_contrib_test)
}

end
