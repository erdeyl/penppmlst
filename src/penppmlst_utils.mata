*! version 0.1.0  03jan2026
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

end
