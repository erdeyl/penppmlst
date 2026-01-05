*! version 0.5.0  05jan2026
*! R-compatible functions for penppmlst package
*! These functions ensure mathematical equivalence with R penppml
*! Stata implementation by Erdey, László (2026)
*!   Faculty of Economics and Business, University of Debrecen, Hungary
*! Based on R penppml by Breinlich, Corradi, Rocha, Ruta, Santos Silva, Zylkin

version 17.0
mata:
mata set matastrict on

// ============================================================================
// POISSON DEVIANCE COMPUTATION (R-COMPATIBLE)
// Following R penppml: temp = -(y * log(y/mu) - (y-mu))
//                      temp[y==0] = -mu[y==0]
//                      deviance = -2 * sum(temp) / n
// ============================================================================

real scalar compute_deviance_r(real colvector y, real colvector mu)
{
    real scalar n
    real colvector temp
    real scalar i

    n = rows(y)
    temp = J(n, 1, 0)

    // R formula: temp = -(y * log(y/mu) - (y-mu))
    // For y > 0: temp = -(y * log(y/mu) - (y - mu))
    // For y = 0: temp = -mu
    for (i = 1; i <= n; i++) {
        if (y[i] == 0) {
            temp[i] = -mu[i]
        } else {
            temp[i] = -(y[i] * ln(y[i] / mu[i]) - (y[i] - mu[i]))
        }
    }

    return(-2 * sum(temp) / n)
}

// ============================================================================
// MU CLAMPING (R-COMPATIBLE)
// R uses 1e-190 to 1e190 bounds for mu
// ============================================================================

real colvector clamp_mu_r(real colvector mu)
{
    real colvector result
    real scalar i, n

    n = rows(mu)
    result = mu

    for (i = 1; i <= n; i++) {
        if (result[i] < 1e-190) result[i] = 1e-190
        else if (result[i] > 1e190) result[i] = 1e190
    }

    return(result)
}

// Stata-compatible mu clamp (narrower bounds for numerical stability)
real colvector clamp_mu_stata(real colvector mu)
{
    real colvector result
    real scalar i, n

    n = rows(mu)
    result = mu

    for (i = 1; i <= n; i++) {
        if (result[i] < 1e-10) result[i] = 1e-10
        else if (result[i] > 1e10) result[i] = 1e10
    }

    return(result)
}

// ============================================================================
// BIC COMPUTATION (R-COMPATIBLE)
// R penppml: bic = deviance + k * log(n) / n
// ============================================================================

real scalar compute_bic_r(real scalar deviance, real scalar k, real scalar n)
{
    return(deviance + k * ln(n) / n)
}

// ============================================================================
// PLUGIN LAMBDA COMPUTATION (Belloni, Chernozhukov, Hansen, Kozbur 2016)
// lambda = c * sqrt(n) * Phi^{-1}(1 - gamma / (2k))
// where gamma = 0.1 / log(n) by default
// ============================================================================

real scalar compute_plugin_lambda(real scalar n, real scalar k,
                                   | real scalar c_val, real scalar gamma_val)
{
    real scalar gamma, quantile, lambda

    // Set defaults following R penppml
    if (args() < 3 | c_val == .) c_val = 1.1
    if (args() < 4 | gamma_val == .) gamma_val = 0.1 / ln(n)

    gamma = gamma_val

    // Compute quantile of standard normal: Phi^{-1}(1 - gamma / (2k))
    quantile = invnormal(1 - gamma / (2 * k))

    // Lambda formula
    lambda = c_val * sqrt(n) * quantile

    return(lambda)
}

// ============================================================================
// CLUSTER MATRIX COMPUTATION FOR PLUGIN PENALTY WEIGHTS
// Following R penppml's cluster_matrix function:
// Computes XeeX matrix for cluster-robust SEs
// ============================================================================

real matrix compute_cluster_matrix(real colvector e, real colvector cluster_id,
                                    real matrix X)
{
    real scalar n, k, n_clusters, g
    real matrix XeeX
    real colvector obs_in_cluster
    real matrix X_cluster
    real colvector e_cluster
    real colvector score_cluster

    n = rows(e)
    k = cols(X)
    n_clusters = max(cluster_id)

    XeeX = J(k, k, 0)

    // Sum scores within clusters, then compute outer product
    for (g = 1; g <= n_clusters; g++) {
        obs_in_cluster = selectindex(cluster_id :== g)

        if (rows(obs_in_cluster) > 0) {
            X_cluster = X[obs_in_cluster, .]
            e_cluster = e[obs_in_cluster]

            // Sum of X_i * e_i within cluster
            score_cluster = quadcross(X_cluster', e_cluster)

            // Outer product
            XeeX = XeeX + score_cluster * score_cluster'
        }
    }

    return(XeeX)
}

// ============================================================================
// PLUGIN PHI COMPUTATION (coefficient-specific penalty weights)
// Following R penppml: phi = sqrt(diag(cluster_matrix(mu * z_resid, cluster, x_resid)) / n)
// ============================================================================

real colvector compute_plugin_phi(real colvector mu, real colvector z_resid,
                                   real colvector cluster_id, real matrix x_resid)
{
    real scalar n, k
    real matrix XeeX
    real colvector e
    real colvector phi

    n = rows(mu)
    k = cols(x_resid)

    // e = mu * z_resid (or mu * residuals)
    e = mu :* z_resid

    // Compute cluster-robust score variance
    XeeX = compute_cluster_matrix(e, cluster_id, x_resid)

    // phi = sqrt(diag(XeeX) / n)
    phi = sqrt(diagonal(XeeX) / n)

    return(phi)
}

// ============================================================================
// GLMNET-STYLE LAMBDA SCALING FOR PLUGIN LASSO
// Following R penppml: lambda_glmnet = lambda / sum(y) * sum(phi) / k
// ============================================================================

real scalar scale_lambda_glmnet(real scalar lambda, real colvector y,
                                 real colvector phi)
{
    real scalar k, sum_y, sum_phi

    k = rows(phi)
    sum_y = sum(y)
    sum_phi = sum(phi)

    return(lambda / sum_y * sum_phi / k)
}

// ============================================================================
// CONVERGENCE CRITERION (R-compatible)
// Following R penppml: crit = abs(delta_deviance) / max(min(deviance, old_deviance), 0.1)
// ============================================================================

real scalar compute_convergence_criterion(real scalar deviance,
                                           real scalar old_deviance)
{
    real scalar delta_deviance, denom_crit

    delta_deviance = old_deviance - deviance

    // Handle case where deviance decreased significantly
    if (deviance < 0.1 * delta_deviance) {
        delta_deviance = deviance
    }

    denom_crit = max((min((deviance, old_deviance)), 0.1))

    return(abs(delta_deviance) / denom_crit)
}

// ============================================================================
// MU INITIALIZATION (R-compatible)
// Following R penppml: mu = (y + mean(y)) / 2
// ============================================================================

real colvector init_mu_r(real colvector y)
{
    return((y :+ mean(y)) / 2)
}

// ============================================================================
// WORKING RESPONSE TRANSFORMATION (R-compatible)
// z = (y - mu) / mu + log(mu)
// ============================================================================

real colvector compute_working_response(real colvector y, real colvector mu)
{
    return((y - mu) :/ mu + ln(mu))
}

// ============================================================================
// END OF R-COMPATIBLE FUNCTIONS
// ============================================================================

end
