*! version 0.5.0  05jan2026
*! Core PenPPML class: Penalized PPML with High-Dimensional Fixed Effects
*! Stata implementation by Erdey, László (2026)
*!   Faculty of Economics and Business, University of Debrecen, Hungary
*! Based on R penppml by Breinlich, Corradi, Rocha, Ruta, Santos Silva, Zylkin

version 17.0
mata:
mata set matastrict on

// ============================================================================
// PenPPML CLASS DEFINITION
// Implements Penalized Poisson Pseudo Maximum Likelihood with HDFE
// ============================================================================

class PenPPML {
    // ===== DATA =====
    real colvector      y               // Dependent variable (counts)
    real matrix         X               // Regressors to be penalized
    real matrix         Z               // Unpenalized controls (optional)
    real colvector      w               // User-specified weights
    real colvector      offset          // Offset term (optional)
    real scalar         n               // Number of observations
    real scalar         p               // Number of penalized regressors
    real scalar         q               // Number of unpenalized controls

    // ===== FIXED EFFECTS =====
    string rowvector    fe_vars         // FE variable names
    real scalar         n_fe            // Number of FE sets
    real matrix         fe_ids          // FE group identifiers
    real colvector      fe_contrib      // FE contribution to linear predictor

    // ===== ALGORITHM PARAMETERS =====
    real scalar         tol             // Outer IRLS convergence tolerance
    real scalar         tol_inner       // Inner coordinate descent tolerance
    real scalar         maxiter         // Max IRLS iterations
    real scalar         maxiter_inner   // Max CD iterations per IRLS step
    real scalar         standardize     // Standardize X before penalization
    real scalar         verbose         // Print iteration info

    // ===== PENALTY PARAMETERS =====
    string scalar       penalty         // "lasso", "ridge", or "elasticnet"
    real scalar         lambda          // Overall penalty level
    real scalar         alpha           // Elastic net mixing (1=lasso, 0=ridge)
    real colvector      psi             // Coefficient-specific penalty loadings

    // ===== RESULTS =====
    real colvector      beta            // Penalized coefficient estimates
    real colvector      beta_unpen      // Unpenalized control coefficients
    real colvector      beta_post       // Post-lasso coefficients (if requested)
    real colvector      mu              // Fitted values exp(eta)
    real colvector      eta             // Linear predictor
    real matrix         V               // Variance-covariance matrix (post-lasso)
    real scalar         deviance        // Final deviance
    real scalar         ll              // Log pseudo-likelihood
    real scalar         converged       // Convergence flag (1=yes)
    real scalar         iterations      // Iterations used
    real colvector      selected        // Indicator for selected variables (1=selected)
    real scalar         n_selected      // Number of selected variables

    // ===== STANDARDIZATION =====
    real rowvector      X_means         // Column means (for unstandardizing)
    real rowvector      X_sds           // Column SDs (for unstandardizing)
    real matrix         X_std           // Standardized X

    // ===== HDFE SETTINGS =====
    string scalar       hdfe_method     // "mata" or "ppmlhdfe"
    real scalar         r_compatible    // Use R-compatible settings

    // ===== METHODS =====
    void                new()
    void                init()
    void                set_data()
    void                set_fe()
    void                set_penalty()
    void                set_options()
    real scalar         solve()
    void                irls_step()
    real colvector      partial_out_fe()
    void                compute_post_lasso()
    void                compute_vcov()
    real scalar         compute_ll()
    void                print_header()
    void                print_iteration()
    void                print_results()
    void                set_hdfe_method()
    real colvector      partial_out_fe_mata()
    real colvector      partial_out_fe_ppmlhdfe()

    // ===== PPMLHDFE INTEGRATION =====
    pointer(class FixedEffects scalar) scalar HDFE_ptr  // Pointer to reghdfe's FixedEffects
    real scalar         hdfe_initialized                 // Has HDFE been initialized?
    void                init_hdfe_object()               // Initialize FixedEffects object
    void                update_hdfe_weights()            // Update HDFE weights
}

// ============================================================================
// CONSTRUCTOR
// ============================================================================

void PenPPML::new()
{
    // Set defaults
    tol = 1e-8
    tol_inner = 1e-7
    maxiter = 1000
    maxiter_inner = 1000
    standardize = 1
    verbose = 0

    penalty = "lasso"
    lambda = 0
    alpha = 1

    converged = 0
    iterations = 0
    n_fe = 0
    q = 0

    hdfe_method = "mata"
    r_compatible = 0
    hdfe_initialized = 0
    HDFE_ptr = NULL
}

// ============================================================================
// INITIALIZATION
// ============================================================================

void PenPPML::init()
{
    real scalar i

    n = rows(y)
    p = cols(X)

    // Initialize weights if not provided
    if (rows(w) == 0) {
        w = J(n, 1, 1)
    }

    // Initialize offset if not provided
    if (rows(offset) == 0) {
        offset = J(n, 1, 0)
    }

    // Initialize penalty loadings if not provided (uniform weights)
    if (rows(psi) == 0) {
        psi = J(p, 1, 1)
    }

    // Standardize X if requested
    if (standardize) {
        standardize_weighted(X, w, X_std, X_means, X_sds)
    } else {
        X_std = X
        X_means = J(1, p, 0)
        X_sds = J(1, p, 1)
    }

    // Initialize mu with simple starting values
    // Use (y + mean(y)) / 2 to avoid mu = 0
    mu = (y :+ mean(y)) / 2
    mu = clamp_vec(mu, 1e-10, 1e10)

    // Initialize linear predictor
    eta = ln(mu)

    // Initialize coefficients
    beta = J(p, 1, 0)
    if (q > 0) {
        beta_unpen = J(q, 1, 0)
    }

    // Initialize selected indicators
    selected = J(p, 1, 0)
    n_selected = 0

    // Initialize FE contribution
    if (n_fe > 0) {
        fe_contrib = J(n, 1, 0)
    } else {
        fe_contrib = J(n, 1, 0)
    }
}

// ============================================================================
// SET DATA
// ============================================================================

void PenPPML::set_data(real colvector y_in, real matrix X_in,
                       | real matrix Z_in, real colvector w_in,
                         real colvector offset_in)
{
    y = y_in
    X = X_in

    if (args() >= 3 & rows(Z_in) > 0) {
        Z = Z_in
        q = cols(Z)
    }

    if (args() >= 4 & rows(w_in) > 0) {
        w = w_in
    }

    if (args() >= 5 & rows(offset_in) > 0) {
        offset = offset_in
    }
}

// ============================================================================
// SET FIXED EFFECTS
// Simple implementation using within-transformation
// For production, integrate with reghdfe's FixedEffects class
// ============================================================================

void PenPPML::set_fe(real matrix fe_ids_in, | string rowvector fe_names_in)
{
    fe_ids = fe_ids_in
    n_fe = cols(fe_ids)

    if (args() >= 2) {
        fe_vars = fe_names_in
    }
}

// ============================================================================
// SET PENALTY PARAMETERS
// ============================================================================

void PenPPML::set_penalty(string scalar penalty_in, real scalar lambda_in,
                          | real scalar alpha_in, real colvector psi_in)
{
    penalty = penalty_in
    lambda = lambda_in

    if (args() >= 3) {
        alpha = alpha_in
    }

    if (args() >= 4 & rows(psi_in) > 0) {
        psi = psi_in
    }
}

// ============================================================================
// SET ALGORITHM OPTIONS
// ============================================================================

void PenPPML::set_options(| real scalar tol_in, real scalar maxiter_in,
                            real scalar standardize_in, real scalar verbose_in)
{
    if (args() >= 1 & tol_in != .) tol = tol_in
    if (args() >= 2 & maxiter_in != .) maxiter = maxiter_in
    if (args() >= 3 & standardize_in != .) standardize = standardize_in
    if (args() >= 4 & verbose_in != .) verbose = verbose_in
}

// ============================================================================
// PARTIAL OUT FIXED EFFECTS (DISPATCHER)
// Routes to appropriate backend based on hdfe_method setting
// ============================================================================

real colvector PenPPML::partial_out_fe(real colvector v, real colvector wts)
{
    // Dispatch based on hdfe_method
    if (hdfe_method == "ppmlhdfe") {
        return(partial_out_fe_ppmlhdfe(v, wts))
    } else {
        // Default: pure Mata implementation (R-compatible)
        return(partial_out_fe_mata(v, wts))
    }
}

// ============================================================================
// PARTIAL OUT FE VIA MATA (R-COMPATIBLE)
// Implements alternating projections (Gaure 2013) for multiple FE sets
// Each variable is demeaned within groups, iterating until convergence
// This is the default method, compatible with R penppml's fixest approach
// ============================================================================

real colvector PenPPML::partial_out_fe_mata(real colvector v, real colvector wts)
{
    real colvector v_demean
    real scalar fe_iter, max_fe_iter, fe_converged
    real colvector v_old
    real scalar k, g, n_groups
    real colvector group_sums, group_weights, group_means
    real scalar max_change

    if (n_fe == 0) {
        return(v)
    }

    max_fe_iter = 100
    v_demean = v

    fe_converged = 0
    for (fe_iter = 1; fe_iter <= max_fe_iter; fe_iter++) {
        v_old = v_demean

        // Demean within each FE set
        for (k = 1; k <= n_fe; k++) {
            n_groups = max(fe_ids[., k])
            group_sums = J(n_groups, 1, 0)
            group_weights = J(n_groups, 1, 0)

            // Accumulate weighted sums and weights by group
            for (g = 1; g <= n; g++) {
                group_sums[fe_ids[g, k]] = group_sums[fe_ids[g, k]] + wts[g] * v_demean[g]
                group_weights[fe_ids[g, k]] = group_weights[fe_ids[g, k]] + wts[g]
            }

            // Compute group means
            group_means = group_sums :/ (group_weights :+ (group_weights :== 0))

            // Subtract group means
            for (g = 1; g <= n; g++) {
                v_demean[g] = v_demean[g] - group_means[fe_ids[g, k]]
            }
        }

        // Check convergence
        max_change = max(abs(v_demean - v_old))
        if (max_change < 1e-10) {
            fe_converged = 1
            break
        }
    }

    return(v_demean)
}

// ============================================================================
// MAIN SOLVER: IRLS WITH PENALTY
// ============================================================================

real scalar PenPPML::solve()
{
    real scalar iter
    real colvector z, irls_w, sqrt_w
    real matrix X_tilde, Z_tilde
    real colvector z_tilde
    real colvector beta_new
    real scalar deviance_old, deviance_new
    real scalar eps, max_coef_change
    real scalar j

    // Initialize
    init()

    if (verbose) {
        print_header()
    }

    // Initial deviance
    deviance = compute_deviance(y, mu)

    // ===== OUTER IRLS LOOP =====
    for (iter = 1; iter <= maxiter; iter++) {

        // ----- Step 1: Working variable transformation -----
        // z = (y - mu) / mu + log(mu)  [working response]
        // w_irls = mu                   [Poisson variance = mu]
        z = (y - mu) :/ mu + ln(mu)
        irls_w = w :* mu

        // Numerical safeguards
        irls_w = clamp_vec(irls_w, 1e-10, 1e10)
        sqrt_w = sqrt(irls_w)

        // ----- Step 2: Partial out fixed effects -----
        if (n_fe > 0) {
            z_tilde = partial_out_fe(z, irls_w)

            // Partial out X columns
            X_tilde = J(n, p, .)
            for (j = 1; j <= p; j++) {
                X_tilde[., j] = partial_out_fe(X_std[., j], irls_w)
            }

            // Partial out Z columns if present
            if (q > 0) {
                Z_tilde = J(n, q, .)
                for (j = 1; j <= q; j++) {
                    Z_tilde[., j] = partial_out_fe(Z[., j], irls_w)
                }
            }
        } else {
            z_tilde = z
            X_tilde = X_std
            if (q > 0) Z_tilde = Z
        }

        // ----- Step 3: Penalized weighted least squares -----
        if (penalty == "ridge" | alpha == 0) {
            // Ridge regression: analytical solution
            beta_new = ridge_solve(X_tilde, z_tilde, irls_w, lambda)
        } else if (penalty == "lasso" | penalty == "elasticnet") {
            // Lasso or elastic net: coordinate descent
            beta_new = coordinate_descent(X_tilde, z_tilde, irls_w,
                                          lambda, psi, alpha,
                                          tol_inner, maxiter_inner, 0)
        } else {
            // Unknown penalty - use lasso as default
            beta_new = coordinate_descent(X_tilde, z_tilde, irls_w,
                                          lambda, psi, 1,
                                          tol_inner, maxiter_inner, 0)
        }

        beta = beta_new

        // ----- Step 4: Update linear predictor and mu -----
        // Unstandardize coefficients for prediction
        eta = X * unstandardize_coefs(beta, X_means, X_sds) + offset

        // Add unpenalized controls
        if (q > 0) {
            // Fit unpenalized controls via WLS on residual
            // (simplified: could integrate into coordinate descent)
            beta_unpen = ridge_solve(Z_tilde, z_tilde - X_tilde * beta, irls_w, 0)
            eta = eta + Z * beta_unpen
        }

        // Add FE contribution (recompute from working variable residuals)
        if (n_fe > 0) {
            fe_contrib = z - z_tilde - (X_std * beta)
            if (q > 0) fe_contrib = fe_contrib - (Z * beta_unpen)
            eta = eta + fe_contrib
        }

        // Update mu = exp(eta) with numerical bounds
        mu = exp(eta)
        mu = clamp_vec(mu, 1e-10, 1e10)

        // ----- Step 5: Check convergence -----
        deviance_old = deviance
        deviance_new = compute_deviance(y, mu)

        // Relative deviance change
        eps = abs(deviance_new - deviance_old) / max((abs(deviance_old), 0.1))

        if (verbose) {
            print_iteration(iter, deviance_new, eps, count_selected(beta))
        }

        deviance = deviance_new

        if (eps < tol) {
            converged = 1
            iterations = iter
            break
        }

        iterations = iter
    }

    // ===== POST-PROCESSING =====

    // Unstandardize final coefficients
    beta = unstandardize_coefs(beta, X_means, X_sds)

    // Identify selected variables
    selected = get_selected(beta)
    n_selected = sum(selected)

    // Compute log-likelihood
    ll = compute_ll()

    if (verbose) {
        print_results()
    }

    return(converged)
}

// ============================================================================
// COMPUTE POST-LASSO ESTIMATES
// Refit unpenalized PPML on selected variables only
// ============================================================================

void PenPPML::compute_post_lasso()
{
    real matrix X_selected
    real colvector beta_sel
    real scalar iter, post_converged
    real colvector z, irls_w, z_tilde
    real matrix X_sel_tilde
    real scalar j, k, idx
    real colvector sel_idx
    real scalar n_sel
    real scalar deviance_old, deviance_new, eps

    n_sel = n_selected
    if (n_sel == 0) {
        beta_post = J(p, 1, 0)
        return
    }

    // Get indices of selected variables
    sel_idx = select(1::p, selected)

    // Extract selected columns
    X_selected = X[., sel_idx]

    // Initialize from current estimates
    beta_sel = beta[sel_idx]
    mu = exp(X_selected * beta_sel + offset)
    if (q > 0) mu = mu :* exp(Z * beta_unpen)
    if (n_fe > 0) mu = mu :* exp(fe_contrib)
    mu = clamp_vec(mu, 1e-10, 1e10)

    deviance = compute_deviance(y, mu)

    // IRLS for unpenalized estimation on selected variables
    post_converged = 0
    for (iter = 1; iter <= maxiter; iter++) {
        z = (y - mu) :/ mu + ln(mu)
        irls_w = w :* mu
        irls_w = clamp_vec(irls_w, 1e-10, 1e10)

        if (n_fe > 0) {
            z_tilde = partial_out_fe(z, irls_w)
            X_sel_tilde = J(n, n_sel, .)
            for (k = 1; k <= n_sel; k++) {
                X_sel_tilde[., k] = partial_out_fe(X_selected[., k], irls_w)
            }
        } else {
            z_tilde = z
            X_sel_tilde = X_selected
        }

        // Unpenalized WLS
        beta_sel = ridge_solve(X_sel_tilde, z_tilde, irls_w, 0)

        // Update eta and mu
        eta = X_selected * beta_sel + offset
        if (q > 0) eta = eta + Z * beta_unpen
        if (n_fe > 0) {
            fe_contrib = z - z_tilde - (X_sel_tilde * beta_sel)
            eta = eta + fe_contrib
        }

        mu = exp(eta)
        mu = clamp_vec(mu, 1e-10, 1e10)

        deviance_old = deviance
        deviance_new = compute_deviance(y, mu)
        eps = abs(deviance_new - deviance_old) / max((abs(deviance_old), 0.1))
        deviance = deviance_new

        if (eps < tol) {
            post_converged = 1
            break
        }
    }

    // Store post-lasso coefficients
    beta_post = J(p, 1, 0)
    for (k = 1; k <= n_sel; k++) {
        beta_post[sel_idx[k]] = beta_sel[k]
    }
}

// ============================================================================
// COMPUTE VARIANCE-COVARIANCE MATRIX
// Sandwich estimator for post-lasso coefficients
// ============================================================================

void PenPPML::compute_vcov()
{
    real matrix X_selected
    real colvector sel_idx
    real scalar n_sel, k
    real matrix H, G, V_sel
    real colvector score_i
    real scalar i

    n_sel = n_selected
    if (n_sel == 0) {
        V = J(p, p, 0)
        return
    }

    sel_idx = select(1::p, selected)
    X_selected = X[., sel_idx]

    // Hessian (expected information): H = X'WX where W = diag(mu)
    H = quadcross(X_selected, mu :* w, X_selected)

    // Score outer product: G = sum_i (y_i - mu_i)^2 * x_i * x_i'
    G = J(n_sel, n_sel, 0)
    for (i = 1; i <= n; i++) {
        score_i = (y[i] - mu[i]) * X_selected[i, .]'
        G = G + w[i] * score_i * score_i'
    }

    // Sandwich: V = H^{-1} G H^{-1}
    V_sel = invsym(H) * G * invsym(H)

    // Expand to full coefficient vector
    V = J(p, p, 0)
    for (i = 1; i <= n_sel; i++) {
        for (k = 1; k <= n_sel; k++) {
            V[sel_idx[i], sel_idx[k]] = V_sel[i, k]
        }
    }
}

// ============================================================================
// COMPUTE LOG PSEUDO-LIKELIHOOD
// ll = sum[ y*log(mu) - mu - log(y!) ]
// ============================================================================

real scalar PenPPML::compute_ll()
{
    real colvector ll_i

    ll_i = y :* ln(mu) - mu - lngamma(y :+ 1)
    ll_i = ll_i :* w

    return(sum(ll_i))
}

// ============================================================================
// PRINT FUNCTIONS
// ============================================================================

void PenPPML::print_header()
{
    printf("\n")
    printf("{txt}Penalized PPML Estimation\n")
    printf("{hline 60}\n")
    printf("Penalty: %s", penalty)
    if (penalty == "elasticnet") {
        printf(" (alpha = %g)", alpha)
    }
    printf("\n")
    printf("Lambda: %g\n", lambda)
    printf("N = %g, p = %g", n, p)
    if (n_fe > 0) {
        printf(", FE sets = %g", n_fe)
    }
    printf("\n")
    printf("{hline 60}\n")
    printf("{txt}%5s  %14s  %12s  %8s\n", "Iter", "Deviance", "Change", "Selected")
    printf("{hline 60}\n")
}

void PenPPML::print_iteration(real scalar iter, real scalar dev,
                              real scalar change, real scalar nsel)
{
    printf("{txt}%5.0f  %14.6f  %12.2e  %8.0f\n", iter, dev, change, nsel)
}

void PenPPML::print_results()
{
    printf("{hline 60}\n")
    if (converged) {
        printf("{txt}Converged in %g iterations\n", iterations)
    } else {
        printf("{err}Warning: Did not converge in %g iterations\n", iterations)
    }
    printf("{txt}Final deviance: %g\n", deviance)
    printf("Log pseudo-likelihood: %g\n", ll)
    printf("Variables selected: %g / %g\n", n_selected, p)
    printf("{hline 60}\n")
}

// ============================================================================
// SET HDFE METHOD
// ============================================================================

void PenPPML::set_hdfe_method(string scalar method, real scalar r_compat)
{
    hdfe_method = method
    r_compatible = r_compat
}

// ============================================================================
// PARTIAL OUT FE VIA PPMLHDFE/REGHDFE
// Uses reghdfe's FixedEffects Mata class (used by ppmlhdfe)
// This provides optimized HDFE computation via alternating projections
// with acceleration methods (Cimmino, symmetric Kaczmarz, etc.)
// ============================================================================

real colvector PenPPML::partial_out_fe_ppmlhdfe(real colvector v, real colvector wts)
{
    real matrix data
    real colvector v_demean

    if (n_fe == 0) {
        return(v)
    }

    // Initialize HDFE object if not done yet
    if (!hdfe_initialized) {
        init_hdfe_object()
    }

    // Update weights in the HDFE object
    update_hdfe_weights(wts)

    // Create data matrix for partial_out (column vector)
    data = v

    // Call reghdfe's _partial_out method
    // Parameters: data, standardize=0, drop_singletons=0, demean=0, flush_aux=1
    (*HDFE_ptr)._partial_out(data, 0, 0, 0, 1)

    v_demean = data

    return(v_demean)
}

// ============================================================================
// INITIALIZE HDFE OBJECT
// Creates reghdfe's FixedEffects object for use in partial-out operations
// ============================================================================

void PenPPML::init_hdfe_object()
{
    string scalar absorb_spec
    string scalar touse_var
    real scalar k
    
    // Build absorb specification from fe_vars
    absorb_spec = ""
    for (k = 1; k <= n_fe; k++) {
        if (k > 1) absorb_spec = absorb_spec + " "
        absorb_spec = absorb_spec + fe_vars[k]
    }
    
    // Get touse variable name (all obs in current sample)
    touse_var = st_tempname()
    (void) st_addvar("byte", touse_var)
    st_store(., touse_var, J(n, 1, 1))
    
    // Create FixedEffects object via fixed_effects() function
    // Parameters: absorb, touse, weight_type, depvar, ., verbose
    HDFE_ptr = &fixed_effects(absorb_spec, touse_var, "aweight", "", ., 0)
    
    // Load initial weights (will be updated during IRLS)
    (*HDFE_ptr).load_weights("aweight", "", w, 0)
    
    hdfe_initialized = 1
}

// ============================================================================
// UPDATE HDFE WEIGHTS
// Updates the weights in the FixedEffects object for weighted demeaning
// ============================================================================

void PenPPML::update_hdfe_weights(real colvector wts)
{
    if (!hdfe_initialized) {
        errprintf("HDFE object not initialized
")
        exit(error(198))
    }
    
    // Update sorted weights in the HDFE object
    (*HDFE_ptr).update_sorted_weights(wts)
    
    // Update cached objects that depend on weights
    (*HDFE_ptr).update_cvar_objects()
}

// ============================================================================
// END OF PenPPML CLASS
// ============================================================================

end
