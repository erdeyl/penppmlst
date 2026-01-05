*! version 0.5.0  05jan2026
*! Plugin lasso routines for penppmlst package
*! Implements penalty weight computation following Belloni, Chernozhukov, Hansen, Kozbur (2016)
*! Stata implementation by Erdey, László (2026)
*!   Faculty of Economics and Business, University of Debrecen, Hungary
*! Authors: Based on R penppml by Breinlich, Corradi, Rocha, Ruta, Santos Silva, Zylkin

version 17.0
mata:
mata set matastrict on

// ============================================================================
// PLUGIN LASSO CLASS
// Computes data-driven penalty weights and lambda for robust variable selection
// ============================================================================

class PenPPML_Plugin {
    // Data
    real colvector      y
    real matrix         X
    real matrix         fe_ids
    real colvector      w
    real colvector      cluster_id
    real scalar         n
    real scalar         p
    real scalar         n_clusters

    // Initial estimates (from unpenalized fit)
    real colvector      mu_init
    real colvector      residuals
    real matrix         X_tilde         // Partialled-out X

    // Plugin parameters
    real scalar         c               // Slack parameter (default: 1.1)
    real scalar         gamma           // Significance level parameter

    // Results
    real colvector      psi             // Penalty loadings (p x 1)
    real scalar         lambda          // Overall penalty level
    real colvector      score_var       // Score variance by coefficient

    // Methods
    void                new()
    void                set_data()
    void                set_cluster()
    void                set_parameters()
    void                compute_initial_fit()
    void                compute_penalty_loadings()
    void                compute_lambda()
    real colvector      get_psi()
    real scalar         get_lambda()
    void                print_summary()
}

// ============================================================================
// CONSTRUCTOR
// ============================================================================

void PenPPML_Plugin::new()
{
    c = 1.1             // Belloni et al. recommend c > 1 for finite samples
    gamma = .           // Will be computed based on n, p
    n_clusters = 0
}

// ============================================================================
// SET DATA
// ============================================================================

void PenPPML_Plugin::set_data(real colvector y_in, real matrix X_in,
                               | real matrix fe_ids_in, real colvector w_in)
{
    y = y_in
    X = X_in
    n = rows(y)
    p = cols(X)

    if (args() >= 3 & rows(fe_ids_in) > 0) {
        fe_ids = fe_ids_in
    }
    else {
        fe_ids = J(0, 0, .)
    }

    if (args() >= 4 & rows(w_in) > 0) {
        w = w_in
    }
    else {
        w = J(n, 1, 1)
    }
}

// ============================================================================
// SET CLUSTER VARIABLE
// ============================================================================

void PenPPML_Plugin::set_cluster(real colvector cluster_id_in)
{
    cluster_id = cluster_id_in
    n_clusters = max(cluster_id)
}

// ============================================================================
// SET PLUGIN PARAMETERS
// ============================================================================

void PenPPML_Plugin::set_parameters(| real scalar c_in, real scalar gamma_in)
{
    if (args() >= 1 & c_in != .) c = c_in
    if (args() >= 2 & gamma_in != .) gamma = gamma_in
}

// ============================================================================
// COMPUTE INITIAL FIT (UNPENALIZED PPML)
// ============================================================================

void PenPPML_Plugin::compute_initial_fit()
{
    class PenPPML scalar M
    real scalar j

    // Fit unpenalized PPML to get residuals
    M = PenPPML()
    M.set_data(y, X, J(n, 0, .), w)
    if (rows(fe_ids) > 0) {
        M.set_fe(fe_ids)
    }
    // Use ridge with lambda=0 (effectively OLS initialization)
    M.set_penalty("ridge", 1e-6, 0, J(p, 1, 1))
    M.set_options(1e-8, 200, 1, 0)
    M.solve()

    mu_init = M.mu
    residuals = y - mu_init

    // Get partialled-out X (without FE)
    // For simplicity, use centered X if no FE
    if (rows(fe_ids) > 0) {
        X_tilde = J(n, p, .)
        for (j = 1; j <= p; j++) {
            X_tilde[., j] = M.partial_out_fe(X[., j], w :* mu_init)
        }
    }
    else {
        X_tilde = X :- mean(X, w)
    }
}

// ============================================================================
// COMPUTE PENALTY LOADINGS (PSI)
// Following Belloni, Chernozhukov, Hansen, Kozbur (2016)
// ============================================================================

void PenPPML_Plugin::compute_penalty_loadings()
{
    real scalar j, g
    real colvector score_j
    real scalar var_j
    real colvector cluster_sum
    real matrix cluster_scores

    // Compute initial fit if not done
    if (rows(mu_init) == 0) {
        compute_initial_fit()
    }

    psi = J(p, 1, .)
    score_var = J(p, 1, .)

    // For each coefficient, compute score variance
    for (j = 1; j <= p; j++) {
        // Score for coefficient j: X_j * (y - mu) / sqrt(mu)
        // Under Poisson, Var(y) = mu, so score is X_j * residual
        score_j = X_tilde[., j] :* residuals :* sqrt(w)

        if (n_clusters > 0) {
            // Cluster-robust variance
            // Sum scores within clusters, then compute variance
            cluster_sum = J(n_clusters, 1, 0)

            for (g = 1; g <= n_clusters; g++) {
                cluster_sum[g] = sum(select(score_j, cluster_id :== g))
            }

            var_j = sum(cluster_sum :^ 2) / n
        }
        else {
            // Heteroskedasticity-robust variance (no clustering)
            var_j = sum(score_j :^ 2) / n
        }

        score_var[j] = var_j
        psi[j] = sqrt(var_j)
    }

    // Normalize penalty loadings
    // Option 1: Divide by max (keeps loadings between 0 and 1)
    psi = psi / max(psi)

    // Ensure minimum loading (prevents division issues)
    psi = rowmax((psi, J(p, 1, 0.01)))
}

// ============================================================================
// COMPUTE OVERALL LAMBDA
// Lambda formula from Belloni et al.:
// lambda = c * sqrt(n) * Phi^{-1}(1 - gamma / (2p))
// where gamma is the significance level (probability of false selection)
// ============================================================================

void PenPPML_Plugin::compute_lambda()
{
    real scalar quantile

    // Set gamma if not specified
    // Default: gamma = 0.1 / log(n) following Belloni et al.
    if (gamma == .) {
        gamma = 0.1 / ln(n)
    }

    // Compute quantile of standard normal
    // Phi^{-1}(1 - gamma / (2p))
    quantile = invnormal(1 - gamma / (2 * p))

    // Lambda formula
    lambda = c * sqrt(n) * quantile

    // Scale by average penalty loading for interpretability
    // This ensures lambda is on a comparable scale to CV-selected lambda
    lambda = lambda * mean(psi)
}

// ============================================================================
// GETTERS
// ============================================================================

real colvector PenPPML_Plugin::get_psi()
{
    if (rows(psi) == 0) {
        compute_penalty_loadings()
    }
    return(psi)
}

real scalar PenPPML_Plugin::get_lambda()
{
    if (lambda == .) {
        compute_lambda()
    }
    return(lambda)
}

// ============================================================================
// PRINT SUMMARY
// ============================================================================

void PenPPML_Plugin::print_summary()
{
    real scalar j
    real colvector sorted_psi
    real colvector sort_idx

    printf("\n{txt}Plugin Lasso Summary\n")
    printf("{hline 60}\n")
    printf("Observations: %g\n", n)
    printf("Coefficients: %g\n", p)
    if (n_clusters > 0) {
        printf("Clusters: %g\n", n_clusters)
    }
    printf("Slack parameter (c): %g\n", c)
    printf("Gamma: %g\n", gamma)
    printf("Lambda: %12.6f\n", lambda)
    printf("{hline 60}\n")
    printf("Penalty loadings summary:\n")
    printf("  Min:    %8.4f\n", min(psi))
    printf("  Median: %8.4f\n", mm_median(psi))
    printf("  Max:    %8.4f\n", max(psi))
    printf("  Mean:   %8.4f\n", mean(psi))
    printf("{hline 60}\n")
}

// ============================================================================
// ITERATED PLUGIN LASSO
// Refines penalty weights iteratively using post-lasso residuals
// ============================================================================

class PenPPML_IteratedPlugin {
    // Inherits from Plugin but iterates
    class PenPPML_Plugin scalar plugin

    // Iteration settings
    real scalar         maxiter
    real scalar         tol
    real scalar         verbose

    // Results
    real colvector      psi_final
    real scalar         lambda_final
    real scalar         iterations
    real scalar         converged

    // Methods
    void                new()
    void                run()
}

void PenPPML_IteratedPlugin::new()
{
    maxiter = 5
    tol = 0.1           // Relative change in psi
    verbose = 0
}

void PenPPML_IteratedPlugin::run(real colvector y, real matrix X,
                                  | real matrix fe_ids, real colvector w,
                                    real colvector cluster_id)
{
    real scalar iter
    real colvector psi_old, psi_new
    real scalar max_change
    class PenPPML scalar M
    real scalar n, p

    n = rows(y)
    p = cols(X)

    // Initialize plugin
    plugin = PenPPML_Plugin()
    plugin.set_data(y, X, fe_ids, w)
    if (rows(cluster_id) > 0) {
        plugin.set_cluster(cluster_id)
    }

    // First iteration: compute initial penalty loadings
    plugin.compute_penalty_loadings()
    plugin.compute_lambda()
    psi_old = plugin.psi

    converged = 0
    for (iter = 1; iter <= maxiter; iter++) {
        // Fit penalized model with current weights
        M = PenPPML()
        M.set_data(y, X, J(n, 0, .), w)
        if (rows(fe_ids) > 0) {
            M.set_fe(fe_ids)
        }
        M.set_penalty("lasso", plugin.lambda, 1, plugin.psi)
        M.set_options(1e-8, 100, 1, 0)
        M.solve()

        // Update residuals for next iteration
        plugin.mu_init = M.mu
        plugin.residuals = y - M.mu

        // Recompute penalty loadings
        plugin.compute_penalty_loadings()
        psi_new = plugin.psi

        // Check convergence
        max_change = max(abs(psi_new - psi_old) :/ (psi_old :+ 0.01))

        if (verbose) {
            printf("Iteration %g: max psi change = %g\n", iter, max_change)
        }

        if (max_change < tol) {
            converged = 1
            iterations = iter
            break
        }

        psi_old = psi_new
        iterations = iter
    }

    // Final results
    psi_final = plugin.psi
    lambda_final = plugin.lambda
}

// ============================================================================
// BOOTSTRAP LASSO
// Runs plugin lasso on bootstrap samples and selects variables
// that appear in a specified fraction of replications
// ============================================================================

class PenPPML_BootstrapLasso {
    // Data
    real colvector      y
    real matrix         X
    real matrix         fe_ids
    real colvector      w
    real colvector      cluster_id
    real scalar         n
    real scalar         p

    // Bootstrap settings
    real scalar         nboot           // Number of bootstrap replications
    real scalar         threshold       // Selection threshold (0.5 = 50%)

    // Results
    real matrix         selection_freq  // p x 1: fraction of times each var selected
    real colvector      selected        // Final selection indicator
    real scalar         n_selected

    // Methods
    void                new()
    void                set_data()
    void                run()
    void                print_summary()
}

void PenPPML_BootstrapLasso::new()
{
    nboot = 100
    threshold = 0.5
}

void PenPPML_BootstrapLasso::set_data(real colvector y_in, real matrix X_in,
                                       | real matrix fe_ids_in, real colvector w_in,
                                         real colvector cluster_id_in)
{
    y = y_in
    X = X_in
    n = rows(y)
    p = cols(X)

    if (args() >= 3) fe_ids = fe_ids_in
    if (args() >= 4) w = w_in
    if (args() >= 5) cluster_id = cluster_id_in
}

void PenPPML_BootstrapLasso::run()
{
    real scalar b, j
    real colvector boot_idx
    real colvector y_boot, w_boot
    real matrix X_boot, fe_boot
    class PenPPML_Plugin scalar plugin
    class PenPPML scalar M
    real matrix selections

    printf("{txt}Bootstrap Lasso: %g replications, threshold = %g\n",
           nboot, threshold)

    selections = J(p, nboot, 0)

    for (b = 1; b <= nboot; b++) {
        if (mod(b, 10) == 0) {
            printf(".")
            displayflush()
        }

        // Bootstrap sample (with replacement)
        if (rows(cluster_id) > 0) {
            // Cluster bootstrap
            boot_idx = cluster_bootstrap_idx(cluster_id)
        }
        else {
            // Simple bootstrap
            boot_idx = ceil(n * runiform(n, 1))
        }

        y_boot = y[boot_idx]
        X_boot = X[boot_idx, .]
        w_boot = w[boot_idx]
        if (rows(fe_ids) > 0) {
            fe_boot = fe_ids[boot_idx, .]
        }
        else {
            fe_boot = J(0, 0, .)
        }

        // Compute plugin weights for this sample
        plugin = PenPPML_Plugin()
        plugin.set_data(y_boot, X_boot, fe_boot, w_boot)
        plugin.compute_penalty_loadings()
        plugin.compute_lambda()

        // Fit penalized model
        M = PenPPML()
        M.set_data(y_boot, X_boot, J(n, 0, .), w_boot)
        if (rows(fe_boot) > 0) {
            M.set_fe(fe_boot)
        }
        M.set_penalty("lasso", plugin.lambda, 1, plugin.psi)
        M.set_options(1e-7, 50, 1, 0)
        M.solve()

        // Record selections
        selections[., b] = M.selected
    }

    printf("\n")

    // Compute selection frequency
    selection_freq = mean(selections')

    // Final selection based on threshold
    selected = (selection_freq :>= threshold)
    n_selected = sum(selected)
}

real colvector cluster_bootstrap_idx(real colvector cluster_id)
{
    real scalar n, n_clusters, g
    real colvector boot_clusters, boot_idx
    real colvector cluster_obs

    n = rows(cluster_id)
    n_clusters = max(cluster_id)

    // Sample clusters with replacement
    boot_clusters = ceil(n_clusters * runiform(n_clusters, 1))

    // Collect all observations from sampled clusters
    boot_idx = J(0, 1, .)
    for (g = 1; g <= n_clusters; g++) {
        cluster_obs = selectindex(cluster_id :== boot_clusters[g])
        boot_idx = boot_idx \ cluster_obs
    }

    return(boot_idx)
}

void PenPPML_BootstrapLasso::print_summary()
{
    printf("\n{txt}Bootstrap Lasso Results\n")
    printf("{hline 60}\n")
    printf("Bootstrap replications: %g\n", nboot)
    printf("Selection threshold: %g\n", threshold)
    printf("Variables selected: %g / %g\n", n_selected, p)
    printf("{hline 60}\n")
    printf("Selection frequency distribution:\n")
    printf("  [0.0-0.2]: %g variables\n", sum(selection_freq :< 0.2))
    printf("  [0.2-0.4]: %g variables\n", sum((selection_freq :>= 0.2) :& (selection_freq :< 0.4)))
    printf("  [0.4-0.6]: %g variables\n", sum((selection_freq :>= 0.4) :& (selection_freq :< 0.6)))
    printf("  [0.6-0.8]: %g variables\n", sum((selection_freq :>= 0.6) :& (selection_freq :< 0.8)))
    printf("  [0.8-1.0]: %g variables\n", sum(selection_freq :>= 0.8))
    printf("{hline 60}\n")
}

// ============================================================================
// ICEBERG LASSO
// Two-step procedure:
// 1. Run plugin lasso to get initial selection
// 2. Run individual lasso regressions on each selected variable
// ============================================================================

class PenPPML_IcebergLasso {
    // Data
    real colvector      y
    real matrix         X
    real matrix         fe_ids
    real colvector      w
    real scalar         n
    real scalar         p

    // Results
    real colvector      selected_step1  // First-step selection
    real colvector      selected_final  // Final selection
    real colvector      beta            // Final coefficients
    real scalar         n_selected

    // Methods
    void                new()
    void                run()
}

void PenPPML_IcebergLasso::new()
{
    // Empty constructor
}

void PenPPML_IcebergLasso::run(real colvector y_in, real matrix X_in,
                                | real matrix fe_ids_in, real colvector w_in)
{
    class PenPPML_Plugin scalar plugin
    class PenPPML scalar M
    real scalar j, n_step1
    real colvector sel_idx
    real matrix X_sel

    y = y_in
    X = X_in
    n = rows(y)
    p = cols(X)

    if (args() >= 3) fe_ids = fe_ids_in
    else fe_ids = J(0, 0, .)

    if (args() >= 4) w = w_in
    else w = J(n, 1, 1)

    printf("{txt}Iceberg Lasso: Step 1 - Plugin selection\n")

    // Step 1: Plugin lasso for initial selection
    plugin = PenPPML_Plugin()
    plugin.set_data(y, X, fe_ids, w)
    plugin.compute_penalty_loadings()
    plugin.compute_lambda()

    M = PenPPML()
    M.set_data(y, X, J(n, 0, .), w)
    if (rows(fe_ids) > 0) {
        M.set_fe(fe_ids)
    }
    M.set_penalty("lasso", plugin.lambda, 1, plugin.psi)
    M.set_options(1e-8, 100, 1, 0)
    M.solve()

    selected_step1 = M.selected
    n_step1 = sum(selected_step1)

    printf("Step 1 selected: %g variables\n", n_step1)

    if (n_step1 == 0) {
        selected_final = J(p, 1, 0)
        beta = J(p, 1, 0)
        n_selected = 0
        return
    }

    printf("{txt}Iceberg Lasso: Step 2 - Individual variable lasso\n")

    // Step 2: For each selected variable, run individual lasso
    // This refines the selection by checking if each variable
    // survives when considered individually
    sel_idx = selectindex(selected_step1)
    selected_final = J(p, 1, 0)

    for (j = 1; j <= n_step1; j++) {
        // Run lasso with only this variable
        X_sel = X[., sel_idx[j]]

        M = PenPPML()
        M.set_data(y, X_sel, J(n, 0, .), w)
        if (rows(fe_ids) > 0) {
            M.set_fe(fe_ids)
        }
        M.set_penalty("lasso", plugin.lambda * plugin.psi[sel_idx[j]], 1, J(1, 1, 1))
        M.set_options(1e-8, 50, 1, 0)
        M.solve()

        if (M.n_selected > 0) {
            selected_final[sel_idx[j]] = 1
        }
    }

    n_selected = sum(selected_final)

    // Compute final coefficients via unpenalized fit on selected
    if (n_selected > 0) {
        sel_idx = selectindex(selected_final)
        X_sel = X[., sel_idx]

        M = PenPPML()
        M.set_data(y, X_sel, J(n, 0, .), w)
        if (rows(fe_ids) > 0) {
            M.set_fe(fe_ids)
        }
        M.set_penalty("ridge", 0, 0, J(n_selected, 1, 1))
        M.set_options(1e-8, 100, 1, 0)
        M.solve()

        beta = J(p, 1, 0)
        for (j = 1; j <= n_selected; j++) {
            beta[sel_idx[j]] = M.beta[j]
        }
    }
    else {
        beta = J(p, 1, 0)
    }

    printf("Final selection: %g variables\n", n_selected)
}

// ============================================================================
// END OF PLUGIN LASSO ROUTINES
// ============================================================================

end
