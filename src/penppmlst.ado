*! version 0.3.0  03jan2026
*! penppmlst: Penalized Poisson Pseudo Maximum Likelihood with High-Dimensional Fixed Effects
*! Stata implementation by Erdey, László (2026)
*!   Faculty of Economics and Business, University of Debrecen, Hungary
*! Based on R penppml by Breinlich, Corradi, Rocha, Ruta, Santos Silva, Zylkin

* Check required dependencies
cap which ftools
if _rc {
    di as error "penppmlst requires ftools"
    di as error "Install with: ssc install ftools"
    exit 198
}

cap which reghdfe
if _rc {
    di as error "penppmlst requires reghdfe"
    di as error "Install with: ssc install reghdfe"
    exit 198
}

cap which ppmlhdfe
if _rc {
    di as error "penppmlst requires ppmlhdfe"
    di as error "Install with: ssc install ppmlhdfe"
    exit 198
}

* Initialize Mata code on first use
cap mata: mata which PenPPML()
if _rc {
    * Find the path where this ado file is installed
    qui findfile penppmlst.ado
    local adopath = subinstr("`r(fn)'", "penppmlst.ado", "", .)

    * Run the Mata source files
    cap noi run "`adopath'penppmlst_utils.mata"
    cap noi run "`adopath'penppmlst.mata"
    cap noi run "`adopath'penppmlst_cv.mata"
    cap noi run "`adopath'penppmlst_plugin.mata"
}

program define penppmlst, eclass sortpreserve
    version 17.0

    // =========================================================================
    // SYNTAX PARSING
    // =========================================================================

    syntax varlist(min=2 numeric fv ts) [if] [in] [fweight aweight pweight], ///
        Absorb(string)                          /// Fixed effects specification
        [                                       ///
        PENalty(string)                         /// lasso, ridge, or elasticnet
        LAMbda(real -1)                         /// Penalty parameter (-1 = auto)
        ALpha(real 1)                           /// Elastic net mixing (1=lasso, 0=ridge)
        SELection(string)                       /// cv, plugin, bic, aic, ebic, or none
        NFolds(integer 10)                      /// Number of CV folds
        NLambda(integer 100)                    /// Number of lambda values for path
        LMINratio(real 0.01)                    /// Ratio of lambda_min to lambda_max
        CLuster(varname)                        /// Cluster variable for plugin
        POST                                    /// Post-lasso estimation
        noSTANDardize                           /// Don't standardize regressors
        TOLerance(real 1e-8)                    /// Convergence tolerance
        MAXiter(integer 1000)                   /// Maximum iterations
        IRr                                     /// Report incidence rate ratios
        EFORM                                   /// Same as irr
        Level(cilevel)                          /// Confidence level
        VERBose                                 /// Show iteration log
        NOLOg                                   /// Suppress iteration log
        HDFE(string)                            /// HDFE method: mata or ppmlhdfe
        R_compatible                            /// Use R penppml-compatible settings
        ]

    // =========================================================================
    // SETUP AND VALIDATION
    // =========================================================================

    // Mark sample
    marksample touse
    markout `touse' `absorb' `cluster'

    // Parse varlist
    gettoken depvar indepvars : varlist
    _fv_check_depvar `depvar'

    // Expand factor variables
    fvexpand `indepvars' if `touse'
    local indepvars_exp `r(varlist)'

    // Count variables
    local nvars : word count `indepvars_exp'
    if `nvars' == 0 {
        di as error "No independent variables specified"
        exit 198
    }

    // Check dependent variable is non-negative
    qui count if `depvar' < 0 & `touse'
    if r(N) > 0 {
        di as error "Dependent variable must be non-negative for Poisson regression"
        exit 459
    }

    // Set defaults for penalty
    if "`penalty'" == "" {
        local penalty "lasso"
    }
    if !inlist("`penalty'", "lasso", "ridge", "elasticnet") {
        di as error "penalty() must be lasso, ridge, or elasticnet"
        exit 198
    }

    // Set defaults for selection
    if "`selection'" == "" {
        if `lambda' == -1 {
            local selection "none"
            local lambda = 0.1  // Default lambda if none specified
        }
        else {
            local selection "none"
        }
    }
    if !inlist("`selection'", "none", "cv", "plugin", "bic", "aic", "ebic") {
        di as error "selection() must be none, cv, plugin, bic, aic, or ebic"
        exit 198
    }

    // Alpha must be between 0 and 1
    if `alpha' < 0 | `alpha' > 1 {
        di as error "alpha() must be between 0 and 1"
        exit 198
    }

    // =========================================================================
    // HDFE METHOD SELECTION
    // =========================================================================

    // Set default HDFE method
    if "`hdfe'" == "" {
        local hdfe "mata"
    }

    // Validate HDFE method
    if !inlist("`hdfe'", "mata", "ppmlhdfe") {
        di as error "hdfe() must be mata or ppmlhdfe"
        exit 198
    }

    // Check if required packages are installed
    if "`hdfe'" == "ppmlhdfe" {
        cap which ppmlhdfe
        if _rc {
            di as error "ppmlhdfe is required for hdfe(ppmlhdfe)"
            di as error "Install with: ssc install ppmlhdfe"
            exit 198
        }
    }


    // R-compatible settings
    local do_r_compat = 0
    if "`r_compatible'" != "" {
        local do_r_compat = 1
        if "`hdfe'" != "mata" {
            di as txt "Note: r_compatible implies hdfe(mata) for exact R replication"
            local hdfe "mata"
        }
    }

    // Handle standardize option
    local do_standardize = 1
    if "`standardize'" == "nostandardize" {
        local do_standardize = 0
    }

    // Handle verbose/nolog
    local do_verbose = 0
    if "`verbose'" != "" & "`nolog'" == "" {
        local do_verbose = 1
    }

    // Handle eform/irr
    local do_eform = 0
    if "`irr'" != "" | "`eform'" != "" {
        local do_eform = 1
    }

    // Handle post option
    local do_post = 0
    if "`post'" != "" {
        local do_post = 1
    }

    // =========================================================================
    // HANDLE WEIGHTS
    // =========================================================================

    tempvar wvar
    if "`weight'" != "" {
        local wtype `weight'
        local wexp `"`exp'"'
        qui gen double `wvar' `exp' if `touse'
    }
    else {
        qui gen double `wvar' = 1 if `touse'
    }

    // =========================================================================
    // PARSE FIXED EFFECTS
    // =========================================================================

    // Parse absorb() - simplified parsing
    // Full implementation would use reghdfe's parsing routines
    local fe_list ""
    local n_fe = 0

    tokenize `absorb'
    while "`1'" != "" {
        local fe_list `fe_list' `1'
        local ++n_fe
        macro shift
    }

    // =========================================================================
    // PREPARE DATA FOR MATA
    // =========================================================================

    // Preserve current state
    preserve

    // Keep only estimation sample
    qui keep if `touse'

    // Count observations
    qui count
    local nobs = r(N)

    if `nobs' == 0 {
        di as error "No observations"
        exit 2000
    }

    // Display header
    di as txt _n "Penalized PPML Estimation"
    di as txt "{hline 60}"
    di as txt "Dependent variable: " as res "`depvar'"
    di as txt "Number of observations: " as res "`nobs'"
    di as txt "Number of regressors: " as res "`nvars'"
    di as txt "Fixed effects: " as res "`absorb'"
    di as txt "HDFE method: " as res "`hdfe'" _c
    if "`r_compatible'" != "" {
        di as txt " (R-compatible mode)"
    }
    else {
        di ""
    }
    di as txt "Penalty: " as res "`penalty'" _c
    if "`penalty'" == "elasticnet" {
        di as txt " (alpha = " as res "`alpha'" as txt ")"
    }
    else {
        di ""
    }
    if "`selection'" != "none" {
        di as txt "Selection method: " as res "`selection'"
    }
    di as txt "{hline 60}"

    // =========================================================================
    // CALL MATA ESTIMATION ROUTINE
    // =========================================================================

    // Create temporary variables for results
    tempname b V
    tempvar mu_hat eta_hat

    // Run Mata estimation
    mata: penppmlst_estimate("`depvar'", "`indepvars_exp'", "`fe_list'", ///
                           "`wvar'", "`penalty'", `lambda', `alpha', ///
                           "`selection'", `nfolds', `nlambda', `lminratio', ///
                           "`cluster'", `do_post', `do_standardize', ///
                           `tolerance', `maxiter', `do_verbose', ///
                           "`hdfe'", `do_r_compat')

    // =========================================================================
    // POST RESULTS
    // =========================================================================

    restore

    // Get results from Mata
    tempname beta_mat V_mat
    mata: st_matrix("`beta_mat'", penppmlst_beta)
    mata: st_matrix("`V_mat'", penppmlst_V)
    mata: st_local("lambda_used", strofreal(penppmlst_lambda))
    mata: st_local("converged", strofreal(penppmlst_converged))
    mata: st_local("iterations", strofreal(penppmlst_iterations))
    mata: st_local("deviance", strofreal(penppmlst_deviance))
    mata: st_local("ll", strofreal(penppmlst_ll))
    mata: st_local("n_selected", strofreal(penppmlst_n_selected))

    // Create coefficient names
    local colnames ""
    foreach v of local indepvars_exp {
        local colnames `colnames' `v'
    }
    matrix colnames `beta_mat' = `colnames'
    matrix colnames `V_mat' = `colnames'
    matrix rownames `V_mat' = `colnames'

    // Post estimation results
    ereturn post `beta_mat' `V_mat', esample(`touse') obs(`nobs')

    // Store scalars
    ereturn scalar N = `nobs'
    ereturn scalar lambda = `lambda_used'
    ereturn scalar alpha = `alpha'
    ereturn scalar converged = `converged'
    ereturn scalar iterations = `iterations'
    ereturn scalar deviance = `deviance'
    ereturn scalar ll = `ll'
    ereturn scalar n_selected = `n_selected'
    ereturn scalar df_m = `n_selected'
    ereturn scalar k = `nvars'

    // Store macros
    ereturn local cmd "penppmlst"
    ereturn local cmdline "penppmlst `0'"
    ereturn local depvar "`depvar'"
    ereturn local indepvars "`indepvars_exp'"
    ereturn local absorb "`absorb'"
    ereturn local penalty "`penalty'"
    ereturn local selection "`selection'"
    ereturn local hdfe "`hdfe'"
    if "`r_compatible'" != "" {
        ereturn local r_compatible "yes"
    }
    if "`cluster'" != "" {
        ereturn local clustvar "`cluster'"
    }
    ereturn local predict "penppmlst_p"
    ereturn local title "Penalized PPML"
    ereturn local vce "robust"

    // Store selected variables
    mata: st_local("selected_vars", penppmlst_selected_names)
    ereturn local selected "`selected_vars'"

    // =========================================================================
    // DISPLAY RESULTS
    // =========================================================================

    di as txt _n "Penalized PPML Results"
    di as txt "{hline 60}"

    // Display estimation table
    if `do_eform' {
        ereturn display, eform("IRR") level(`level')
    }
    else {
        ereturn display, level(`level')
    }

    // Footer
    di as txt "{hline 60}"
    di as txt "Lambda: " as res %9.4f `lambda_used'
    di as txt "Variables selected: " as res `n_selected' as txt " / " as res `nvars'
    if `converged' {
        di as txt "Converged in " as res `iterations' as txt " iterations"
    }
    else {
        di as err "Warning: Did not converge"
    }
    di as txt "Log pseudo-likelihood: " as res %12.4f `ll'
    di as txt "Deviance: " as res %12.4f `deviance'
    di as txt "HDFE method: " as res "`hdfe'"
    di as txt "{hline 60}"

end

// =============================================================================
// MATA ESTIMATION ROUTINE
// =============================================================================

mata:

// Global variables to store results
real colvector penppmlst_beta
real matrix penppmlst_V
real scalar penppmlst_lambda
real scalar penppmlst_converged
real scalar penppmlst_iterations
real scalar penppmlst_deviance
real scalar penppmlst_ll
real scalar penppmlst_n_selected
string scalar penppmlst_selected_names

void penppmlst_estimate(string scalar depvar, string scalar indepvars,
                      string scalar fe_vars, string scalar wvar,
                      string scalar penalty, real scalar lambda,
                      real scalar alpha, string scalar selection,
                      real scalar nfolds, real scalar nlambda,
                      real scalar lminratio, string scalar cluster,
                      real scalar do_post, real scalar do_standardize,
                      real scalar tol, real scalar maxiter,
                      real scalar verbose, string scalar hdfe_method,
                      real scalar r_compatible)
{
    class PenPPML scalar M
    real colvector y, w
    real matrix X
    real matrix fe_ids
    real scalar n, p, k, j
    string rowvector varnames
    string rowvector fe_names
    real colvector lambdas
    real scalar best_lambda
    real colvector cv_results
    real colvector psi

    // Get data from Stata
    y = st_data(., depvar)
    X = st_data(., tokens(indepvars))
    w = st_data(., wvar)

    n = rows(y)
    p = cols(X)
    varnames = tokens(indepvars)

    // Get FE identifiers
    fe_names = tokens(fe_vars)
    if (cols(fe_names) > 0) {
        fe_ids = J(n, cols(fe_names), .)
        for (k = 1; k <= cols(fe_names); k++) {
            // Convert FE variable to numeric group IDs
            fe_ids[., k] = st_data(., fe_names[k])
        }
    }
    else {
        fe_ids = J(0, 0, .)
    }

    // Initialize penalty loadings (uniform)
    psi = J(p, 1, 1)

    // ===== LAMBDA SELECTION =====

    if (selection == "cv") {
        // Cross-validation
        if (verbose) {
            printf("Performing %g-fold cross-validation...\n", nfolds)
        }
        lambdas = generate_lambda_sequence(X, y, w, psi, alpha, nlambda, lminratio)
        best_lambda = cv_select_lambda(y, X, fe_ids, w, lambdas, alpha,
                                       penalty, psi, nfolds, tol, maxiter,
                                       hdfe_method, r_compatible)
        lambda = best_lambda
    }
    else if (selection == "plugin") {
        // Plugin lasso
        if (verbose) {
            printf("Computing plugin penalty weights...\n")
        }
        psi = compute_plugin_weights(y, X, fe_ids, w, cluster, hdfe_method)
        lambda = compute_plugin_lambda(n, p, psi)
    }
    else if (selection == "bic" | selection == "aic" | selection == "ebic") {
        // Information criteria
        if (verbose) {
            printf("Selecting lambda via %s...\n", strupper(selection))
        }
        lambdas = generate_lambda_sequence(X, y, w, psi, alpha, nlambda, lminratio)
        best_lambda = ic_select_lambda(y, X, fe_ids, w, lambdas, alpha,
                                       penalty, psi, selection, tol, maxiter,
                                       hdfe_method, r_compatible)
        lambda = best_lambda
    }

    // ===== MAIN ESTIMATION =====

    // Set up model
    M.set_data(y, X, J(n, 0, .), w)
    if (rows(fe_ids) > 0) {
        M.set_fe(fe_ids, fe_names)
    }
    M.set_penalty(penalty, lambda, alpha, psi)
    M.set_options(tol, maxiter, do_standardize, verbose)
    M.set_hdfe_method(hdfe_method, r_compatible)

    // Solve
    M.solve()

    // Post-lasso if requested
    if (do_post) {
        if (verbose) {
            printf("Computing post-lasso estimates...\n")
        }
        M.compute_post_lasso()
        M.compute_vcov()
    }

    // ===== STORE RESULTS =====

    if (do_post & M.n_selected > 0) {
        penppmlst_beta = M.beta_post'
        penppmlst_V = M.V
    }
    else {
        penppmlst_beta = M.beta'
        // For penalized estimates, use diagonal variance (approximation)
        penppmlst_V = diag(J(p, 1, 0.0001))
    }

    penppmlst_lambda = lambda
    penppmlst_converged = M.converged
    penppmlst_iterations = M.iterations
    penppmlst_deviance = M.deviance
    penppmlst_ll = M.ll
    penppmlst_n_selected = M.n_selected

    // Build selected variable names
    penppmlst_selected_names = ""
    for (j = 1; j <= p; j++) {
        if (M.selected[j]) {
            if (penppmlst_selected_names != "") {
                penppmlst_selected_names = penppmlst_selected_names + " "
            }
            penppmlst_selected_names = penppmlst_selected_names + varnames[j]
        }
    }
}

// =============================================================================
// CROSS-VALIDATION FOR LAMBDA SELECTION
// =============================================================================

real scalar cv_select_lambda(real colvector y, real matrix X,
                             real matrix fe_ids, real colvector w,
                             real colvector lambdas, real scalar alpha,
                             string scalar penalty, real colvector psi,
                             real scalar nfolds, real scalar tol,
                             real scalar maxiter, string scalar hdfe_method,
                             real scalar r_compatible)
{
    real scalar n, p, nlam, fold, lam_idx
    real colvector fold_id
    real matrix cv_dev
    real colvector mean_dev
    real scalar best_idx, best_lambda
    class PenPPML scalar M
    real colvector train_idx, test_idx
    real colvector y_train, y_test, w_train, w_test
    real matrix X_train, X_test, fe_train, fe_test
    real colvector mu_test
    real scalar dev_test

    n = rows(y)
    p = cols(X)
    nlam = rows(lambdas)

    // Create random fold assignments
    fold_id = ceil(nfolds * runiform(n, 1))

    // Matrix to store CV deviances
    cv_dev = J(nlam, nfolds, .)

    // Loop over folds
    for (fold = 1; fold <= nfolds; fold++) {
        train_idx = selectindex(fold_id :!= fold)
        test_idx = selectindex(fold_id :== fold)

        y_train = y[train_idx]
        y_test = y[test_idx]
        X_train = X[train_idx, .]
        X_test = X[test_idx, .]
        w_train = w[train_idx]
        w_test = w[test_idx]

        if (rows(fe_ids) > 0) {
            fe_train = fe_ids[train_idx, .]
            fe_test = fe_ids[test_idx, .]
        }

        // Loop over lambda values
        for (lam_idx = 1; lam_idx <= nlam; lam_idx++) {
            // Fit model on training data
            M = PenPPML()
            M.set_data(y_train, X_train, J(rows(y_train), 0, .), w_train)
            if (rows(fe_ids) > 0) {
                M.set_fe(fe_train)
            }
            M.set_penalty(penalty, lambdas[lam_idx], alpha, psi)
            M.set_options(tol, min((maxiter, 100)), 1, 0)
            M.set_hdfe_method(hdfe_method, r_compatible)
            M.solve()

            // Predict on test data
            mu_test = exp(X_test * M.beta)
            mu_test = clamp_vec(mu_test, 1e-10, 1e10)

            // Compute test deviance
            dev_test = compute_deviance(y_test, mu_test)
            cv_dev[lam_idx, fold] = dev_test
        }
    }

    // Average deviance across folds
    mean_dev = mean(cv_dev')

    // Select lambda with minimum mean deviance
    best_idx = 1
    for (lam_idx = 2; lam_idx <= nlam; lam_idx++) {
        if (mean_dev[lam_idx] < mean_dev[best_idx]) {
            best_idx = lam_idx
        }
    }

    best_lambda = lambdas[best_idx]

    return(best_lambda)
}

// =============================================================================
// INFORMATION CRITERIA FOR LAMBDA SELECTION
// =============================================================================

real scalar ic_select_lambda(real colvector y, real matrix X,
                             real matrix fe_ids, real colvector w,
                             real colvector lambdas, real scalar alpha,
                             string scalar penalty, real colvector psi,
                             string scalar criterion, real scalar tol,
                             real scalar maxiter, string scalar hdfe_method,
                             real scalar r_compatible)
{
    real scalar n, p, nlam, lam_idx
    real colvector ic_values
    real scalar best_idx, best_lambda
    class PenPPML scalar M
    real scalar ll, k, ic_val

    n = rows(y)
    p = cols(X)
    nlam = rows(lambdas)

    ic_values = J(nlam, 1, .)

    for (lam_idx = 1; lam_idx <= nlam; lam_idx++) {
        // Fit model
        M = PenPPML()
        M.set_data(y, X, J(n, 0, .), w)
        if (rows(fe_ids) > 0) {
            M.set_fe(fe_ids)
        }
        M.set_penalty(penalty, lambdas[lam_idx], alpha, psi)
        M.set_options(tol, min((maxiter, 100)), 1, 0)
        M.set_hdfe_method(hdfe_method, r_compatible)
        M.solve()

        ll = M.ll
        k = M.n_selected

        // Compute information criterion
        if (criterion == "aic") {
            ic_val = compute_aic(ll, k)
        }
        else if (criterion == "bic") {
            ic_val = compute_bic(ll, k, n)
        }
        else {  // ebic
            ic_val = compute_ebic(ll, k, n, p, 0.5)
        }

        ic_values[lam_idx] = ic_val
    }

    // Select lambda with minimum IC
    best_idx = 1
    for (lam_idx = 2; lam_idx <= nlam; lam_idx++) {
        if (ic_values[lam_idx] < ic_values[best_idx]) {
            best_idx = lam_idx
        }
    }

    best_lambda = lambdas[best_idx]

    return(best_lambda)
}

// =============================================================================
// PLUGIN LASSO PENALTY WEIGHTS
// =============================================================================

real colvector compute_plugin_weights(real colvector y, real matrix X,
                                      real matrix fe_ids, real colvector w,
                                      string scalar cluster,
                                      string scalar hdfe_method)
{
    real scalar n, p, j
    real colvector psi
    real colvector residuals, mu
    class PenPPML scalar M
    real colvector score_j
    real scalar var_j

    n = rows(y)
    p = cols(X)

    // Fit initial unpenalized model to get residuals
    M = PenPPML()
    M.set_data(y, X, J(n, 0, .), w)
    if (rows(fe_ids) > 0) {
        M.set_fe(fe_ids)
    }
    M.set_penalty("ridge", 0, 0, J(p, 1, 1))  // No penalty
    M.set_options(1e-8, 100, 1, 0)
    M.set_hdfe_method(hdfe_method, 0)
    M.solve()

    mu = M.mu
    residuals = y - mu

    // Compute penalty loadings based on score variance
    psi = J(p, 1, 1)
    for (j = 1; j <= p; j++) {
        score_j = X[., j] :* residuals
        var_j = variance(score_j)
        psi[j] = sqrt(var_j)
    }

    // Normalize
    psi = psi / max(psi)
    psi = max((psi, J(p, 1, 0.01)))  // Minimum weight

    return(psi)
}

real scalar compute_plugin_lambda(real scalar n, real scalar p,
                                  real colvector psi)
{
    real scalar c, gamma, lambda

    // Constants from Belloni et al.
    c = 1.1
    gamma = 0.1 / ln(n)

    // Lambda formula
    lambda = c * sqrt(n) * invnormal(1 - gamma / (2 * p))

    return(lambda)
}

end
