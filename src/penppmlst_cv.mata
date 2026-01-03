*! version 0.1.0  03jan2026
*! Cross-validation routines for penppmlst package
*! Stata implementation by Erdey, László (2026)
*!   Faculty of Economics and Business, University of Debrecen, Hungary
*! Based on R penppml by Breinlich, Corradi, Rocha, Ruta, Santos Silva, Zylkin

version 17.0
mata:
mata set matastrict on

// ============================================================================
// K-FOLD CROSS-VALIDATION CLASS
// ============================================================================

class PenPPML_CV {
    // Data
    real colvector      y
    real matrix         X
    real matrix         fe_ids
    real colvector      w
    real scalar         n
    real scalar         p

    // CV settings
    real scalar         nfolds
    real colvector      fold_id
    string scalar       criterion       // "deviance", "mse", "mae"

    // Penalty settings
    string scalar       penalty
    real scalar         alpha
    real colvector      psi
    real colvector      lambdas
    real scalar         nlambda

    // Algorithm settings
    real scalar         tol
    real scalar         maxiter
    real scalar         standardize

    // Results
    real matrix         cv_scores       // nlambda x nfolds matrix of scores
    real colvector      mean_score      // Mean score across folds
    real colvector      se_score        // Standard error of scores
    real scalar         lambda_min      // Lambda with minimum mean score
    real scalar         lambda_1se      // Lambda within 1 SE of minimum
    real scalar         idx_min         // Index of lambda_min
    real scalar         idx_1se         // Index of lambda_1se
    real matrix         beta_path       // Coefficient path (nlambda x p)
    real colvector      n_selected_path // Number selected at each lambda

    // Methods
    void                new()
    void                set_data()
    void                set_cv_options()
    void                set_penalty_options()
    void                create_folds()
    void                set_folds()
    void                run()
    void                find_optimal_lambda()
    void                compute_full_path()
    void                print_summary()
}

// ============================================================================
// CONSTRUCTOR
// ============================================================================

void PenPPML_CV::new()
{
    nfolds = 10
    criterion = "deviance"
    penalty = "lasso"
    alpha = 1
    nlambda = 100
    tol = 1e-8
    maxiter = 100
    standardize = 1
}

// ============================================================================
// SET DATA
// ============================================================================

void PenPPML_CV::set_data(real colvector y_in, real matrix X_in,
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
// SET CV OPTIONS
// ============================================================================

void PenPPML_CV::set_cv_options(| real scalar nfolds_in, string scalar criterion_in)
{
    if (args() >= 1 & nfolds_in != .) nfolds = nfolds_in
    if (args() >= 2 & criterion_in != "") criterion = criterion_in
}

// ============================================================================
// SET PENALTY OPTIONS
// ============================================================================

void PenPPML_CV::set_penalty_options(string scalar penalty_in, real scalar alpha_in,
                                     | real colvector psi_in, real scalar nlambda_in,
                                       real colvector lambdas_in)
{
    penalty = penalty_in
    alpha = alpha_in

    if (args() >= 3 & rows(psi_in) > 0) {
        psi = psi_in
    }
    else {
        psi = J(p, 1, 1)
    }

    if (args() >= 4 & nlambda_in != .) nlambda = nlambda_in

    if (args() >= 5 & rows(lambdas_in) > 0) {
        lambdas = lambdas_in
        nlambda = rows(lambdas)
    }
    else {
        // Generate lambda sequence
        lambdas = generate_lambda_sequence(X, y, w, psi, alpha, nlambda, 0.01)
    }
}

// ============================================================================
// CREATE RANDOM FOLD ASSIGNMENTS
// ============================================================================

void PenPPML_CV::create_folds()
{
    real colvector perm
    real scalar i, fold_size, start_idx, end_idx

    // Random permutation for stratified assignment
    perm = jumble(1::n)

    fold_id = J(n, 1, .)
    fold_size = floor(n / nfolds)

    for (i = 1; i <= nfolds; i++) {
        start_idx = (i - 1) * fold_size + 1
        if (i == nfolds) {
            end_idx = n
        }
        else {
            end_idx = i * fold_size
        }
        fold_id[perm[start_idx::end_idx]] = J(end_idx - start_idx + 1, 1, i)
    }
}

// ============================================================================
// SET USER-DEFINED FOLDS
// ============================================================================

void PenPPML_CV::set_folds(real colvector fold_id_in)
{
    fold_id = fold_id_in
    nfolds = max(fold_id)
}

// ============================================================================
// RUN CROSS-VALIDATION
// ============================================================================

void PenPPML_CV::run()
{
    real scalar fold, lam_idx
    real colvector train_idx, test_idx
    real colvector y_train, y_test, w_train, w_test
    real matrix X_train, X_test, fe_train, fe_test
    class PenPPML scalar M
    real colvector mu_test, eta_test
    real scalar score

    // Create folds if not already set
    if (rows(fold_id) == 0) {
        create_folds()
    }

    // Initialize results matrix
    cv_scores = J(nlambda, nfolds, .)

    printf("{txt}Cross-validation: %g folds, %g lambda values\n", nfolds, nlambda)
    printf("{txt}Fold: ")

    // Loop over folds
    for (fold = 1; fold <= nfolds; fold++) {
        printf("%g ", fold)
        displayflush()

        // Split data
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
        else {
            fe_train = J(0, 0, .)
            fe_test = J(0, 0, .)
        }

        // Loop over lambda values (warm start)
        for (lam_idx = 1; lam_idx <= nlambda; lam_idx++) {

            // Fit model on training data
            M = PenPPML()
            M.set_data(y_train, X_train, J(rows(y_train), 0, .), w_train)
            if (rows(fe_train) > 0) {
                M.set_fe(fe_train)
            }
            M.set_penalty(penalty, lambdas[lam_idx], alpha, psi)
            M.set_options(tol, maxiter, standardize, 0)
            M.solve()

            // Predict on test data
            eta_test = X_test * M.beta
            mu_test = exp(eta_test)
            mu_test = clamp_vec(mu_test, 1e-10, 1e10)

            // Compute score
            if (criterion == "deviance") {
                score = compute_deviance(y_test, mu_test) / rows(y_test)
            }
            else if (criterion == "mse") {
                score = mean((y_test - mu_test):^2)
            }
            else {  // mae
                score = mean(abs(y_test - mu_test))
            }

            cv_scores[lam_idx, fold] = score
        }
    }

    printf("\n")

    // Compute summary statistics
    mean_score = mean(cv_scores')'
    se_score = J(nlambda, 1, .)
    for (lam_idx = 1; lam_idx <= nlambda; lam_idx++) {
        se_score[lam_idx] = sqrt(variance(cv_scores[lam_idx, .]')) / sqrt(nfolds)
    }

    // Find optimal lambda
    find_optimal_lambda()
}

// ============================================================================
// FIND OPTIMAL LAMBDA (MIN AND 1SE RULE)
// ============================================================================

void PenPPML_CV::find_optimal_lambda()
{
    real scalar lam_idx
    real scalar min_score, threshold

    // Find minimum
    idx_min = 1
    min_score = mean_score[1]
    for (lam_idx = 2; lam_idx <= nlambda; lam_idx++) {
        if (mean_score[lam_idx] < min_score) {
            min_score = mean_score[lam_idx]
            idx_min = lam_idx
        }
    }
    lambda_min = lambdas[idx_min]

    // 1SE rule: largest lambda within 1 SE of minimum
    threshold = min_score + se_score[idx_min]
    idx_1se = idx_min
    for (lam_idx = 1; lam_idx < idx_min; lam_idx++) {
        if (mean_score[lam_idx] <= threshold) {
            idx_1se = lam_idx
            break
        }
    }
    lambda_1se = lambdas[idx_1se]
}

// ============================================================================
// COMPUTE FULL REGULARIZATION PATH
// ============================================================================

void PenPPML_CV::compute_full_path()
{
    real scalar lam_idx
    class PenPPML scalar M

    beta_path = J(nlambda, p, .)
    n_selected_path = J(nlambda, 1, .)

    printf("{txt}Computing regularization path...\n")

    for (lam_idx = 1; lam_idx <= nlambda; lam_idx++) {
        M = PenPPML()
        M.set_data(y, X, J(n, 0, .), w)
        if (rows(fe_ids) > 0) {
            M.set_fe(fe_ids)
        }
        M.set_penalty(penalty, lambdas[lam_idx], alpha, psi)
        M.set_options(tol, maxiter, standardize, 0)
        M.solve()

        beta_path[lam_idx, .] = M.beta'
        n_selected_path[lam_idx] = M.n_selected
    }
}

// ============================================================================
// PRINT SUMMARY
// ============================================================================

void PenPPML_CV::print_summary()
{
    printf("\n{txt}Cross-Validation Results\n")
    printf("{hline 60}\n")
    printf("Number of folds: %g\n", nfolds)
    printf("Number of lambda values: %g\n", nlambda)
    printf("Criterion: %s\n", criterion)
    printf("{hline 60}\n")
    printf("Lambda (min):  %12.6f  (index %g)\n", lambda_min, idx_min)
    printf("Lambda (1se):  %12.6f  (index %g)\n", lambda_1se, idx_1se)
    printf("Mean CV score at min: %12.6f (SE: %8.6f)\n",
           mean_score[idx_min], se_score[idx_min])
    printf("{hline 60}\n")
}

// ============================================================================
// ROLLING CROSS-VALIDATION FOR TIME SERIES
// ============================================================================

class PenPPML_RollingCV {
    // Data
    real colvector      y
    real matrix         X
    real matrix         fe_ids
    real colvector      w
    real scalar         n
    real scalar         p

    // Rolling CV settings
    real scalar         origin          // First training window size
    real scalar         horizon         // Forecast horizon
    real scalar         step            // Step size between windows

    // Penalty settings
    string scalar       penalty
    real scalar         alpha
    real colvector      psi
    real colvector      lambdas
    real scalar         nlambda

    // Algorithm settings
    real scalar         tol
    real scalar         maxiter

    // Results
    real matrix         cv_scores       // nlambda x n_windows
    real colvector      mean_score
    real scalar         lambda_min
    real scalar         n_windows

    // Methods
    void                new()
    void                set_data()
    void                run()
}

void PenPPML_RollingCV::new()
{
    origin = 0
    horizon = 1
    step = 1
    penalty = "lasso"
    alpha = 1
    nlambda = 100
    tol = 1e-8
    maxiter = 100
}

void PenPPML_RollingCV::set_data(real colvector y_in, real matrix X_in,
                                  | real matrix fe_ids_in, real colvector w_in)
{
    y = y_in
    X = X_in
    n = rows(y)
    p = cols(X)

    if (args() >= 3 & rows(fe_ids_in) > 0) {
        fe_ids = fe_ids_in
    }

    if (args() >= 4 & rows(w_in) > 0) {
        w = w_in
    }
    else {
        w = J(n, 1, 1)
    }

    // Default origin: 50% of data
    if (origin == 0) {
        origin = floor(n / 2)
    }
}

void PenPPML_RollingCV::run()
{
    real scalar window, lam_idx, train_end, test_start, test_end
    real colvector train_idx, test_idx
    real colvector y_train, y_test, w_train
    real matrix X_train, X_test
    class PenPPML scalar M
    real colvector mu_test
    real scalar score

    // Count windows
    n_windows = floor((n - origin - horizon) / step) + 1
    if (n_windows < 1) {
        printf("{err}Not enough observations for rolling CV\n")
        return
    }

    // Generate lambda sequence if not provided
    if (rows(lambdas) == 0) {
        psi = J(p, 1, 1)
        lambdas = generate_lambda_sequence(X, y, w, psi, alpha, nlambda, 0.01)
    }

    cv_scores = J(nlambda, n_windows, .)

    printf("{txt}Rolling CV: %g windows, origin=%g, horizon=%g\n",
           n_windows, origin, horizon)

    for (window = 1; window <= n_windows; window++) {
        train_end = origin + (window - 1) * step
        test_start = train_end + 1
        test_end = min((test_start + horizon - 1, n))

        train_idx = 1::train_end
        test_idx = test_start::test_end

        y_train = y[train_idx]
        y_test = y[test_idx]
        X_train = X[train_idx, .]
        X_test = X[test_idx, .]
        w_train = w[train_idx]

        for (lam_idx = 1; lam_idx <= nlambda; lam_idx++) {
            M = PenPPML()
            M.set_data(y_train, X_train, J(rows(y_train), 0, .), w_train)
            M.set_penalty(penalty, lambdas[lam_idx], alpha, psi)
            M.set_options(tol, maxiter, 1, 0)
            M.solve()

            mu_test = exp(X_test * M.beta)
            mu_test = clamp_vec(mu_test, 1e-10, 1e10)
            score = compute_deviance(y_test, mu_test) / rows(y_test)
            cv_scores[lam_idx, window] = score
        }
    }

    mean_score = mean(cv_scores')'

    // Find minimum
    lambda_min = lambdas[1]
    for (lam_idx = 2; lam_idx <= nlambda; lam_idx++) {
        if (mean_score[lam_idx] < mean_score[1]) {
            lambda_min = lambdas[lam_idx]
        }
    }
}

// ============================================================================
// END OF CROSS-VALIDATION ROUTINES
// ============================================================================

end
