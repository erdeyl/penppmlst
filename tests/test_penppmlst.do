/*******************************************************************************
* Test Suite for penppmlst
* Version 0.1.0
* Date: 03jan2026
*
* This file contains unit tests and integration tests for the penppmlst package.
* Run all tests with: do test_penppmlst.do
*******************************************************************************/

clear all
set more off
cap log close
log using test_penppmlst.log, replace

* Set path to source files
* After installation via net install, the files are in your ado path

* Load Mata source files
run "penppmlst_utils.mata"
run "penppmlst.mata"
run "penppmlst_cv.mata"
run "penppmlst_plugin.mata"

* Load ado files


di as txt _n "{hline 70}"
di as txt "PENPPML TEST SUITE"
di as txt "{hline 70}" _n

local n_tests = 0
local n_passed = 0
local n_failed = 0

********************************************************************************
* TEST 1: Soft-thresholding function
********************************************************************************
di as txt "TEST 1: Soft-thresholding function"

mata {
    // Test cases: (input, threshold, expected output)
    z = (5 \ -3 \ 2 \ 0.5 \ -0.5)
    gamma = (2 \ 2 \ 2 \ 2 \ 2)
    expected = (3 \ -1 \ 0 \ 0 \ 0)

    result = soft_threshold(z, gamma)

    passed = all(abs(result - expected) :< 1e-10)
    if (passed) {
        printf("  PASSED: Soft-thresholding gives correct results\n")
        st_local("test1", "passed")
    } else {
        printf("  FAILED: Soft-thresholding incorrect\n")
        result
        expected
        st_local("test1", "failed")
    }
}

if "`test1'" == "passed" {
    local ++n_passed
}
else {
    local ++n_failed
}
local ++n_tests

********************************************************************************
* TEST 2: Ridge solver
********************************************************************************
di as txt _n "TEST 2: Ridge regression solver"

mata {
    // Simple test case: y = 2*x1 + 3*x2 + noise
    n = 100
    X = rnormal(n, 2, 0, 1)
    y = X[., 1] * 2 + X[., 2] * 3 + rnormal(n, 1, 0, 0.1)
    w = J(n, 1, 1)

    // With lambda=0, should approximate OLS
    beta_ridge = ridge_solve(X, y, w, 0)

    // Check if close to true values
    error = max(abs(beta_ridge - (2 \ 3)))
    passed = (error < 0.5)

    if (passed) {
        printf("  PASSED: Ridge solver (lambda=0) approximates OLS (error = %g)\n", error)
        st_local("test2", "passed")
    } else {
        printf("  FAILED: Ridge solver error too large (%g)\n", error)
        st_local("test2", "failed")
    }
}

if "`test2'" == "passed" {
    local ++n_passed
}
else {
    local ++n_failed
}
local ++n_tests

********************************************************************************
* TEST 3: Coordinate descent for lasso
********************************************************************************
di as txt _n "TEST 3: Coordinate descent for lasso"

mata {
    // Test: sparse true model, high lambda should zero out noise variables
    n = 200
    p = 10

    // True model: only first 2 variables matter
    X = rnormal(n, p, 0, 1)
    beta_true = (3 \ -2 \ J(p-2, 1, 0))
    y = X * beta_true + rnormal(n, 1, 0, 0.5)
    w = J(n, 1, 1)
    psi = J(p, 1, 1)

    // High lambda should select only true variables (approximately)
    lambda = 5
    beta_lasso = coordinate_descent(X, y, w, lambda, psi, 1, 1e-7, 1000, 0)

    // Check sparsity: noise variables should be near zero
    noise_coefs = beta_lasso[3::p]
    max_noise = max(abs(noise_coefs))

    // Check signal: true variables should be non-zero (but shrunk)
    signal_coefs = beta_lasso[1::2]
    min_signal = min(abs(signal_coefs))

    passed = (max_noise < 0.1) & (min_signal > 0.5)

    if (passed) {
        printf("  PASSED: Lasso correctly sparsifies (noise max=%g, signal min=%g)\n",
               max_noise, min_signal)
        st_local("test3", "passed")
    } else {
        printf("  FAILED: Lasso sparsity incorrect\n")
        printf("  Noise max: %g, Signal min: %g\n", max_noise, min_signal)
        st_local("test3", "failed")
    }
}

if "`test3'" == "passed" {
    local ++n_passed
}
else {
    local ++n_failed
}
local ++n_tests

********************************************************************************
* TEST 4: Poisson deviance computation
********************************************************************************
di as txt _n "TEST 4: Poisson deviance computation"

mata {
    // Perfect fit should give deviance = 0
    y = (1 \ 2 \ 3 \ 4 \ 5)
    mu = y  // Perfect prediction

    dev = compute_deviance(y, mu)
    passed = (abs(dev) < 1e-10)

    if (passed) {
        printf("  PASSED: Deviance = 0 for perfect fit\n")
        st_local("test4a", "passed")
    } else {
        printf("  FAILED: Deviance should be 0 for y=mu, got %g\n", dev)
        st_local("test4a", "failed")
    }

    // Deviance should be positive for imperfect fit
    mu2 = (2 \ 2 \ 2 \ 2 \ 2)
    dev2 = compute_deviance(y, mu2)
    passed2 = (dev2 > 0)

    if (passed2) {
        printf("  PASSED: Deviance > 0 for imperfect fit (dev = %g)\n", dev2)
        st_local("test4b", "passed")
    } else {
        printf("  FAILED: Deviance should be positive, got %g\n", dev2)
        st_local("test4b", "failed")
    }
}

if "`test4a'" == "passed" & "`test4b'" == "passed" {
    local ++n_passed
}
else {
    local ++n_failed
}
local ++n_tests

********************************************************************************
* TEST 5: Lambda sequence generation
********************************************************************************
di as txt _n "TEST 5: Lambda sequence generation"

mata {
    n = 100
    p = 5
    X = rnormal(n, p, 0, 1)
    y = exp(X[., 1]) + rnormal(n, 1, 0, 0.1)
    w = J(n, 1, 1)
    psi = J(p, 1, 1)

    lambdas = generate_lambda_sequence(X, y, w, psi, 1, 50, 0.01)

    // Check properties
    nlam = rows(lambdas)
    is_decreasing = all(lambdas[2::nlam] :< lambdas[1::nlam-1])
    is_positive = all(lambdas :> 0)
    correct_count = (nlam == 50)

    passed = is_decreasing & is_positive & correct_count

    if (passed) {
        printf("  PASSED: Lambda sequence is valid (n=%g, max=%g, min=%g)\n",
               nlam, max(lambdas), min(lambdas))
        st_local("test5", "passed")
    } else {
        printf("  FAILED: Lambda sequence invalid\n")
        st_local("test5", "failed")
    }
}

if "`test5'" == "passed" {
    local ++n_passed
}
else {
    local ++n_failed
}
local ++n_tests

********************************************************************************
* TEST 6: PenPPML class - basic fit
********************************************************************************
di as txt _n "TEST 6: PenPPML class - basic estimation"

mata {
    // Generate Poisson data
    n = 200
    p = 5

    X = rnormal(n, p, 0, 1)
    eta_true = X[., 1] * 0.5 + X[., 2] * (-0.3)
    mu_true = exp(eta_true)
    y = rpoisson(1, mu_true)
    w = J(n, 1, 1)

    // Fit penalized model
    class PenPPML scalar M
    M.set_data(y, X, J(n, 0, .), w)
    M.set_penalty("lasso", 0.1, 1, J(p, 1, 1))
    M.set_options(1e-8, 100, 1, 0)
    converged = M.solve()

    // Checks
    check_converged = (converged == 1)
    check_deviance = (M.deviance < 1000)  // Reasonable deviance
    check_coefs = (rows(M.beta) == p)

    passed = check_converged & check_deviance & check_coefs

    if (passed) {
        printf("  PASSED: PenPPML converged (iter=%g, dev=%g, selected=%g)\n",
               M.iterations, M.deviance, M.n_selected)
        st_local("test6", "passed")
    } else {
        printf("  FAILED: PenPPML estimation failed\n")
        printf("  Converged: %g, Deviance: %g, Beta rows: %g\n",
               converged, M.deviance, rows(M.beta))
        st_local("test6", "failed")
    }
}

if "`test6'" == "passed" {
    local ++n_passed
}
else {
    local ++n_failed
}
local ++n_tests

********************************************************************************
* TEST 7: Ridge vs Lasso - Ridge should not zero coefficients
********************************************************************************
di as txt _n "TEST 7: Ridge vs Lasso behavior"

mata {
    n = 100
    p = 5
    X = rnormal(n, p, 0, 1)
    y = rpoisson(1, exp(0.5 :+ X * J(p, 1, 0.1)))
    w = J(n, 1, 1)

    // Lasso with moderate lambda
    class PenPPML scalar M_lasso
    M_lasso.set_data(y, X, J(n, 0, .), w)
    M_lasso.set_penalty("lasso", 1, 1, J(p, 1, 1))
    M_lasso.set_options(1e-8, 100, 1, 0)
    M_lasso.solve()

    // Ridge with same lambda
    class PenPPML scalar M_ridge
    M_ridge.set_data(y, X, J(n, 0, .), w)
    M_ridge.set_penalty("ridge", 1, 0, J(p, 1, 1))
    M_ridge.set_options(1e-8, 100, 1, 0)
    M_ridge.solve()

    // Ridge should have all non-zero, lasso may have zeros
    n_nonzero_lasso = sum(abs(M_lasso.beta) :> 1e-10)
    n_nonzero_ridge = sum(abs(M_ridge.beta) :> 1e-10)

    passed = (n_nonzero_ridge == p) & (n_nonzero_lasso <= p)

    if (passed) {
        printf("  PASSED: Ridge has %g non-zero, Lasso has %g non-zero\n",
               n_nonzero_ridge, n_nonzero_lasso)
        st_local("test7", "passed")
    } else {
        printf("  FAILED: Ridge should keep all coefficients non-zero\n")
        st_local("test7", "failed")
    }
}

if "`test7'" == "passed" {
    local ++n_passed
}
else {
    local ++n_failed
}
local ++n_tests

********************************************************************************
* TEST 8: Cross-validation class
********************************************************************************
di as txt _n "TEST 8: Cross-validation"

mata {
    n = 150
    p = 5
    X = rnormal(n, p, 0, 1)
    y = rpoisson(1, exp(X[., 1] * 0.3))
    w = J(n, 1, 1)

    class PenPPML_CV scalar CV
    CV.set_data(y, X)
    CV.set_cv_options(3, "deviance")  // 3-fold CV
    CV.set_penalty_options("lasso", 1, J(p, 1, 1), 10)  // 10 lambda values
    CV.run()

    // Checks
    has_results = (rows(CV.mean_score) == 10)
    lambda_positive = (CV.lambda_min > 0)
    idx_valid = (CV.idx_min >= 1) & (CV.idx_min <= 10)

    passed = has_results & lambda_positive & idx_valid

    if (passed) {
        printf("  PASSED: CV completed (lambda_min=%g, idx=%g)\n",
               CV.lambda_min, CV.idx_min)
        st_local("test8", "passed")
    } else {
        printf("  FAILED: CV results invalid\n")
        st_local("test8", "failed")
    }
}

if "`test8'" == "passed" {
    local ++n_passed
}
else {
    local ++n_failed
}
local ++n_tests

********************************************************************************
* TEST 9: Information criteria
********************************************************************************
di as txt _n "TEST 9: Information criteria (AIC, BIC, EBIC)"

mata {
    ll = -100
    k = 5
    n = 200
    p = 50

    aic = compute_aic(ll, k)
    bic = compute_bic(ll, k, n)
    ebic = compute_ebic(ll, k, n, p, 0.5)

    // AIC should be smallest (no sample size penalty)
    // BIC adds log(n) penalty
    // EBIC adds additional log(p) penalty
    correct_order = (aic < bic) & (bic < ebic)

    // Check formulas
    aic_correct = abs(aic - (-2*ll + 2*k)) < 1e-10
    bic_correct = abs(bic - (-2*ll + k*ln(n))) < 1e-10

    passed = correct_order & aic_correct & bic_correct

    if (passed) {
        printf("  PASSED: IC values correct (AIC=%g, BIC=%g, EBIC=%g)\n",
               aic, bic, ebic)
        st_local("test9", "passed")
    } else {
        printf("  FAILED: IC computation incorrect\n")
        st_local("test9", "failed")
    }
}

if "`test9'" == "passed" {
    local ++n_passed
}
else {
    local ++n_failed
}
local ++n_tests

********************************************************************************
* TEST 10: Fixed effects partialling (simple case)
********************************************************************************
di as txt _n "TEST 10: Fixed effects partialling"

mata {
    // Create data with group means
    n = 100
    n_groups = 10
    group_id = ceil((1::n) / (n/n_groups))

    // Variable with group effects
    group_means = (1::n_groups) :* 2
    x = J(n, 1, 0)
    for (i = 1; i <= n; i++) {
        x[i] = group_means[group_id[i]] + rnormal(1, 1, 0, 0.1)
    }

    w = J(n, 1, 1)

    // Create PenPPML and partial out
    class PenPPML scalar M
    M.n = n
    M.n_fe = 1
    M.fe_ids = group_id

    x_demean = M.partial_out_fe(x, w)

    // After demeaning, group means should be near zero
    max_group_mean = 0
    for (g = 1; g <= n_groups; g++) {
        group_mean = mean(select(x_demean, group_id :== g))
        if (abs(group_mean) > max_group_mean) {
            max_group_mean = abs(group_mean)
        }
    }

    passed = (max_group_mean < 0.01)

    if (passed) {
        printf("  PASSED: Group means after partialling < 0.01 (max=%g)\n",
               max_group_mean)
        st_local("test10", "passed")
    } else {
        printf("  FAILED: Group means not zeroed (max=%g)\n", max_group_mean)
        st_local("test10", "failed")
    }
}

if "`test10'" == "passed" {
    local ++n_passed
}
else {
    local ++n_failed
}
local ++n_tests

********************************************************************************
* SUMMARY
********************************************************************************
di as txt _n "{hline 70}"
di as txt "TEST SUMMARY"
di as txt "{hline 70}"
di as txt "Total tests:  " as res `n_tests'
di as txt "Passed:       " as res `n_passed' as txt " (" as res %4.1f 100*`n_passed'/`n_tests' as txt "%)"
di as txt "Failed:       " as res `n_failed' as txt " (" as res %4.1f 100*`n_failed'/`n_tests' as txt "%)"
di as txt "{hline 70}"

if `n_failed' == 0 {
    di as txt _n "{bf:ALL TESTS PASSED}"
}
else {
    di as err _n "{bf:SOME TESTS FAILED}"
}

log close
