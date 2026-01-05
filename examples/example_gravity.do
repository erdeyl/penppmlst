/*******************************************************************************
* Example: Penalized PPML for Gravity Model Estimation
* Package: penppmlst
* Version: 0.5.0
*
* This example demonstrates using penppmlst to estimate a gravity model of
* international trade with many potential trade agreement provisions.
*******************************************************************************/

clear all
set more off

* ============================================================================
* SETUP
* ============================================================================

* After net install, penppmlst is available automatically
* net install penppmlst, from("https://raw.githubusercontent.com/erdeyl/penppmlst/main/src") replace


* ============================================================================
* GENERATE SIMULATED TRADE DATA
* ============================================================================

* Number of countries and years
local n_exp = 30       // Exporters
local n_imp = 30       // Importers
local n_year = 10      // Years
local n_prov = 50      // Trade agreement provisions

* Generate bilateral trade data
set seed 12345
set obs `=`n_exp' * `n_imp' * `n_year''

* Country and year identifiers
gen exporter = ceil(_n / (`n_imp' * `n_year'))
gen importer = mod(ceil(_n / `n_year') - 1, `n_imp') + 1
gen year = mod(_n - 1, `n_year') + 2010

* Remove same-country pairs
drop if exporter == importer

* Create pair identifier
egen pair = group(exporter importer)

* Generate gravity variables
gen distance = runiform(100, 10000)
gen ln_dist = ln(distance)
gen contig = (runiform() < 0.1)
gen comlang = (runiform() < 0.15)
gen colony = (runiform() < 0.05)

* Generate trade agreement indicator
gen pta = (runiform() < 0.3)

* Generate provision dummies (only active when PTA exists)
forvalues p = 1/`n_prov' {
    gen prov`p' = pta * (runiform() < 0.3)
}

* True data generating process:
* Only provisions 1-5 actually matter
gen ln_trade_star = 10 ///
    - 0.8 * ln_dist ///
    + 0.5 * contig ///
    + 0.3 * comlang ///
    + 0.2 * colony ///
    + 0.15 * prov1 ///
    + 0.12 * prov2 ///
    - 0.10 * prov3 ///
    + 0.08 * prov4 ///
    + 0.05 * prov5

* Add exporter-year and importer-year effects
bysort exporter year: gen exp_year_fe = rnormal(0, 1) if _n == 1
bysort exporter year: replace exp_year_fe = exp_year_fe[1]

bysort importer year: gen imp_year_fe = rnormal(0, 1) if _n == 1
bysort importer year: replace imp_year_fe = imp_year_fe[1]

replace ln_trade_star = ln_trade_star + exp_year_fe + imp_year_fe

* Generate Poisson trade flows
gen trade = rpoisson(exp(ln_trade_star))

* ============================================================================
* ESTIMATION EXAMPLES
* ============================================================================

di as txt _n "{hline 70}"
di as txt "EXAMPLE 1: Standard PPML (Unpenalized)"
di as txt "{hline 70}" _n

* First, estimate with ppmlhdfe for comparison (if installed)
cap which ppmlhdfe
if _rc == 0 {
    ppmlhdfe trade ln_dist contig comlang colony prov1-prov5, ///
        absorb(i.exporter#i.year i.importer#i.year) d
}

di as txt _n "{hline 70}"
di as txt "EXAMPLE 2: Penalized PPML with Lasso (Fixed Lambda)"
di as txt "{hline 70}" _n

* Estimate with lasso penalty and fixed lambda
penppmlst trade prov1-prov`n_prov', ///
    absorb(exporter importer year) ///
    penalty(lasso) lambda(0.1) ///
    verbose

* Display results
di as txt _n "Selected provisions:"
di as txt "`e(selected)'"

di as txt _n "{hline 70}"
di as txt "EXAMPLE 3: Penalized PPML with Cross-Validation"
di as txt "{hline 70}" _n

* Use cross-validation to select lambda
* Note: d(fe_contrib) stores FE contribution for predictions
penppmlst trade prov1-prov`n_prov', ///
    absorb(exporter importer year) ///
    penalty(lasso) ///
    selection(cv) nfolds(5) nlambda(20) ///
    post d(fe_contrib_cv)

* Show results
di as txt _n "Lambda selected by CV: " as res e(lambda)
di as txt "Variables selected: " as res e(n_selected) " / `n_prov'"
di as txt "Selected: `e(selected)'"

* Prediction (requires d() option in estimation)
predict mu_cv, mu
predict resid_cv, residuals

summarize mu_cv resid_cv fe_contrib_cv

di as txt _n "{hline 70}"
di as txt "EXAMPLE 4: Plugin Lasso"
di as txt "{hline 70}" _n

* Use plugin method with pair clustering
penppmlst trade prov1-prov`n_prov', ///
    absorb(exporter importer year) ///
    penalty(lasso) ///
    selection(plugin) cluster(pair) ///
    post

di as txt _n "Lambda (plugin): " as res e(lambda)
di as txt "Variables selected: " as res e(n_selected)
di as txt "Selected: `e(selected)'"

di as txt _n "{hline 70}"
di as txt "EXAMPLE 5: Ridge Regression"
di as txt "{hline 70}" _n

* Ridge regression (no variable selection, just shrinkage)
penppmlst trade prov1-prov`n_prov', ///
    absorb(exporter importer year) ///
    penalty(ridge) lambda(1)

di as txt "All coefficients non-zero (ridge shrinks but doesn't zero)"
di as txt "Variables with |coef| > 0.01:"
matrix b = e(b)
local nvars = colsof(b)
forvalues i = 1/`nvars' {
    if abs(b[1,`i']) > 0.01 {
        local vname : word `i' of `e(indepvars)'
        di as txt "  `vname': " as res b[1,`i']
    }
}

di as txt _n "{hline 70}"
di as txt "EXAMPLE 6: Elastic Net (alpha = 0.5)"
di as txt "{hline 70}" _n

* Elastic net combines lasso and ridge
penppmlst trade prov1-prov`n_prov', ///
    absorb(exporter importer year) ///
    penalty(elasticnet) alpha(0.5) ///
    selection(cv) nfolds(5)

di as txt "Elastic net (alpha=0.5) selected: " as res e(n_selected) " variables"

di as txt _n "{hline 70}"
di as txt "EXAMPLE 7: Including Unpenalized Controls"
di as txt "{hline 70}" _n

* Some variables should not be penalized (traditional gravity variables)
* Currently, all indepvars are penalized; this shows the typical workflow

* Step 1: Residualize trade on control variables
reghdfe trade ln_dist contig comlang colony, ///
    absorb(i.exporter#i.year i.importer#i.year) residuals(trade_resid)

* Step 2: Apply penalized PPML to residuals (simplified approach)
* Note: For proper implementation, integrate unpenalized controls directly

di as txt _n "{hline 70}"
di as txt "EXAMPLE 8: R-Compatible Mode"
di as txt "{hline 70}" _n

* Use r_compatible for cross-platform reproducibility with R penppml
penppmlst trade prov1-prov`n_prov', ///
    absorb(exporter importer year) ///
    selection(plugin) cluster(pair) ///
    r_compatible d(fe_contrib_r)

di as txt "R-compatible mode uses:"
di as txt "  - R-style mu bounds [1e-190, 1e190]"
di as txt "  - R-style deviance computation"
di as txt "  - Pure Mata HDFE (matching R's collapse::fhdwithin)"

di as txt _n "{hline 70}"
di as txt "SUMMARY"
di as txt "{hline 70}"
di as txt "True provisions with effects: prov1-prov5"
di as txt "Recovery depends on signal strength and sample size"
di as txt "{hline 70}"

* ============================================================================
* CLEANUP
* ============================================================================

* Save results if needed
* save "gravity_results.dta", replace

di as txt _n "Example completed successfully!"
