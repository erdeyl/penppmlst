*! version 0.1.0  03jan2026
*! penppmlst_p: Predict command for penppmlst
*! Stata implementation by Erdey, László (2026)
*!   Faculty of Economics and Business, University of Debrecen, Hungary
*! Based on R penppmlst by Breinlich, Corradi, Rocha, Ruta, Santos Silva, Zylkin

program define penppmlst_p
    version 17.0

    // =========================================================================
    // SYNTAX
    // =========================================================================

    syntax newvarname [if] [in] , [    ///
        MU                             /// Predicted mean (exp(xb)), default
        XB                             /// Linear predictor
        ETa                            /// Same as xb
        RESiduals                      /// Response residuals (y - mu)
        DEViance                       /// Deviance residuals
        ANSCombe                       /// Anscombe residuals
        PEARson                        /// Pearson residuals
        SCOres                         /// Score (first derivative of log-likelihood)
        SELected                       /// Indicator for selected variables
        noOFFset                       /// Ignore offset
    ]

    // =========================================================================
    // CHECK ESTIMATION RESULTS
    // =========================================================================

    if "`e(cmd)'" != "penppmlst" {
        di as error "penppmlst_p only works after penppmlst estimation"
        exit 301
    }

    // =========================================================================
    // SET DEFAULT AND CHECK OPTIONS
    // =========================================================================

    local type_count = ("`mu'" != "") + ("`xb'" != "") + ("`eta'" != "") + ///
                       ("`residuals'" != "") + ("`deviance'" != "") + ///
                       ("`anscombe'" != "") + ("`pearson'" != "") + ///
                       ("`scores'" != "") + ("`selected'" != "")

    if `type_count' > 1 {
        di as error "Only one prediction type allowed"
        exit 198
    }

    if `type_count' == 0 {
        local mu "mu"
        di as txt "(option mu assumed; predicted mean)"
    }

    // Handle eta as alias for xb
    if "`eta'" != "" {
        local xb "xb"
    }

    // =========================================================================
    // MARK SAMPLE
    // =========================================================================

    marksample touse, novarlist
    qui replace `touse' = 0 if !e(sample)

    // =========================================================================
    // GET ESTIMATION INFORMATION
    // =========================================================================

    local depvar `e(depvar)'
    local indepvars `e(indepvars)'

    tempname b
    matrix `b' = e(b)

    local nvars : word count `indepvars'

    // =========================================================================
    // COMPUTE PREDICTIONS
    // =========================================================================

    tempvar xb_temp mu_temp

    // Compute linear predictor
    qui gen double `xb_temp' = 0 if `touse'

    local i = 0
    foreach var of local indepvars {
        local ++i
        qui replace `xb_temp' = `xb_temp' + `b'[1, `i'] * `var' if `touse'
    }

    // Add offset if present and not suppressed
    if "`e(offset)'" != "" & "`offset'" == "" {
        qui replace `xb_temp' = `xb_temp' + `e(offset)' if `touse'
    }

    // Compute mu
    qui gen double `mu_temp' = exp(`xb_temp') if `touse'
    qui replace `mu_temp' = max(`mu_temp', 1e-10) if `touse'
    qui replace `mu_temp' = min(`mu_temp', 1e10) if `touse'

    // =========================================================================
    // GENERATE REQUESTED PREDICTION
    // =========================================================================

    if "`mu'" != "" {
        // Predicted mean
        gen `typlist' `varlist' = `mu_temp' if `touse'
        label var `varlist' "Predicted mean"
    }

    else if "`xb'" != "" {
        // Linear predictor
        gen `typlist' `varlist' = `xb_temp' if `touse'
        label var `varlist' "Linear predictor"
    }

    else if "`residuals'" != "" {
        // Response residuals: y - mu
        gen `typlist' `varlist' = `depvar' - `mu_temp' if `touse'
        label var `varlist' "Response residuals"
    }

    else if "`deviance'" != "" {
        // Deviance residuals
        // d_i = sign(y - mu) * sqrt(2 * [y*log(y/mu) - (y - mu)])
        tempvar dev_term
        qui gen double `dev_term' = . if `touse'

        // Handle y = 0 case: 0 * log(0) = 0
        qui replace `dev_term' = 2 * (-(`depvar' - `mu_temp')) if `depvar' == 0 & `touse'
        qui replace `dev_term' = 2 * (`depvar' * ln(`depvar' / `mu_temp') - (`depvar' - `mu_temp')) ///
            if `depvar' > 0 & `touse'

        gen `typlist' `varlist' = sign(`depvar' - `mu_temp') * sqrt(max(`dev_term', 0)) if `touse'
        label var `varlist' "Deviance residuals"
    }

    else if "`anscombe'" != "" {
        // Anscombe residuals
        // a_i = (3/2) * (y^(2/3) - mu^(2/3)) / mu^(1/6)
        gen `typlist' `varlist' = 1.5 * (`depvar'^(2/3) - `mu_temp'^(2/3)) / `mu_temp'^(1/6) if `touse'
        label var `varlist' "Anscombe residuals"
    }

    else if "`pearson'" != "" {
        // Pearson residuals: (y - mu) / sqrt(mu)
        gen `typlist' `varlist' = (`depvar' - `mu_temp') / sqrt(`mu_temp') if `touse'
        label var `varlist' "Pearson residuals"
    }

    else if "`scores'" != "" {
        // Score: d(log L)/d(eta) = y - mu
        gen `typlist' `varlist' = `depvar' - `mu_temp' if `touse'
        label var `varlist' "Score"
    }

    else if "`selected'" != "" {
        // Selected variable indicator
        // This creates a variable indicating whether each covariate was selected

        // Parse selected variables from e(selected)
        local selected_vars `e(selected)'

        // For now, just report a message
        di as txt "Selected variables: `selected_vars'"
        di as txt "Number selected: " as res e(n_selected) as txt " / " as res `nvars'

        // Generate indicator for observations with non-zero predictions
        gen `typlist' `varlist' = (`xb_temp' != 0) if `touse'
        label var `varlist' "Has non-zero predicted value"
    }

end

// ============================================================================
// EXTENDED PREDICT FOR MARGINAL EFFECTS (FUTURE)
// ============================================================================

program define penppmlst_margins
    version 17.0

    syntax [varlist(default=none)] [if] [in] , [    ///
        DYdx(varlist)                  /// Marginal effects
        EYdx(varlist)                  /// Elasticities
        ATmeans                        /// Evaluate at means
        ]

    di as error "margins after penppmlst not yet implemented"
    di as txt "Use post option with penppmlst for valid inference"
    exit 198

end
