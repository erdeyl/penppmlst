*! version 0.5.0  05jan2026
*! penppmlst_p: Predict command for penppmlst
*! Stata implementation by Erdey, Laszlo (2026)
*!   Faculty of Economics and Business, University of Debrecen, Hungary
*! Based on R penppml by Breinlich, Corradi, Rocha, Ruta, Santos Silva, Zylkin

program define penppmlst_p
    version 17.0

    syntax newvarname [if] [in] , [    ///
        MU                             /// Predicted mean (exp(xb+d)), default
        XB                             /// Linear predictor (without FE)
        XBD                            /// Linear predictor including FE
        ETa                            /// Same as xbd
        D                              /// Fixed effect contribution only
        RESiduals                      /// Response residuals (y - mu)
        DEViance                       /// Deviance residuals
        ANSCombe                       /// Anscombe residuals
        PEARson                        /// Pearson residuals
        SCOres                         /// Score (first derivative of log-likelihood)
        SELected                       /// Indicator for selected variables
        noOFFset                       /// Ignore offset
    ]

    if "`e(cmd)'" != "penppmlst" {
        di as error "penppmlst_p only works after penppmlst estimation"
        exit 301
    }

    local type_count = ("`mu'" != "") + ("`xb'" != "") + ("`xbd'" != "") + ///
                       ("`eta'" != "") + ("`d'" != "") + ///
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

    if "`eta'" != "" {
        local xbd "xbd"
    }

    local needs_fe = ("`xbd'" != "" | "`mu'" != "" | "`d'" != "" | ///
                      "`residuals'" != "" | "`deviance'" != "" | ///
                      "`anscombe'" != "" | "`pearson'" != "" | "`scores'" != "")

    if `needs_fe' & "`e(d)'" == "" & "`e(absorb)'" != "" & "`e(absorb)'" != "_cons" {
        di as error "predict requires the -d()- option in penppmlst when FEs are used"
        di as txt "Re-run penppmlst with d(varname) option to store FE contribution"
        exit 198
    }

    marksample touse, novarlist
    qui replace `touse' = 0 if !e(sample)

    local depvar `e(depvar)'
    local indepvars `e(indepvars)'

    tempname b
    matrix `b' = e(b)

    local nvars : word count `indepvars'

    tempvar xb_temp d_temp mu_temp

    qui gen double `xb_temp' = 0 if `touse'

    local i = 0
    foreach var of local indepvars {
        local ++i
        qui replace `xb_temp' = `xb_temp' + `b'[1, `i'] * `var' if `touse'
    }

    if "`e(offset)'" != "" & "`offset'" == "" {
        qui replace `xb_temp' = `xb_temp' + `e(offset)' if `touse'
    }

    if "`e(d)'" != "" {
        conf var `e(d)', exact
        qui gen double `d_temp' = `e(d)' if `touse'
    }
    else {
        qui gen double `d_temp' = 0 if `touse'
    }

    qui gen double `mu_temp' = exp(`xb_temp' + `d_temp') if `touse'
    qui replace `mu_temp' = max(`mu_temp', 1e-10) if `touse'
    qui replace `mu_temp' = min(`mu_temp', 1e10) if `touse'

    if "`mu'" != "" {
        gen `typlist' `varlist' = `mu_temp' if `touse'
        label var `varlist' "Predicted mean"
    }
    else if "`xb'" != "" {
        gen `typlist' `varlist' = `xb_temp' if `touse'
        label var `varlist' "Linear predictor (xb)"
    }
    else if "`xbd'" != "" {
        gen `typlist' `varlist' = `xb_temp' + `d_temp' if `touse'
        label var `varlist' "Linear predictor (xb + d)"
    }
    else if "`d'" != "" {
        gen `typlist' `varlist' = `d_temp' if `touse'
        label var `varlist' "Fixed effect contribution (d)"
    }
    else if "`residuals'" != "" {
        gen `typlist' `varlist' = `depvar' - `mu_temp' if `touse'
        label var `varlist' "Response residuals"
    }
    else if "`deviance'" != "" {
        gen `typlist' `varlist' = 2 * cond(`depvar' > 0, ///
            `mu_temp' - `depvar' + `depvar' * ln(`depvar' / `mu_temp'), ///
            `mu_temp') if `touse'
        label var `varlist' "Deviance residuals"
    }
    else if "`anscombe'" != "" {
        gen `typlist' `varlist' = 1.5 * (`depvar'^(2/3) - `mu_temp'^(2/3)) / `mu_temp'^(1/6) if `touse'
        label var `varlist' "Anscombe residuals"
    }
    else if "`pearson'" != "" {
        gen `typlist' `varlist' = (`depvar' - `mu_temp') / sqrt(`mu_temp') if `touse'
        label var `varlist' "Pearson residuals"
    }
    else if "`scores'" != "" {
        gen `typlist' `varlist' = `depvar' - `mu_temp' if `touse'
        label var `varlist' "Score"
    }
    else if "`selected'" != "" {
        local selected_vars `e(selected)'
        di as txt "Selected variables: `selected_vars'"
        di as txt "Number selected: " as res e(n_selected) as txt " / " as res `nvars'
        gen `typlist' `varlist' = (`xb_temp' != 0) if `touse'
        label var `varlist' "Has non-zero predicted value"
    }
end

program define penppmlst_margins
    version 17.0
    syntax [varlist(default=none)] [if] [in] , [DYdx(varlist) EYdx(varlist) ATmeans]
    di as error "margins after penppmlst not yet implemented"
    di as txt "Use post option with penppmlst for valid inference"
    exit 198
end
