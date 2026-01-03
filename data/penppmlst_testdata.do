*! penppmlst_testdata.do
*! Creates test datasets for penppmlst package validation
*! Based on the structure of R penppml's trade dataset

* Create simulated trade panel data similar to penppml's trade dataset
* This is a synthetic dataset for testing purposes

clear all
set seed 12345

* Parameters
local n_exp = 50      // Number of exporters
local n_imp = 50      // Number of importers
local n_time = 5      // Number of time periods
local n_provisions = 16  // Number of provision dummies

* Generate panel structure
set obs `=`n_exp' * `n_imp' * `n_time''

* Create identifiers
gen exp_id = mod(_n - 1, `n_exp') + 1
gen imp_id = mod(floor((_n - 1) / `n_exp'), `n_imp') + 1
gen time = mod(floor((_n - 1) / (`n_exp' * `n_imp')), `n_time') + 2012

* Create string country codes (3-letter)
gen str3 exp = ""
gen str3 imp = ""
forvalues i = 1/`n_exp' {
    local c1 = char(65 + mod(`i'-1, 26))
    local c2 = char(65 + mod(floor((`i'-1)/26), 26))
    local c3 = char(65 + mod(floor((`i'-1)/676), 26))
    replace exp = "`c1'`c2'`c3'" if exp_id == `i'
}
forvalues i = 1/`n_imp' {
    local c1 = char(65 + mod(`i'-1, 26))
    local c2 = char(65 + mod(floor((`i'-1)/26), 26))
    local c3 = char(65 + mod(floor((`i'-1)/676), 26))
    replace imp = "`c1'`c2'`c3'" if imp_id == `i'
}

* Create pair identifier
gen pair_id = exp_id * 1000 + imp_id
egen pair = group(exp imp)

* Generate fixed effect contributions (to create realistic trade flows)
gen fe_exp = rnormal(10, 2) if _n <= `n_exp'
gen fe_imp = rnormal(10, 2) if _n <= `n_imp'
bysort exp_id: replace fe_exp = fe_exp[1]
bysort imp_id: replace fe_imp = fe_imp[1]

* Generate distance-like bilateral component
gen bilateral = abs(exp_id - imp_id) / `n_exp' * 2 + rnormal(0, 0.5)

* Generate trade agreement indicator (about 20% of pairs)
gen fta = runiform() < 0.2

* Generate provision dummies (conditional on FTA)
forvalues p = 1/`n_provisions' {
    gen prov_`p' = (fta == 1) * (runiform() < 0.3 + 0.4 * runiform())
}

* Create realistic export values using gravity-like structure
gen log_export = fe_exp + fe_imp - bilateral + ///
                 0.5 * fta + ///
                 0.1 * prov_1 + 0.15 * prov_2 - 0.05 * prov_3 + ///
                 0.2 * prov_4 + 0.08 * prov_5 + 0.12 * prov_6 + ///
                 rnormal(0, 1)

* Convert to levels (with many zeros for PPML)
gen export = exp(log_export) * (runiform() > 0.3)
replace export = 0 if exp == imp  // No self-trade
replace export = round(export * 1000000, 1)  // Scale to millions

* Add some noise provisions that have no effect
forvalues p = 7/`n_provisions' {
    replace prov_`p' = (fta == 1) * (runiform() < 0.25)
}

* Create cluster variable
gen cluster = pair

* Label variables
label var exp "Exporter (3-letter code)"
label var imp "Importer (3-letter code)"
label var time "Year"
label var export "Merchandise exports (USD)"
label var fta "Free Trade Agreement indicator"
label var pair "Pair identifier"
label var cluster "Cluster for robust SEs"

forvalues p = 1/`n_provisions' {
    label var prov_`p' "Trade provision `p'"
}

* Drop auxiliary variables
drop exp_id imp_id pair_id fe_exp fe_imp bilateral log_export

* Order variables
order exp imp time export fta pair cluster prov_*

* Summary statistics
describe
sum export, detail
tab fta

* Save dataset
compress
save "trade_test.dta", replace

* Create smaller Americas subset for quick testing
preserve
keep if inlist(substr(exp, 1, 1), "A", "B", "C") & inlist(substr(imp, 1, 1), "A", "B", "C")
save "trade_test_americas.dta", replace
restore

* Create countries lookup table
clear
set obs 50
gen str3 iso = ""
gen str50 name = ""
gen str20 region = ""
gen str30 subregion = ""

forvalues i = 1/50 {
    local c1 = char(65 + mod(`i'-1, 26))
    local c2 = char(65 + mod(floor((`i'-1)/26), 26))
    local c3 = char(65 + mod(floor((`i'-1)/676), 26))
    replace iso = "`c1'`c2'`c3'" in `i'
    replace name = "Country `i'" in `i'

    * Assign regions
    if `i' <= 10 {
        replace region = "Americas" in `i'
        replace subregion = "North America" in `i'
    }
    else if `i' <= 20 {
        replace region = "Americas" in `i'
        replace subregion = "South America" in `i'
    }
    else if `i' <= 30 {
        replace region = "Europe" in `i'
        replace subregion = "Western Europe" in `i'
    }
    else if `i' <= 40 {
        replace region = "Asia" in `i'
        replace subregion = "Eastern Asia" in `i'
    }
    else {
        replace region = "Africa" in `i'
        replace subregion = "Sub-Saharan Africa" in `i'
    }
}

label var iso "ISO 3166 code"
label var name "Country name"
label var region "Continent"
label var subregion "Sub-region"

save "countries_test.dta", replace

di _n "Test datasets created successfully:"
di "  - trade_test.dta (full panel)"
di "  - trade_test_americas.dta (subset for quick tests)"
di "  - countries_test.dta (country lookup)"
