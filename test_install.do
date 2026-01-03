* Test installation of penppmlst from GitHub
* Run this in Stata to verify the package installs correctly

clear all
set more off

* Display Stata version
di "Stata version: " c(stata_version)
di "Current date: " c(current_date)

* ============================================================================
* Test 1: net install from GitHub
* ============================================================================
di _n(2) "=============================================="
di "TEST 1: net install from GitHub"
di "=============================================="

* First, remove any existing installation
cap ado uninstall penppmlst

* Try to install from GitHub raw content
net from "https://raw.githubusercontent.com/erdeyl/penppmlst/main/src"
net describe penppmlst
net install penppmlst, replace

* Check if files were installed
which penppmlst

* ============================================================================
* Test 2: Check help file
* ============================================================================
di _n(2) "=============================================="
di "TEST 2: Check help file loads"
di "=============================================="

help penppmlst

* ============================================================================
* Test 3: Basic syntax check (no data, should give informative error)
* ============================================================================
di _n(2) "=============================================="
di "TEST 3: Basic syntax check"
di "=============================================="

cap noi penppmlst
if _rc {
    di "Expected error (no variables specified): " _rc
}

* ============================================================================
* Summary
* ============================================================================
di _n(2) "=============================================="
di "INSTALLATION TEST COMPLETE"
di "=============================================="
di "If you see this message, the package installed successfully."
di "Run 'help penppmlst' to see documentation."
