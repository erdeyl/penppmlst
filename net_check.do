* Net install check for penppmlst
clear all
set more off

display "=============================================="
display "Testing net install from GitHub"
display "=============================================="

* Check net describe from GitHub raw
capture noisily net from "https://raw.githubusercontent.com/erdeyl/penppmlst/main/src"
if _rc == 0 {
    display "net from successful"
    capture noisily net describe penppmlst
}
else {
    display "net from failed with error: " _rc
}

display ""
display "=============================================="
display "SSC check (package not on SSC yet)"
display "=============================================="

* SSC would require submission - just check if command works
capture noisily ssc describe penppmlst
if _rc != 0 {
    display "penppmlst not found on SSC (expected - not yet submitted)"
}
