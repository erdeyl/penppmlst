* Test help file format
clear all
set more off

* Set adopath to include local src folder
adopath + "C:/Users/erdey/Claude/penppmlst/src"

* Try to view help file
capture noisily help penppmlst

* Check for any syntax errors
display "Help file test complete"
