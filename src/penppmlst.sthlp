{smcl}
{* *! version 0.3.0  03jan2026}{...}
{viewerjumpto "Syntax" "penppmlst##syntax"}{...}
{viewerjumpto "Description" "penppmlst##description"}{...}
{viewerjumpto "Options" "penppmlst##options"}{...}
{viewerjumpto "Examples" "penppmlst##examples"}{...}
{viewerjumpto "Stored results" "penppmlst##results"}{...}
{viewerjumpto "Methods" "penppmlst##methods"}{...}
{viewerjumpto "References" "penppmlst##references"}{...}
{viewerjumpto "Authors" "penppmlst##authors"}{...}

{title:Title}

{p2colset 5 18 20 2}{...}
{p2col:{bf:penppmlst} {hline 2}}Penalized Poisson Pseudo Maximum Likelihood with High-Dimensional Fixed Effects{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:penppmlst}
{depvar}
{indepvars}
{ifin}
{weight}{cmd:,}
{opt absorb(absvars)}
[{it:options}]

{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{opt absorb(absvars)}}categorical variables identifying fixed effects; required{p_end}
{synopt:{opt pen:alty(string)}}penalty type: {bf:lasso}, {bf:ridge}, or {bf:elasticnet}; default is {bf:lasso}{p_end}
{synopt:{opt lam:bda(#)}}penalty parameter; if not specified, uses selection method{p_end}
{synopt:{opt al:pha(#)}}elastic net mixing parameter; 1=lasso, 0=ridge; default is 1{p_end}

{syntab:Lambda Selection}
{synopt:{opt sel:ection(string)}}method: {bf:cv}, {bf:plugin}, {bf:bic}, {bf:aic}, {bf:ebic}, or {bf:none}{p_end}
{synopt:{opt nf:olds(#)}}number of cross-validation folds; default is 10{p_end}
{synopt:{opt nl:ambda(#)}}number of lambda values for regularization path; default is 100{p_end}
{synopt:{opt lminratio(#)}}ratio of lambda_min to lambda_max; default is 0.01{p_end}
{synopt:{opt cl:uster(varname)}}cluster variable for plugin penalty weights{p_end}

{syntab:Estimation}
{synopt:{opt post}}compute post-lasso estimates (unpenalized on selected variables){p_end}
{synopt:{opt nostand:ardize}}do not standardize regressors before penalization{p_end}
{synopt:{opt tol:erance(#)}}convergence tolerance; default is 1e-8{p_end}
{synopt:{opt max:iter(#)}}maximum IRLS iterations; default is 1000{p_end}

{syntab:HDFE Method}
{synopt:{opt hdfe(string)}}HDFE backend: {bf:mata} or {bf:ppmlhdfe}; default is {bf:mata}{p_end}
{synopt:{opt r_compatible}}use R penppml-compatible settings for reproducibility{p_end}

{syntab:Reporting}
{synopt:{opt irr}}report incidence rate ratios{p_end}
{synopt:{opt eform}}same as {opt irr}{p_end}
{synopt:{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt:{opt verb:ose}}display iteration log{p_end}
{synopt:{opt nolog}}suppress iteration log{p_end}
{synoptline}
{p2colreset}{...}

{p 4 6 2}
{cmd:fweight}s, {cmd:aweight}s, and {cmd:pweight}s are allowed; see {help weight}.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:penppmlst} fits penalized Poisson Pseudo Maximum Likelihood (PPML) regressions
with high-dimensional fixed effects (HDFE). The command supports lasso, ridge,
and elastic net penalties, with automatic penalty parameter selection via
cross-validation, plugin methods, or information criteria.

{pstd}
This command is a Stata implementation of the R package {bf:penppmlst} by
Breinlich, Corradi, Rocha, Ruta, Santos Silva, and Zylkin (2021).
The primary application is estimating gravity models of international trade
with many potential trade agreement provisions or other covariates.

{pstd}
The command combines:

{p 8 12 2}
1. Iteratively Reweighted Least Squares (IRLS) for Poisson regression

{p 8 12 2}
2. Alternating projections (Gaure 2013) for high-dimensional fixed effects

{p 8 12 2}
3. Coordinate descent (Friedman, Hastie, Tibshirani 2010) for lasso penalties

{p 8 12 2}
4. Plugin penalty (Belloni et al. 2016) for data-driven regularization


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opt absorb(absvars)} specifies the categorical variables that identify
the fixed effects to be absorbed. Multiple fixed effects can be specified,
and interactions are supported using the {cmd:#} operator. This option is required.

{pmore}
Example: {cmd:absorb(i.exporter i.importer i.year)}

{pmore}
Example with interactions: {cmd:absorb(i.exporter#i.year i.importer#i.year)}

{phang}
{opt penalty(string)} specifies the type of penalty:

{p 12 16 2}
{bf:lasso} - L1 penalty (absolute value of coefficients); produces sparse solutions

{p 12 16 2}
{bf:ridge} - L2 penalty (squared coefficients); shrinks but does not set to zero

{p 12 16 2}
{bf:elasticnet} - combination of L1 and L2, controlled by {opt alpha()}

{phang}
{opt lambda(#)} specifies the penalty parameter. Larger values result in more
shrinkage. If not specified, lambda is selected using the method in {opt selection()}.

{phang}
{opt alpha(#)} controls the elastic net mixing. {cmd:alpha(1)} is pure lasso,
{cmd:alpha(0)} is pure ridge, and values in between give elastic net.
Default is {cmd:alpha(1)}.


{dlgtab:HDFE Method}

{phang}
{opt hdfe(string)} specifies the backend for absorbing high-dimensional fixed effects:

{p 12 16 2}
{bf:mata} - Pure Mata implementation using alternating projections (Gaure 2013).
This is the default and is fully compatible with the R penppml package, producing
numerically comparable results. Recommended when reproducibility with R is required.

{p 12 16 2}
{bf:ppmlhdfe} - Uses reghdfe's optimized FixedEffects Mata class (the same
infrastructure used by ppmlhdfe). Significantly faster for large datasets with
many fixed effects due to acceleration methods (Cimmino, symmetric Kaczmarz).
Requires {cmd:ssc install ftools} and {cmd:ssc install reghdfe}.
Results may differ slightly from R due to algorithmic differences.

{phang}
{opt r_compatible} forces settings that produce numerically comparable results
to the R penppml package. Implies {cmd:hdfe(mata)} and uses R-compatible
lambda scaling. Use this option when you need to replicate R results or
compare estimates across platforms.


{marker examples}{...}
{title:Examples}

{pstd}Basic penalized PPML with lasso and fixed lambda{p_end}
{phang2}{cmd:. penppmlst trade tariff distance, absorb(i.exporter i.importer i.year) lambda(0.1)}{p_end}

{pstd}Cross-validation for lambda selection with post-lasso{p_end}
{phang2}{cmd:. penppmlst trade tariff distance contig comlang colony, absorb(i.exp#i.year i.imp#i.year) selection(cv) nfolds(5) post}{p_end}

{pstd}Plugin lasso with cluster-robust penalties{p_end}
{phang2}{cmd:. penppmlst trade provision1-provision100, absorb(i.pair i.year) selection(plugin) cluster(pair) post}{p_end}

{pstd}Ridge regression{p_end}
{phang2}{cmd:. penppmlst trade gravity_vars, absorb(i.exp#i.imp) penalty(ridge) lambda(1)}{p_end}

{pstd}Using ppmlhdfe backend for faster HDFE computation{p_end}
{phang2}{cmd:. penppmlst trade tariff distance, absorb(i.exp i.imp i.year) hdfe(ppmlhdfe) selection(cv)}{p_end}

{pstd}R-compatible estimation for cross-platform reproducibility{p_end}
{phang2}{cmd:. penppmlst trade provision1-provision100, absorb(i.pair i.year) selection(plugin) r_compatible}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:penppmlst} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(lambda)}}penalty parameter used{p_end}
{synopt:{cmd:e(n_selected)}}number of selected variables{p_end}
{synopt:{cmd:e(ll)}}log pseudo-likelihood{p_end}
{synopt:{cmd:e(deviance)}}final deviance{p_end}
{synopt:{cmd:e(converged)}}1 if converged, 0 otherwise{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:penppmlst}{p_end}
{synopt:{cmd:e(depvar)}}dependent variable name{p_end}
{synopt:{cmd:e(selected)}}names of selected variables{p_end}
{synopt:{cmd:e(penalty)}}penalty type used{p_end}
{synopt:{cmd:e(hdfe)}}HDFE method used (mata or ppmlhdfe){p_end}
{synopt:{cmd:e(r_compatible)}}1 if R-compatible mode, 0 otherwise{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix (if {opt post} specified){p_end}


{marker references}{...}
{title:References}

{phang}
Belloni, A., V. Chernozhukov, C. Hansen, and D. Kozbur. 2016.
Inference in high-dimensional panel models with an application to gun control.
{it:Journal of Business & Economic Statistics} 34: 590-605.

{phang}
Breinlich, H., V. Corradi, N. Rocha, M. Ruta, J.M.C. Santos Silva, and T. Zylkin. 2021.
Machine learning in international trade research: Evaluating the impact of trade agreements.
{it:Policy Research Working Paper} 9629. World Bank.

{phang}
Correia, S., P. Guimaraes, and T. Zylkin. 2020.
Fast Poisson estimation with high-dimensional fixed effects.
{it:Stata Journal} 20: 95-115.


{marker authors}{...}
{title:Authors}

{pstd}
{bf:Stata implementation:}

{pstd}
Erdey, László (2026){break}
Faculty of Economics and Business{break}
University of Debrecen{break}
Debrecen, Hungary

{pstd}
{bf:Original R package penppmlst:}

{pstd}
Diego Ferreras Garrucho{break}
Tom Zylkin{break}
Joao Cruz{break}
Nicolas Apfel

{pstd}
R package: {browse "https://github.com/tomzylkin/penppmlst"}


{marker also}{...}
{title:Also see}

{psee}
{space 2}Help:  {help ppmlhdfe}, {help reghdfe}, {help lasso}
{p_end}
