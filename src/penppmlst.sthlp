{smcl}
{* *! version 0.5.0  03jan2026}{...}
{viewerjumpto "Syntax" "penppmlst##syntax"}{...}
{viewerjumpto "Description" "penppmlst##description"}{...}
{viewerjumpto "Options" "penppmlst##options"}{...}
{viewerjumpto "Selection Methods" "penppmlst##selection"}{...}
{viewerjumpto "Examples" "penppmlst##examples"}{...}
{viewerjumpto "Stored results" "penppmlst##results"}{...}
{viewerjumpto "Methods" "penppmlst##methods"}{...}
{viewerjumpto "R vs Stata Differences" "penppmlst##differences"}{...}
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
{synopt:{opt sel:ection(string)}}method: {bf:cv}, {bf:plugin}, {bf:bic}, {bf:aic}, {bf:ebic}, {bf:or {bf:none}{p_end}
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
This command is a Stata implementation of the R package {bf:penppml} by
Breinlich, Corradi, Rocha, Ruta, Santos Silva, and Zylkin (2021).
The primary application is estimating gravity models of international trade
with many potential trade agreement provisions or other covariates.

{pstd}
The command combines:

{p 8 12 2}
1. Iteratively Reweighted Least Squares (IRLS) for Poisson regression

{p 8 12 2}
2. Alternating projections (Gaure 2013) or reghdfe's accelerated methods for high-dimensional fixed effects

{p 8 12 2}
3. Coordinate descent (Friedman, Hastie, Tibshirani 2010) for lasso penalties

{p 8 12 2}
4. Plugin penalty (Belloni et al. 2016) for data-driven regularization

{p 8 12 2}
5. Bootstrap lasso for robust variable selection across resampled datasets

{p 8 12 2}
6. Iceberg lasso for multi-outcome estimation without fixed effects


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


{dlgtab:Lambda Selection}

{phang}
{opt selection(string)} specifies the method for selecting the penalty parameter:

{p 12 16 2}
{bf:none} - use the value specified in {opt lambda()}

{p 12 16 2}
{bf:cv} - K-fold cross-validation; minimizes out-of-sample deviance

{p 12 16 2}
{bf:plugin} - data-driven penalty following Belloni et al. (2016)

{p 12 16 2}
{bf:bic} - Bayesian Information Criterion

{p 12 16 2}
{bf:aic} - Akaike Information Criterion

{p 12 16 2}
{bf:ebic} - Extended Bayesian Information Criterion

{p 12 16 2}
{bf:bootstrap} - bootstrap lasso; runs plugin lasso on B resampled datasets and

{p 12 16 2}
{bf:iceberg} - multi-outcome plugin lasso without fixed effects; useful for
second-stage estimation in the iceberg procedure

{phang}
{opt cluster(varname)} specifies the cluster variable for computing
cluster-robust penalty weights in plugin and bootstrap methods. When specified,
bootstrap resampling is performed at the cluster level.

{phang}
{opt selection(bootstrap)}. Default is 250. More repetitions provide more
stable variable selection but increase computation time.

{phang}
in which a variable must be selected to be included in the final model.
Default is 0.01 (1%). Higher thresholds produce sparser models.


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
lambda scaling, mu bounds, and deviance computation. Use this option when you
need to replicate R results or compare estimates across platforms.


{marker selection}{...}
{title:Selection Methods}

{dlgtab:Plugin Lasso}

{pstd}
The plugin lasso method (Belloni et al. 2016) computes penalty weights based on
the heteroskedasticity structure of the data. The penalty parameter is:

{p 8 8 2}
lambda = c * sqrt(n) * Phi^-1(1 - gamma/(2k))

{pstd}
where c = 1.1, gamma = 0.1/ln(n), and k is the number of regressors. The penalty
loadings psi_j are computed from the score variance:

{p 8 8 2}
psi_j = sqrt( (1/n) * sum_i (x_ij * e_i)^2 )

{pstd}
where e_i are residuals from an initial unpenalized fit.

{dlgtab:Bootstrap Lasso}

{pstd}
Bootstrap lasso provides robust variable selection by:

{p 8 12 2}
1. Drawing B cluster bootstrap samples with replacement

{p 8 12 2}
2. Running plugin lasso on each bootstrap sample

{p 8 12 2}

{pstd}
This procedure is more robust to individual outliers and data perturbations
than single-sample lasso. See Breinlich et al. (2021) for details.

{dlgtab:Iceberg Lasso}

{pstd}
Iceberg lasso is a simplified plugin lasso without fixed effects, designed for
multi-outcome estimation. It is useful in the second stage of the "iceberg"
procedure where first-stage selected variables become dependent variables
in separate regressions. The method uses glmnet-style penalized weighted
least squares within an IRLS loop.


{marker examples}{...}
{title:Examples}

{pstd}Basic penalized PPML with lasso and fixed lambda{p_end}
{phang2}{cmd:. penppmlst trade tariff distance, absorb(i.exporter i.importer i.year) lambda(0.1)}{p_end}

{pstd}Cross-validation for lambda selection with post-lasso{p_end}
{phang2}{cmd:. penppmlst trade tariff distance contig comlang colony, absorb(i.exp#i.year i.imp#i.year) selection(cv) nfolds(5) post}{p_end}

{pstd}Plugin lasso with cluster-robust penalties{p_end}
{phang2}{cmd:. penppmlst trade provision1-provision100, absorb(i.pair i.year) selection(plugin) cluster(pair) post}{p_end}

{pstd}Bootstrap lasso for robust variable selection{p_end}

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
{synopt:{cmd:e(alpha)}}elastic net mixing parameter{p_end}
{synopt:{cmd:e(n_selected)}}number of selected variables{p_end}
{synopt:{cmd:e(ll)}}log pseudo-likelihood{p_end}
{synopt:{cmd:e(deviance)}}final deviance{p_end}
{synopt:{cmd:e(converged)}}1 if converged, 0 otherwise{p_end}
{synopt:{cmd:e(iterations)}}number of IRLS iterations{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:penppmlst}{p_end}
{synopt:{cmd:e(depvar)}}dependent variable name{p_end}
{synopt:{cmd:e(indepvars)}}independent variable names{p_end}
{synopt:{cmd:e(absorb)}}fixed effects specification{p_end}
{synopt:{cmd:e(selected)}}names of selected variables{p_end}
{synopt:{cmd:e(penalty)}}penalty type used{p_end}
{synopt:{cmd:e(selection)}}lambda selection method{p_end}
{synopt:{cmd:e(hdfe)}}HDFE method used (mata or ppmlhdfe){p_end}
{synopt:{cmd:e(r_compatible)}}yes if R-compatible mode{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix (if {opt post} specified){p_end}


{marker methods}{...}
{title:Methods and Formulas}

{pstd}
{cmd:penppmlst} estimates the penalized Poisson pseudo-likelihood:

{p 8 8 2}
max_beta  sum_i [ y_i * (X_i * beta) - exp(X_i * beta + alpha_FE) ] - lambda * P(beta)

{pstd}
where P(beta) is the penalty function:

{p 8 12 2}
{bf:Lasso}: P(beta) = sum_j psi_j |beta_j|

{p 8 12 2}
{bf:Ridge}: P(beta) = (1/2) sum_j psi_j beta_j^2

{p 8 12 2}
{bf:Elastic net}: P(beta) = alpha * sum_j psi_j |beta_j| + (1-alpha)/2 * sum_j psi_j beta_j^2

{pstd}
The algorithm uses an outer IRLS loop with an inner coordinate descent loop:

{p 8 12 2}
1. Initialize mu = 0.5 * (y + mean(y))

{p 8 12 2}
2. Compute IRLS weights: w = mu

{p 8 12 2}
3. Compute working variable: z = log(mu) + (y - mu)/mu

{p 8 12 2}
4. Partial out fixed effects from (z, X) using weighted projections

{p 8 12 2}
5. Solve penalized weighted least squares via coordinate descent

{p 8 12 2}
6. Update eta = X*beta + FE, mu = exp(eta)

{p 8 12 2}
7. Check convergence on deviance; repeat from step 2


{marker differences}{...}
{title:R vs Stata Differences}

{pstd}
{cmd:penppmlst} is designed to produce results comparable to the R {bf:penppml} package.
However, some differences may occur due to implementation details:

{dlgtab:Sources of Numerical Differences}

{p 4 8 2}
{bf:1. HDFE Algorithm}

{p 8 8 2}
R uses {cmd:collapse::fhdwithin} for fixed effects absorption, which uses
simple alternating projections. The default Stata implementation ({cmd:hdfe(mata)})
matches this approach. However, {cmd:hdfe(ppmlhdfe)} uses accelerated methods
(symmetric Kaczmarz, Cimmino) which may converge to slightly different solutions
due to numerical precision.

{p 8 8 2}
{it:Expected magnitude}: < 1e-6 relative difference in coefficients when using
{cmd:hdfe(mata)} with {cmd:r_compatible}; up to 1e-4 with {cmd:hdfe(ppmlhdfe)}.

{p 4 8 2}
{bf:2. Mu Bounds (Fitted Values)}

{p 8 8 2}
R penppml bounds mu in [1e-5, 1e10] while Stata's default (following ppmlhdfe)
uses [1e-10, 1e10]. The {cmd:r_compatible} option enforces R's bounds.

{p 8 8 2}
{it:Expected magnitude}: Affects convergence behavior, rarely affects final
coefficients by more than 1e-6.

{p 4 8 2}
{bf:3. Lambda Scaling}

{p 8 8 2}
R glmnet scales lambda by 1/n (lambda_glmnet = lambda/n), while Stata's
coordinate descent uses unscaled lambda. The internal implementation adjusts
for this, but cross-validation may select slightly different lambda values.

{p 8 8 2}
{it:Expected magnitude}: Lambda values may differ, but selected variables
typically match when using the same relative position on the lambda path.

{p 4 8 2}
{bf:4. Deviance Computation}

{p 8 8 2}
R computes Poisson deviance as 2*sum(y*log(y/mu) - (y-mu)), handling y=0
by setting 0*log(0)=0. Stata's default may use slightly different numerical
handling. Use {cmd:r_compatible} for exact matching.

{p 8 8 2}
{it:Expected magnitude}: < 1e-8 relative difference with {cmd:r_compatible}.

{p 4 8 2}
{bf:5. Random Number Generation}

{p 8 8 2}
Cross-validation fold assignments and bootstrap samples use different RNG
implementations. Set random seed in both R and Stata before running, but
expect different fold assignments even with the same seed.

{p 8 8 2}
{it:Expected magnitude}: Different CV-selected lambdas; use {cmd:selection(plugin)}
for deterministic results.

{p 4 8 2}
{bf:6. Coordinate Descent Implementation}

{p 8 8 2}
R uses glmnet's Fortran-optimized coordinate descent. Stata uses pure Mata.
Both implement the same algorithm but may have different numerical precision
in the soft-thresholding and weighted inner product computations.

{p 8 8 2}
{it:Expected magnitude}: < 1e-7 relative difference per coefficient.

{dlgtab:Achieving Maximum Compatibility}

{pstd}
To maximize reproducibility with R penppml:

{phang2}{cmd:. penppmlst ... , r_compatible selection(plugin) hdfe(mata)}{p_end}

{pstd}
This combination:

{p 8 12 2}
- Uses R-compatible mu bounds [1e-5, 1e10]

{p 8 12 2}
- Uses pure alternating projections for HDFE (matching collapse::fhdwithin)

{p 8 12 2}
- Uses R-compatible deviance computation

{p 8 12 2}
- Uses deterministic plugin lambda (no RNG dependence)

{dlgtab:When Differences Matter}

{pstd}
Small numerical differences typically do not affect:

{p 8 12 2}
- Which variables are selected (same sparsity pattern)

{p 8 12 2}
- Signs of coefficients

{p 8 12 2}
- Qualitative conclusions about effect sizes

{pstd}
Differences may be meaningful when:

{p 8 12 2}
- Coefficients are very close to zero (on the selection boundary)

{p 8 12 2}
- The model is near-singular or poorly conditioned

{p 8 12 2}
- Comparing exact numerical values for replication purposes


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

{phang}
Friedman, J., T. Hastie, and R. Tibshirani. 2010.
Regularization paths for generalized linear models via coordinate descent.
{it:Journal of Statistical Software} 33(1): 1-22.

{phang}
Gaure, S. 2013.
OLS with multiple high dimensional category variables.
{it:Computational Statistics & Data Analysis} 66: 8-18.


{marker authors}{...}
{title:Authors}

{pstd}
{bf:Stata implementation:}

{pstd}
Erdey, Laszlo (2026){break}
Faculty of Economics and Business{break}
University of Debrecen{break}
Debrecen, Hungary

{pstd}
{bf:Original R package penppml:}

{pstd}
Diego Ferreras Garrucho{break}
Tom Zylkin{break}
Joao Cruz{break}
Nicolas Apfel

{pstd}
Based on methodology by:{break}
Holger Breinlich, Valentina Corradi, Nadia Rocha, Michele Ruta,{break}
J.M.C. Santos Silva, Tom Zylkin

{pstd}
R package: {browse "https://github.com/tomzylkin/penppml"}


{marker also}{...}
{title:Also see}

{psee}
{space 2}Help:  {help ppmlhdfe}, {help reghdfe}, {help lasso}, {help elasticnet}
{p_end}
