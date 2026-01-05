# penppmlst

Penalized Poisson Pseudo Maximum Likelihood with High-Dimensional Fixed Effects for Stata

## Description

`penppmlst` is a Stata package for estimating penalized PPML regressions with lasso, ridge, or elastic net penalties for models with high-dimensional fixed effects.

This is a Stata implementation based on the R package [penppml](https://github.com/tomzylkin/penppml) by Breinlich, Corradi, Rocha, Ruta, Santos Silva, and Zylkin (2021).

## Features

- **Lasso, Ridge, and Elastic Net penalties** for variable selection and regularization
- **High-Dimensional Fixed Effects** via alternating projections (Gaure 2013) or ppmlhdfe backend
- **Cross-Validation** for automatic penalty parameter selection with FE-aware scoring
- **Plugin Lasso** with heteroskedasticity-robust and cluster-robust penalty weights (Belloni et al. 2016)
- **Post-lasso estimation** with valid standard errors
- **Information Criteria** selection (AIC, BIC, EBIC)
- **R-compatible mode** for exact replication of R penppml results
- **Comprehensive predict** with fitted values, residuals, and FE contributions

## Installation

### From GitHub

```stata
net install penppmlst, from("https://raw.githubusercontent.com/erdeyl/penppmlst/main/src") replace
```

### Dependencies

Requires Stata 17 or higher, plus:

```stata
ssc install ftools
```

Optional (for `hdfe(ppmlhdfe)` backend):

```stata
ssc install reghdfe
ssc install ppmlhdfe
```

## Syntax

```stata
penppmlst depvar indepvars, absorb(absvars) [options]
```

### Required Options

| Option | Description |
|--------|-------------|
| `absorb(absvars)` | Fixed effects specification (same syntax as reghdfe/ppmlhdfe) |

### Model Options

| Option | Description |
|--------|-------------|
| `penalty(lasso\|ridge\|elasticnet)` | Penalty type (default: lasso) |
| `lambda(#)` | Penalty parameter |
| `alpha(#)` | Elastic net mixing: 1=lasso, 0=ridge (default: 1) |

### Lambda Selection

| Option | Description |
|--------|-------------|
| `selection(cv\|plugin\|bic\|aic\|ebic\|none)` | Selection method |
| `nfolds(#)` | Number of CV folds (default: 10) |
| `nlambda(#)` | Number of lambda values (default: 100) |
| `cluster(varname)` | Cluster variable for plugin weights |

### Estimation Options

| Option | Description |
|--------|-------------|
| `post` | Compute post-lasso estimates with standard errors |
| `nostandardize` | Do not standardize regressors |
| `tolerance(#)` | Convergence tolerance (default: 1e-8) |
| `maxiter(#)` | Maximum iterations (default: 1000) |
| `hdfe(mata\|ppmlhdfe)` | HDFE backend (default: mata) |
| `r_compatible` | Use R penppml-compatible settings |
| `d(varname)` | Store FE contribution for predict |

### Reporting

| Option | Description |
|--------|-------------|
| `irr` / `eform` | Report incidence rate ratios |
| `level(#)` | Confidence level |
| `verbose` | Show iteration log |

## Postestimation

### Predict

After estimation, use `predict` to generate fitted values and residuals:

```stata
predict newvar [if] [in] , [mu xb xbd d residuals deviance pearson anscombe scores]
```

| Option | Description |
|--------|-------------|
| `mu` | Predicted mean exp(xb+d) (default) |
| `xb` | Linear predictor (without FE) |
| `xbd` | Linear predictor including FE |
| `d` | Fixed effect contribution only |
| `residuals` | Response residuals (y - mu) |
| `deviance` | Deviance residuals |
| `pearson` | Pearson residuals |
| `anscombe` | Anscombe residuals |

**Note:** Predictions requiring FE (mu, xbd, d, residuals) need the `d()` option in the estimation command.

## Examples

```stata
* Basic lasso with fixed lambda
penppmlst trade tariff distance, absorb(i.exporter i.importer i.year) ///
    penalty(lasso) lambda(0.1)

* Cross-validation for lambda selection with post-lasso
penppmlst trade x1-x100, absorb(i.pair i.year) ///
    selection(cv) nfolds(5) post

* Plugin lasso with cluster-robust penalties
penppmlst trade provisions, absorb(i.exp#i.year i.imp#i.year) ///
    selection(plugin) cluster(pair) post

* Ridge regression
penppmlst trade gravity_vars, absorb(i.exp#i.imp) ///
    penalty(ridge) lambda(1)

* Elastic net
penppmlst trade vars, absorb(i.exp i.imp i.year) ///
    penalty(elasticnet) alpha(0.5) selection(cv)

* With FE contribution for predictions
penppmlst trade tariff distance, absorb(i.exp i.imp i.year) ///
    selection(plugin) d(fe_contrib)
predict fitted_values, mu
predict residuals, residuals

* R-compatible mode for cross-platform reproducibility
penppmlst trade provisions, absorb(i.pair i.year) ///
    selection(plugin) r_compatible d(fe_contrib)
```

## Stored Results

### Scalars

| Result | Description |
|--------|-------------|
| `e(N)` | Number of observations |
| `e(lambda)` | Penalty parameter used |
| `e(n_selected)` | Number of selected variables |
| `e(ll)` | Log pseudo-likelihood |
| `e(deviance)` | Final deviance |
| `e(converged)` | Convergence indicator |

### Macros

| Result | Description |
|--------|-------------|
| `e(selected)` | Names of selected variables |
| `e(penalty)` | Penalty type used |
| `e(selection)` | Selection method |
| `e(hdfe)` | HDFE method used |
| `e(r_compatible)` | "yes" if R-compatible mode |
| `e(d)` | FE contribution variable name |

### Matrices

| Result | Description |
|--------|-------------|
| `e(b)` | Coefficient vector |
| `e(V)` | Variance-covariance matrix (if `post` specified) |

## Methodology

The algorithm combines:

1. **IRLS** (Iteratively Reweighted Least Squares) for Poisson regression
2. **Alternating projections** (Gaure 2013) for high-dimensional fixed effects
3. **Coordinate descent** (Friedman, Hastie, Tibshirani 2010) for lasso penalties
4. **Plugin penalty** (Belloni et al. 2016) for data-driven regularization

## References

Belloni, A., V. Chernozhukov, C. Hansen, and D. Kozbur. 2016. "Inference in high-dimensional panel models with an application to gun control." *Journal of Business & Economic Statistics* 34: 590-605.

Breinlich, H., V. Corradi, N. Rocha, M. Ruta, J.M.C. Santos Silva, and T. Zylkin. 2021. "Machine learning in international trade research: Evaluating the impact of trade agreements." *Policy Research Working Paper* 9629. World Bank.

Correia, S., P. Guimaraes, and T. Zylkin. 2020. "Fast Poisson estimation with high-dimensional fixed effects." *Stata Journal* 20: 95-115.

Friedman, J., T. Hastie, and R. Tibshirani. 2010. "Regularization paths for generalized linear models via coordinate descent." *Journal of Statistical Software* 33: 1-22.

## Author

Erdey, László  
Faculty of Economics and Business  
University of Debrecen, Hungary

## License

MIT License
