# Survival

[![Run tests](https://github.com/RaphaelS1/Survival.jl/actions/workflows/run-tests.yml/badge.svg)](https://github.com/RaphaelS1/Survival.jl/actions/workflows/run-tests.yml)
![Experimental](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)

Just used to develop my own Julia skills and to experiment with implementing survival analysis. Builds on existing experience with [mlr3proba](https://github.com/mlr-org/mlr3proba) and [distr6](https://github.com/alan-turing-institute/distr6).

**Don't confuse this package with [JuliaStats/Survival.jl](https://github.com/JuliaStats/Survival.jl).** (I may change the name of this one) - Features here are experimental and sub-optimally coded.

Current features:

* [x] Kaplan Meier estimator (+plotting)
* [x] Nelson Aalen estimator (+plotting)
* [x] Parametric PH models (Exponential and Weibull)
* [x] Parametric AFT models (Exponential and Weibull)
* [x] Discrete and continuous survival prediction objects, including distribution, linear predictor and general risk return types
* [x] Plotting for non parametric estimators
* [x] Surv object for unified censoring indicator (functionality for left, right, interval)
* [x] Extended formula interface for survival objects

Planned features:

* [ ] Cox PH
* [ ] Log-rank tests
* [ ] Residuals (Schoenfeld etc.)
* [ ] Predict type transformations
* [ ] Generic plotting functionality
* [ ] Proportional odds functionality
* [ ] Gompertz, Tobit, Gaussian, Rayleigh, Loglogistic, Lognormal
* [ ] Analytical optimisation
