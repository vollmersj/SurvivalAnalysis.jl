# SurvivalAnalysis

[![CI](https://github.com/RaphaelS1/SurvivalAnalysis.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/RaphaelS1/SurvivalAnalysis.jl/actions/workflows/CI.yml)
![Experimental](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)
[![codecov](https://codecov.io/gh/RaphaelS1/SurvivalAnalysis.jl/branch/main/graph/badge.svg?token=R1QK5X4RVP)](https://codecov.io/gh/RaphaelS1/SurvivalAnalysis.jl)
[![MIT](https://img.shields.io/badge/License-MIT-yelllow)](https://opensource.org/licenses/MIT)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

Just used to develop my own Julia skills and to experiment with implementing survival analysis. Builds on existing experience with [mlr3proba](https://github.com/mlr-org/mlr3proba) and [distr6](https://github.com/alan-turing-institute/distr6).

**Don't confuse this package with [JuliaStats/Survival.jl](https://github.com/JuliaStats/Survival.jl).** - SurvivalAnalysis.jl is experimental and sub-optimally coded.

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
* [ ] Discrimination measures
* [ ] Scoring rules
* [ ] Calibration measures
* [ ] Bayesian interface
* [ ] Documentation
