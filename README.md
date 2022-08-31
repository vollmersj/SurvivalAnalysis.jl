# SurvivalAnalysis.jl: A survival analysis interface for Julia

| Status | Testing | Coverage | Docs | Style | Lifecycle | Licence |
| -------| ------- | -------- | ---- | ----- | --------- | ------- |
| [![deps](https://juliahub.com/docs/SurvivalAnalysis/deps.svg)](https://juliahub.com/ui/Packages/SurvivalAnalysis/N9zkY?t=2) <br> [![version](https://juliahub.com/docs/SurvivalAnalysis/version.svg)](https://juliahub.com/ui/Packages/SurvivalAnalysis/N9zkY) <br> [![pkgeval](https://juliahub.com/docs/SurvivalAnalysis/pkgeval.svg)](https://juliahub.com/ui/Packages/SurvivalAnalysis/N9zkY) | [![CI](https://github.com/RaphaelS1/SurvivalAnalysis.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/RaphaelS1/SurvivalAnalysis.jl/actions/workflows/CI.yml) | [![codecov](https://codecov.io/gh/RaphaelS1/SurvivalAnalysis.jl/branch/main/graph/badge.svg?token=R1QK5X4RVP)](https://codecov.io/gh/RaphaelS1/SurvivalAnalysis.jl) | [![](https://img.shields.io/badge/docs-stable-darkblue.svg)](https://raphaels1.github.io/SurvivalAnalysis.jl/stable) <br> [![](https://img.shields.io/badge/docs-dev-lightblue.svg)](https://raphaels1.github.io/SurvivalAnalysis.jl/dev) | [![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle) | ![Experimental](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg) | [![MIT](https://img.shields.io/badge/License-MIT-yelllow)](https://opensource.org/licenses/MIT) |

## About

Survival analysis interface in Julia, still very experimental. Tries to build on my R experience with [mlr3proba](https://github.com/mlr-org/mlr3proba) and [distr6](https://github.com/alan-turing-institute/distr6).

## Related Packages

* [JuliaStats/Survival.jl](https://github.com/JuliaStats/Survival.jl)

## Features

### Current

* [x] Kaplan Meier estimator (+plotting)
* [x] Nelson Aalen estimator (+plotting)
* [x] Parametric PH models (Exponential and Weibull)
* [x] Parametric AFT models (Exponential and Weibull)
* [x] Discrete and continuous survival prediction objects, including distribution, linear predictor and general risk return types
* [x] Plotting for non parametric estimators
* [x] Surv object for unified censoring indicator (functionality for left, right, interval)
* [x] Extended formula interface for survival objects

### Planned

* [ ] [CoxPH](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/8)
* [ ] [Log-rank tests](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/10)
* [ ] [Residuals (Schoenfeld etc.)](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/11)
* [ ] [Predict type transformations](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/12)
* [ ] [Generic plotting functionality](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/13)
* [ ] [Proportional odds functionality](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/14)
* [ ] [More parametric AFTS](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/15)
* [ ] [Analytical optimisation](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/16)
* [ ] [Discrimination measures](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/17)
* [ ] [Scoring rules](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/19)
* [ ] [Calibration measures](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/18)
* [ ] [Bayesian interface](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/20)
* [ ] [mlr3proba integration](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/21)
* [ ] [MLJ integration](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/22)
* [ ] [Documentation](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/9)
