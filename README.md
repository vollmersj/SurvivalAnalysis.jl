# SurvivalAnalysis.jl: A survival analysis interface for Julia
<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-1-orange.svg?style=flat-square)](#contributors-)
<!-- ALL-CONTRIBUTORS-BADGE:END -->

| Status | Testing | Coverage | Docs | Codebase | Contributions |
| -------| ------- | -------- | ---- | ----- | ------------- |
| [![deps](https://juliahub.com/docs/SurvivalAnalysis/deps.svg)](https://juliahub.com/ui/Packages/SurvivalAnalysis/N9zkY?t=2) <br> [![version](https://juliahub.com/docs/SurvivalAnalysis/version.svg)](https://juliahub.com/ui/Packages/SurvivalAnalysis/N9zkY) <br> [![pkgeval](https://juliahub.com/docs/SurvivalAnalysis/pkgeval.svg)](https://juliahub.com/ui/Packages/SurvivalAnalysis/N9zkY) | [![CI](https://github.com/RaphaelS1/SurvivalAnalysis.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/RaphaelS1/SurvivalAnalysis.jl/actions/workflows/CI.yml) | [![codecov](https://codecov.io/gh/RaphaelS1/SurvivalAnalysis.jl/branch/main/graph/badge.svg?token=R1QK5X4RVP)](https://codecov.io/gh/RaphaelS1/SurvivalAnalysis.jl) | [![](https://img.shields.io/badge/docs-stable-darkblue.svg)](https://raphaels1.github.io/SurvivalAnalysis.jl/stable) <br> [![](https://img.shields.io/badge/docs-dev-lightblue.svg)](https://raphaels1.github.io/SurvivalAnalysis.jl/dev) | [![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/RaphaelS1/SurvivalAnalysis.jl/blob/main/code_style_blue.md) <br> ![Experimental](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg) <br> [![MIT](https://img.shields.io/badge/License-MIT-yelllow)](https://opensource.org/licenses/MIT) | <!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->[![All Contributors](https://img.shields.io/badge/all_contributors-13-orange.svg?style=flat-square)](#contributors)<!-- ALL-CONTRIBUTORS-BADGE:END --> |

## About

Survival analysis interface in Julia, still very experimental. Tries to build on my R experience with [mlr3proba](https://github.com/mlr-org/mlr3proba) and [distr6](https://github.com/alan-turing-institute/distr6).

## Related Packages

* [JuliaStats/Survival.jl](https://github.com/JuliaStats/Survival.jl) implements right-censored time-to-event types, Cox PH, Kaplan-Meier and Nelson-Aalen estimators. Last commit 2022.
* [Testispuncher/Survival.jl](https://github.com/Testispuncher/Survival.jl) implements non-parametric estimators and some plots. Last commit 2017.
* [kkholst/EventHistory.jl](https://github.com/kkholst/EventHistory.jl) implements different censoring and truncation types and Cox PH. Last commit 2018.

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
* [x] [Log-rank tests](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/10)

### Planned

* [ ] [CoxPH](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/8)
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

## Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center"><a href="http://www.raphaelsonabend.co.uk"><img src="https://avatars.githubusercontent.com/u/25639974?v=4?s=100" width="100px;" alt="Raphael Sonabend"/><br /><sub><b>Raphael Sonabend</b></sub></a><br /><a href="https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues?q=author%3ARaphaelS1" title="Bug reports">🐛</a> <a href="https://github.com/RaphaelS1/SurvivalAnalysis.jl/commits?author=RaphaelS1" title="Code">💻</a> <a href="#content-RaphaelS1" title="Content">🖋</a> <a href="https://github.com/RaphaelS1/SurvivalAnalysis.jl/commits?author=RaphaelS1" title="Documentation">📖</a> <a href="#design-RaphaelS1" title="Design">🎨</a> <a href="#example-RaphaelS1" title="Examples">💡</a> <a href="#ideas-RaphaelS1" title="Ideas, Planning, & Feedback">🤔</a> <a href="#maintenance-RaphaelS1" title="Maintenance">🚧</a> <a href="#projectManagement-RaphaelS1" title="Project Management">📆</a> <a href="#question-RaphaelS1" title="Answering Questions">💬</a> <a href="#research-RaphaelS1" title="Research">🔬</a> <a href="https://github.com/RaphaelS1/SurvivalAnalysis.jl/pulls?q=is%3Apr+reviewed-by%3ARaphaelS1" title="Reviewed Pull Requests">👀</a> <a href="https://github.com/RaphaelS1/SurvivalAnalysis.jl/commits?author=RaphaelS1" title="Tests">⚠️</a> <a href="#tutorial-RaphaelS1" title="Tutorials">✅</a></td>
    </tr>
  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
