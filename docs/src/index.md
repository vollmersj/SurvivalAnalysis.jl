# SurvivalAnalysis.jl Manual

Survival analysis in Julia.

## Installation

Install the latest version from JuliaRegistries with

```julia
Pkg.add("SurvivalAnalysis")
```

Or install the latest stable version from GitHub with

```julia
Pkg.add(url="https://github.com/RaphaelS1/SurvivalAnalysis.jl.git")
```

## Using this manual - Diátaxis

Throughout this package we have tried to adhere to the [Diátaxis](https://diataxis.fr/) documentation style. This means clearly separating:

* [Tutorials](tutorials.md) - Learning-oriented lessons guiding reader through steps required to complete a project.
* [How-to guides](howto.md) - Goal-oriented directions taking reader through steps required to solve a real-world problem.
* [Explanations](explanations.md) - Understanding-oriented discussion to clarify and illuminate a particular topic.
* [Reference (API)](api.md) - Information-oriented technical descriptions of the software and how to operate it.

[![Diataxis](https://diataxis.fr/_images/diataxis.png)](https://diataxis.fr/)

As well as having specific pages for these areas, we will also aim to include some element of 'explanation' and 'how-to' into the function and and type reference documentation were sensible.

## Package features

* Kaplan Meier estimator (+plotting)
* Nelson Aalen estimator (+plotting)
* Parametric PH models (Exponential and Weibull)
* Parametric AFT models (Exponential and Weibull)
* Discrete and continuous survival prediction objects, including distribution, linear predictor and general risk return types
* Plotting for non parametric estimators
* Surv object for unified censoring indicator (functionality for left, right, interval)
* Extended formula interface for survival objects

### Upcoming features

* [CoxPH](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/8)
* [Log-rank tests](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/10)
* [Residuals (Schoenfeld etc.)](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/11)
* [Predict type transformations](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/12)
* [Generic plotting functionality](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/13)
* [Proportional odds functionality](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/14)
* [More parametric AFTS](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/15)
* [Analytical optimisation](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/16)
* [Discrimination measures](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/17)
* [Scoring rules](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/19)
* [Calibration measures](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/18)
* [Bayesian interface](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/20)
* [mlr3proba integration](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/21)
* [MLJ integration](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/22)
* [Documentation](https://github.com/RaphaelS1/SurvivalAnalysis.jl/issues/9)
