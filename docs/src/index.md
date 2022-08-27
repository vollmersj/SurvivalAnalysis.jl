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

* [Tutorials](tutorials.md)
* [How-to guides](howto.md)
* [Explanations](explanations.md)
* [Reference](api.md)

[![Diataxis](https://diataxis.fr/_images/diataxis.png)](https://diataxis.fr/)

In this package the function references can be found at [api.md](api.md). We
will later be adding detailed tutorials, how-to guides, and explanations - for
now we just have placeholder files as well as the documentation found within
function references, which all aim to include some element of 'explanation'
and 'how-to'.

## Package Features
- Parametric models
- Semi-parametric models
- Nonparametric estimators

## Upcoming features
- Bayesian methods
- Survival metrics (discrimination, calibration, scoring rules)
