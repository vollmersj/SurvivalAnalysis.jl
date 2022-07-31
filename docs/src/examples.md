# Examples

```@meta
DocTestSetup = quote
    using DataFrames, Distributions, SurvivalAnalysis
end
```

## Fitting

Models can be fit in one of four ways but we only recommend the first.

### 1. Function + Formula

```jldoctest
julia> using Random

julia> f = kaplan_meier(@formula(Srv(Y, D) ~ 1), DataFrame(Y = randn(10), D = trues(10)))

julia> typeof(f)
StatsModels.TableStatisticalModel{KaplanMeier, Matrix{Float64}}
```

### 2. `fit` + Formula

```jldoctest
julia> using Random

julia>  f = fit(KaplanMeier, @formula(Srv(Y, D) ~ 1), DataFrame(Y = randn(10), D = trues(10)))

julia> typeof(f)
StatsModels.TableStatisticalModel{KaplanMeier, Matrix{Float64}}
```

### 3. Function + Data

```jldoctest
julia> using Random

julia> f = kaplan_meier(hcat(ones(10), 1:10), Surv(randn(MersenneTwister(42), 10)))

julia> typeof(f)
KaplanMeier
```

### 4. `fit` + Data

```jldoctest
julia> using Random

julia> f = fit(KaplanMeier, hcat(ones(10), 1:10), Surv(randn(MersenneTwister(42), 10)))

julia> typeof(f)
KaplanMeier
```

## Predicting

If fitting method (1) or (2) are selected then new data must be given as a DataFrame, otherwise a Matrix is sufficient. We strongly recommend the formula method as this ensures the same covariates and predictors are used in fitting and predicting.


```jldoctest
julia> using Random

julia> f = kaplan_meier(@formula(Srv(Y, D) ~ 1), DataFrame(Y = randn(10), D = trues(10)))

julia> typeof(f)
StatsModels.TableStatisticalModel{KaplanMeier, Matrix{Float64}}

julia> predict(f, DataFrame(Y = randn(10), D = trues(10)))
```
