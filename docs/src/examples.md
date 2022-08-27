# Examples

```@meta
DocTestSetup = quote
    using DataFrames, Distributions, SurvivalAnalysis, Random
end
```

# Fitting and predicting

## Fitting

Models can be fit in one of four ways but we only recommend the first.

```jldoctest data
julia> data = DataFrame(Y = [1,1,4,6,8,4,9,4,5,10], D = [true, false, false, false, true, false, false, true, true, false])
10×2 DataFrame
 Row │ Y      D
     │ Int64  Bool
─────┼──────────────
   1 │     1   true
   2 │     1  false
   3 │     4  false
   4 │     6  false
   5 │     8   true
   6 │     4  false
   7 │     9  false
   8 │     4   true
   9 │     5   true
  10 │    10  false
```

### 1. Function + Formula

```jldoctest data
julia> f = kaplan_meier(@formula(Srv(Y, D) ~ 1), data)
StatsModels.TableStatisticalModel{KaplanMeier, Matrix{Float64}}

(Y,D;+) ~ 1

Coefficients:
  n  ncens  nevents
 10      6        4
```

### 2. `fit` + Formula

```jldoctest data
julia>  f = fit(KaplanMeier, @formula(Srv(Y, D) ~ 1), data)
StatsModels.TableStatisticalModel{KaplanMeier, Matrix{Float64}}

(Y,D;+) ~ 1

Coefficients:
  n  ncens  nevents
 10      6        4
```

### 3. Function + Data

```jldoctest data
julia> f = kaplan_meier(hcat(ones(10), 1:10), Surv(data.Y, data.D, :r))
KaplanMeier

Coefficients:
  n  ncens  nevents
 10      6        4
```

### 4. `fit` + Data

```jldoctest data
julia> f = fit(KaplanMeier, hcat(ones(10), 1:10), Surv(data.Y, data.D, :right))
KaplanMeier

Coefficients:
  n  ncens  nevents
 10      6        4
```

## Predicting

If fitting method (1) or (2) are selected then new data must be given as a DataFrame, otherwise a Matrix is sufficient. We strongly recommend the formula method as this ensures the same covariates and predictors are used in fitting and predicting.


```jldoctest data
julia> f = kaplan_meier(@formula(Srv(Y, D) ~ 1), data);

julia> predict(f, DataFrame(Y = randn(10), D = trues(10)));
```
