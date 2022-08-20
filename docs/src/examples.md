# Examples

```@meta
DocTestSetup = quote
    using DataFrames, Distributions, SurvivalAnalysis, Random
end
```

## Fitting

Models can be fit in one of four ways but we only recommend the first.

```jldoctest data
julia> Random.seed!(42);

julia> data = DataFrame(Y = randn(10), D = rand(Bernoulli(), 10))
10×2 DataFrame
 Row │ Y          D
     │ Float64    Bool
─────┼──────────────────
   1 │ -0.363357   true
   2 │  0.251737   true
   3 │ -0.314988  false
   4 │ -0.311252   true
   5 │  0.816307  false
   6 │  0.476738   true
   7 │ -0.859555  false
   8 │ -1.46929    true
   9 │ -0.206613  false
  10 │ -0.310744   true
```

### 1. Function + Formula

```jldoctest data
julia> f = kaplan_meier(@formula(Srv(Y, D) ~ 1), data)
StatsModels.TableStatisticalModel{KaplanMeier, Matrix{Float64}}

(Y,D;+) ~ 1

Coefficients:
KaplanMeier(n = 10, ncens = 4, nevents = 6)
```

### 2. `fit` + Formula

```jldoctest data
julia>  f = fit(KaplanMeier, @formula(Srv(Y, D) ~ 1), data)
StatsModels.TableStatisticalModel{KaplanMeier, Matrix{Float64}}

(Y,D;+) ~ 1

Coefficients:
KaplanMeier(n = 10, ncens = 4, nevents = 6)
```

### 3. Function + Data

```jldoctest data
julia> f = kaplan_meier(hcat(ones(10), 1:10), Surv(data.Y, data.D, "right"))
KaplanMeier(n = 10, ncens = 4, nevents = 6)
```

### 4. `fit` + Data

```jldoctest data
julia> f = fit(KaplanMeier, hcat(ones(10), 1:10), Surv(data.Y, data.D, "right"))
KaplanMeier(n = 10, ncens = 4, nevents = 6)
```

## Predicting

If fitting method (1) or (2) are selected then new data must be given as a DataFrame, otherwise a Matrix is sufficient. We strongly recommend the formula method as this ensures the same covariates and predictors are used in fitting and predicting.


```jldoctest data
julia> Random.seed!(24);

julia> f = kaplan_meier(@formula(Srv(Y, D) ~ 1), data);

julia> predict(f, DataFrame(Y = randn(10), D = trues(10)));
```
