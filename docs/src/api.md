# API

```@meta
DocTestSetup = quote
    using DataFrames, Distributions, SurvivalAnalysis
end
```


```@docs
predict(::SurvivalEstimator, ::DataFrame)
predict(::ParametricPH, ::AbstractMatrix{<:Real})
predict(::ParametricAFT, ::AbstractMatrix{<:Real})
SurvivalPrediction
fit(::Type{<:SurvivalEstimator}, ::AbstractMatrix{<:Real}, ::SurvivalAnalysis.RCSurv}
fit(::Type{<:ParametricSurvival}, ::AbstractMatrix{<:Real}, ::SurvivalAnalysis.RCSurv,
    d::Type{<:ContinuousUnivariateDistribution}; ::Number = 1)
```

## Measures

```@docs
ConcordanceWeights
concordance
```

## Surv

```@docs
Surv
Surv(start::Union{Vector{T}, T} where T <: Number, stop::Union{Vector{T}, T} where T <: Number)
total_censored
total_risk
total_outcomes
total_events
threshold_risk
length
outcome_times
event_times
merge
outcome_status
unique_outcome_times
unique_event_times
survival(::SurvivalEstimator)
time
Srv
reverse(::SurvivalAnalysis.OneSidedSurv)
```

## SurvivalEstimator

```@docs
SurvivalEstimator
NelsonAalen
KaplanMeier
distr
nelson_aalen
kaplan_meier
std
surv_stats
confint(::KaplanMeier, ::Number)
confint(::NelsonAalen, ::Number)
```

## SurvivalModel

```@docs
SurvivalModel
```

## ParametricSurvival

```@docs
ParametricSurvival
ParametricAFT
ParametricPH
ph
aft
scale
baseline
```

## Survival distribution functions

```@docs
cum_hazard
survival
hazard
```

## Plotting

```@autodocs
Modules = [SurvivalAnalysis]
Pages   = ["plots.jl"]
```

# Index

```@index
```
