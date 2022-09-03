# API

```@meta
DocTestSetup = quote
    using DataFrames, Distributions, SurvivalAnalysis
end
```

### Measures

```@docs
concordance
```

### Surv

```@docs
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
apply_recipe
survival(::SurvivalEstimator)
time
Srv
```

### SurvivalEstimator

```@docs
distr
nelson_aalen
kaplan_meier
std
surv_stats
predict(::SurvivalEstimator, ::DataFrame)
confint(::KaplanMeier, ::Number)
confint(::NelsonAalen, ::Number)
fit(::Type{<:SurvivalEstimator}, ::AbstractMatrix{<:Real}, ::RCSurv)
```

### ParametricSurvival

```@docs
ParametricSurvival
ParametricAFT
ParametricPH
ph
aft
scale
baseline
fit(::Type{<:ParametricSurvival}, ::AbstractMatrix{<:Real}, ::RCSurv, ::Type{T}; ::Number = 1)
predict(::ParametricAFT, ::AbstractMatrix{<:Real})
predict(::ParametricPH, ::AbstractMatrix{<:Real})
```

### Survival distribution functions

```@docs
cum_hazard
survival
hazard
```
