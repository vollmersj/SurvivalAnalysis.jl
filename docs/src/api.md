# API

```@meta
DocTestSetup = quote
    using DataFrames, Distributions, SurvivalAnalysis
end
```

### Measures

```@docs
ConcordanceWeights
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
apply_recipe(npe::SurvivalEstimator, plot_confint::Bool = true; level = 0.95)
survival(::SurvivalEstimator)
time
Srv
reverse(::SurvivalAnalysis.OneSidedSurv)
```

### SurvivalEstimator

```@docs
SurvivalEstimator
NelsonAalen
KaplanMeier
distr
nelson_aalen
kaplan_meier
std
surv_stats
predict(::SurvivalEstimator, ::DataFrame)
confint(::KaplanMeier, ::Number)
confint(::NelsonAalen, ::Number)
fit(obj::Type{<:SurvivalEstimator}, X::AbstractMatrix{<:Real}, Y::RCSurv)
```

### SurvivalModel

```@docs
SurvivalModel
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
StatsBase.fit(t::Type{<:ParametricSurvival}, X::AbstractMatrix{<:Real}, Y::RCSurv, d::Type{T}; init::Number = 1) where {T <: ContinuousUnivariateDistribution}
predict(fit::ParametricAFT, X::AbstractMatrix{<:Real})
predict(fit::ParametricPH, X::AbstractMatrix{<:Real})
```

### Survival distribution functions

```@docs
cum_hazard
survival
hazard
```

### Prediction

```@docs
SurvivalPrediction
```
