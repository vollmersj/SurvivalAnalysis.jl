module SurvivalAnalysis

    using DataFrames: DataFrame, groupby
    using Distributions
    using LinearAlgebra: diag, pinv
    using NLSolversBase: hessian!
    using Optim
    using RecipesBase
    using StatsBase
    using StatsModels
    using PrettyTables: pretty_table
    using Tables
    using Statistics: mean

    ## documented
    export Surv, outcome_times, event_times, outcome_status, unique_outcome_times, threshold_risk
    export unique_event_times
    export total_events, total_censored, total_outcomes, total_risk, surv_stats
    export Srv
    export SurvivalModel, SurvivalEstimator
    export kaplan_meier, nelson_aalen, ph, aft, KaplanMeier, NelsonAalen
    export survival, distr
    export ParametricSurvival, ParametricPH, ParametricAFT
    export baseline
    export survival, hazard, cum_hazard
    export SurvivalPrediction
    export SurvivalMeasure, concordance, ConcordanceWeights
    export SurvivalTimeMeasure, MSE, RMSE, MAE
    export logrank

    ## undocumented
    # reexports
    export DataFrame # DataFrames
    export scale, shape, params, Exponential, Weibull # Distributions
    export coef, confint, stderror, vcov, predict, fit, fit!, std # StatsBase
    export @formula # StatsModels
    export time, reverse, length # Base

    # aliases
    export Fₜ, fₜ, pₜ, Sₜ, hₜ, Hₜ
    export km, na, kaplan, nelson
    export cindex

    include("tools.jl")
    include("Surv.jl")
    include("SurvTerm.jl")
    include("SurvivalModel.jl")
    include("SurvivalPrediction.jl")
    include("ContinuousUnivariateDistribution.jl")
    include("ParametricSurvival.jl")
    include("SurvivalEstimator.jl")
    include("plots.jl")
    include("SurvivalMeasure.jl")
    include("SurvivalTimeMeasure.jl")
    include("concordance.jl")
end
