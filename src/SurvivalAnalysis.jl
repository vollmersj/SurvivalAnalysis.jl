module SurvivalAnalysis

    using DataFrames
    using Distributions
    using LinearAlgebra: diag
    using NLSolversBase: hessian!
    using Optim
    using RecipesBase
    using StatsBase
    using StatsModels
    using PrettyTables: pretty_table
    using Crayons

    # reexports
    export DataFrame # DataFrames
    export scale, shape, params, Exponential, Weibull # Distributions
    export coef, confint, stderror, vcov, predict, fit, fit!, std # StatsBase
    export @formula # StatsModels
    export time, reverse, length # Base

    # other exports
    export Surv, outcome_times, event_times, outcome_status, unique_times, threshold_risk
    export unique_event_times
    export total_events, total_censored, total_outcomes, total_risk, surv_stats
    export Srv
    export SurvivalModel
    export kaplan_meier, nelson_aalen, ph, aft, KaplanMeier, NelsonAalen
    export survival, distr
    export ParametricPH, ParametricAFT
    export baseline
    export survival, hazard, cum_hazard, Fₜ, fₜ, pₜ, Sₜ, hₜ, Hₜ
    export SurvivalPrediction
    export concordance, cindex, ConcordanceWeights

    include("tools.jl")
    include("Surv.jl")
    include("SurvTerm.jl")
    include("SurvivalModel.jl")
    include("SurvivalPrediction.jl")
    include("ContinuousUnivariateDistribution.jl")
    include("ParametricSurvival.jl")
    include("SurvivalEstimator.jl")
    include("plot_survivalestimator.jl")
    include("concordance.jl")
end
