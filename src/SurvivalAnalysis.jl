module SurvivalAnalysis

    using DataFrames
    using Distributions
    using LinearAlgebra: diag
    using NLSolversBase: hessian!
    using Optim
    using RecipesBase
    using StatsBase
    using StatsModels

    export DataFrame # DataFrames
    export scale, shape, params, Exponential, Weibull # Distributions
    export coef, confint, stderror, vcov, predict, fit, fit! # StatsBase
    export @formula, terms, coefnames # StatsModels
    export Surv, outcome_times, event_times, outcome_status, unique_times, unique_event_times
    export total_events, total_censored, total_outcomes, total_risk, surv_stats
    export Srv
    export kaplan, nelson, ph, aft
    export ParametricPH, ParametricAFT
    export baseline

    include("utils.jl")
    include("Surv.jl")
    include("SurvTerm.jl")
    include("SurvivalModel.jl")
    include("SurvivalPrediction.jl")
    include("ContinuousPHDistribution.jl")
    include("ContinuousAFTDistribution.jl")
    include("ParametricSurvival.jl")
    include("ParametricAFT.jl")
    include("ParametricPH.jl")
    include("NonParametricEstimator.jl")
    include("KaplanMeier.jl")
    include("NelsonAalen.jl")
end
