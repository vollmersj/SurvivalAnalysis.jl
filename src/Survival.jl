module Survival

    using RecipesBase
    using DataFrames
    using LinearAlgebra: diag
    using Optim
    using NLSolversBase: hessian!
    using Reexport
    @reexport using StatsModels
    using Distributions
    using StatsBase
    export coef, confint, stderror, vcov, predict, fit, fit!

    export
        Surv,
        outcomeTimes,
        eventTimes,
        outcomeStatus,
        uniqueTimes,
        uniqueEventTimes,
        totalEvents,
        totalCensored,
        totalOutcomes,
        totalRisk,
        survStats,
        kaplan,
        nelson,
        ParametricPH,
        ParametricAFT,
        ph,
        aft

    abstract type FittedParametric <: StatisticalModel end

    include("utils.jl")
    include("Surv.jl")
    include("SurvivalPrediction.jl")
    include("ContinuousPHDistribution.jl")
    include("ContinuousAFTDistribution.jl")
    include("FittedParametricAFT.jl")
    include("FittedParametricPH.jl")
    include("NonParametricEstimator.jl")
    include("KaplanMeier.jl")
    include("NelsonAalen.jl")

end
