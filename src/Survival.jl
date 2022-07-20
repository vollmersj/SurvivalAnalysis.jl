module Survival

    using Distributions
    using StatsBase
    using RecipesBase
    using DataFrames
    using LinearAlgebra: diag
    using Optim
    # using LineSearches, Optim, NLSolversBase
    # using LineSearches, Optim, NLSolversBase
    # using Plots

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
        nelson

    abstract type FittedParametric end

    include("utils.jl")
    include("Surv.jl")
    include("SurvivalPrediction.jl")
    include("ParametricPH.jl")
    include("ParametricAFT.jl")
    include("FittedParametricAFT.jl")
    include("FittedParametricPH.jl")
    include("NonParametricEstimator.jl")
    include("KaplanMeier.jl")
    include("NelsonAalen.jl")

end
