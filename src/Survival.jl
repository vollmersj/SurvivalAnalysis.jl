module Survival

    using Distributions
    using StatsBase
    using RecipesBase
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

    include("utils.jl")
    include("Surv.jl")
    include("NonParametricEstimator.jl")
    include("KaplanMeier.jl")
    include("NelsonAalen.jl")

end
