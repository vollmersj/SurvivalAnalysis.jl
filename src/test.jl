using StatsBase
using RecipesBase
using BenchmarkTools

using Distributions

include("src/utils.jl")
include("src/Surv.jl")
include("src/NonParametricEstimator.jl")
include("src/KaplanMeier.jl")
include("src/NelsonAalen.jl")

using Survival

n = 1000
T = round.(rand(Uniform(1, 10), n));
Δ = rand(Binomial(), n) .== 1;
et = EventTime.(T, Δ);
srv = Surv(T, Δ, "right");

@btime EventTime.(T, Δ);
@btime Surv(T, Δ, "right"); ## longer as expected due to postprocessing

@btime (k = kaplan(srv); confint(k));
@btime (k = fit(Survival.KaplanMeier, et); confint(k));

using Plots
plot(kaplan(srv))
