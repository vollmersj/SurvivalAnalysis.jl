include("src/Surv.jl")

using Random: seed!
using BenchmarkTools
using Survival

seed!(1)
n = 1000
T = round.(rand(Uniform(1, 10), n));
Δ = rand(Binomial(), n) .== 1;
ot = Surv(T, Δ, "right");
Surv(T, Δ .== 1, "right")
et = EventTime.(T, Δ);

@benchmark confint.(Ref(kaplan(ot)), uniqueEventTimes(ot))
@benchmark confint(fit(Survival.KaplanMeier, et))

@benchmark confint.(Ref(nelson(ot)), uniqueEventTimes(ot))
@benchmark confint(fit(Survival.NelsonAalen, et))
