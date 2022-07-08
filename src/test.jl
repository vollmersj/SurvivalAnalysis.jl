include("src/Surv.jl")

using Random: seed!
using BenchmarkTools
using Survival

seed!(1)
n = 1000
T = round.(rand(Uniform(1, 10), n));
Δ = rand(Binomial(), n) .== 1;
# ot = Surv(T, Δ, "right");
# Surv(T, Δ .== 1, "right")
et = EventTime.(T, Δ);
ot = rcSurv2(T, Δ);

@btime EventTime.(T, Δ);
@btime rcSurv2(T, Δ); ## longer as expected due to postprocessing

@btime Ref(kaplan2(ot)), uniqueEventTimes(ot);
@btime fit(Survival.KaplanMeier, et);
