include("src/Surv.jl")

using Random: seed!
using BenchmarkTools
using Survival

seed!(1)
n = 100
T = round.(rand(Uniform(1, 10), n));
Δ = rand(Binomial(), n) .== 1;
ot = Surv(T, Δ, "right");
Surv(T, Δ .== 1, "right")
et = EventTime.(T, Δ);

k1 = kaplan(ot)
k2 = kaplan2(ot)

km = fit(Survival.KaplanMeier, et)

@assert k1.survs == k2.survs
km.survival


m1 = median(@benchmark kaplan($ot))
m2 = median(@benchmark kaplan2($ot))
m3 = median(@benchmark fit(Survival.KaplanMeier, $et))
judge(m1, m3)
judge(m2, m3)

@btime m1
@btime m2
@btime m3
