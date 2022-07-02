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

@btime confint.(Ref(kaplan(ot)), uniqueEventTimes(ot));
@btime confint(fit(Survival.KaplanMeier, et));

map(v -> map(x -> -log(x), v), confint.(Ref(nelson(ot)), uniqueEventTimes(ot)))
confint(fit(Survival.NelsonAalen, et))

@btime nelson(ot);
@btime nelson2(ot);

@assert kaplan(ot).survs == kaplan2(ot).survs
@assert kaplan(ot).sd == kaplan2(ot).sd
@assert nelson2(ot).survs == nelson2(ot).survs
@assert nelson(ot).sd == nelson2(ot).sd

m1 = median(@benchmark kaplan(ot));
m2 = median(@benchmark );
judge(m1, m2)

km =

@assert k1.survs == k2.survs

km.survival

confint(km)
confint.(Ref(k2), uniqueEventTimes(ot))

@btime k1 = kaplan(ot);
@btime k2 = kaplan2(ot);

n1 = @btime nelson(ot);
n2 = @btime nelson2(ot);
@assert n1.survs == n2.survs
@assert n1.sd == n2.sd
n1.sd

m1 = median(@benchmark kaplan($ot))
m2 = median(@benchmark kaplan2($ot))

m3 = median(@benchmark fit(Survival.KaplanMeier, $et))
judge(m1, m3)
judge(m2, m3)

@btime m1
@btime m2
@btime m3
