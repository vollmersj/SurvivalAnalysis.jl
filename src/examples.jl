include("src/Surv.jl")

n = 10
T = round.(rand(Uniform(1, 10), n));
Δ = rand(Binomial(), n) .== 1;
srv = Surv(T, Δ, "right");

outcomeTimes(srv)

Surv(1.0, false, "left")
Surv(1.0, true, "left")
Surv(1.0, false, "right")
Surv(1.0, false, "interval")
Surv(1.0, 2.0)

using Distributions

Surv.(rand(Uniform(), 10), rand(Binomial(1, 0.5), 10) .== 1, "right")
Surv.(rand(Uniform(), 10), rand(Binomial(1, 0.5), 10) .== 1, "left")
Surv.(rand(Uniform(), 10), rand(Uniform(), 10))


Surv.(round.(rand(Uniform(1, 5), 10)), round.(rand(Uniform(1, 5), 10)))
survs = Surv.(round.(rand(Uniform(1, 5), 50)), rand(Binomial(1, 0.5), 50) .== 1, "right")

outcomeTimes(survs)
outcomeStatus(survs)
uniqueTimes(survs)

km = kaplan(survs)
confit(km, 2)
na = nelson(survs)

using Plots
plot(km.times, km.survs, linetype=:steppost, ylims = (0, 1))
plot(na.times, na.survs, linetype=:steppost, ylims = (0, 1))

survival.(Ref(km), [1, 4, 7])
cdf.(Ref(km), [1, 4, 7])
chf.(Ref(km), [1, 4, 7])





plot(km)
