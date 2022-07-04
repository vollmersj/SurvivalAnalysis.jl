include("src/Surv.jl")

n = 10
T = round.(rand(Uniform(1, 10), n));
Δ = rand(Binomial(), n) .== 1;
srv = Surv(T, Δ, "right");
srv2 = rcSurv2(T, Δ);

@btime rcSurv(T, Δ); kaplan(srv);
@btime rcSurv2(T, Δ); kaplan2(srv2);

@btime uniqueTimes(srv);
@btime uniqueTimes(srv2);

@btime totalEvents(srv, 8);
@btime totalEvents(srv2, 8);

@btime totalCensored(srv, 8);
@btime totalCensored(srv2, 8);

@btime totalRisk(srv, 8);
@btime totalRisk(srv2, 8);

@btime rcSurv(T, Δ); tabulateRisk(srv, events = false);
@btime rcSurv2(T, Δ); tabulateRisk(srv2, events = false);

@btime kaplan(srv);
@btime kaplan2(srv2);

kaplan(srv)
kaplan2(srv)
@btime kaplan(srv);
@btime kaplan2(srv);
@assert kaplan(srv) == kaplan2(srv)

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
