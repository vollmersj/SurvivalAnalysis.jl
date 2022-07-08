include("src/Surv.jl")

n = 10
T = round.(rand(Uniform(1, 10), n));
Δ = rand(Binomial(), n) .== 1;
## usually will initialise vectors
srv = Surv(T, Δ, "right");
# srv2 = Surv(T .+ 12, Δ, "right");
# isrv = Surv(T, T .+ 21);
# isrv2 = Surv(T, T .+ 10);

## can merge into new object
# merge(srv, srv2)
# merge(isrv, isrv2)

# ## but can also initialise as scalar (but this is a weird use-case!)
# Surv(1, 1, "right")
# Surv(5, 0, "right")

## get cohort statistics for all outcomes
survStats(srv, events_only = false)
## get cohort statistics for events only
survStats(srv, events_only = true)

outcomeTimes(srv) # times at which outcomes occur (event or censoring)
uniqueTimes(srv) # unique times at which outcomes occur (event or censoring)
eventTimes(srv) # (unsorted) times at which events occur
uniqueEventTimes(srv) # unique times at which events occur
outcomeStatus(srv) # status at outcome (event (1) or censoring (0))

totalOutcomes(srv, 6) # total deaths at time T
totalOutcomes(srv, 0) # outcomes/censored/events return 0 if T = 0 for outcomes
totalCensored(srv, 6) # total deaths at time T
totalEvents(srv, 6) # total deaths at time T
totalEvents(srv, 10) # all return 0 if T > max(T)
totalRisk(srv, 6) # total deaths at time T
totalRisk(srv, 0) # return risk set at min(T) if T=0

km = kaplan(srv) # kaplan-Meier
na = nelson(srv) # Nelson-Aalen


confint(km)
# na = nelson(survs)

using Plots
plot(km)
survival.(Ref(km), [1, 4, 7, 10])
