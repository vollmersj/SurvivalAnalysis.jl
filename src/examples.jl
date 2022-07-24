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
surv_stats(srv, events_only = false)
## get cohort statistics for events only
surv_stats(srv, events_only = true)

outcome_times(srv) # times at which outcomes occur (event or censoring)
unique_times(srv) # unique times at which outcomes occur (event or censoring)
event_times(srv) # (unsorted) times at which events occur
unique_event_times(srv) # unique times at which events occur
outcome_status(srv) # status at outcome (event (1) or censoring (0))

total_outcomes(srv, 6) # total deaths at time T
total_outcomes(srv, 0) # outcomes/censored/events return 0 if T = 0 for outcomes
total_censored(srv, 6) # total deaths at time T
total_events(srv, 6) # total deaths at time T
total_events(srv, 10) # all return 0 if T > max(T)
total_risk(srv, 6) # total deaths at time T
total_risk(srv, 0) # return risk set at min(T) if T=0

km = kaplan(srv) # kaplan-Meier
na = nelson(srv) # Nelson-Aalen


confint(km)
# na = nelson(survs)

using Plots
plot(km)
survival.(Ref(km), [1, 4, 7, 10])
