using RCall
using Random: seed!
using Distributions
using BenchmarkTools
using SurvivalAnalysis
using StatsBase
seed!(1)
n = 1000
T = round.(rand(Uniform(1, 10), n));
Δ = rand(Binomial(), n) .== 1;

surv = Surv(T, Δ, "right");

R"
    library(survival)
    time = $T
    status = $Δ
    surv = Surv(time, status)
    km = survfit(surv ~ 1)
    paste(round(km$lower, 2),
    round(km$upper, 2), sep = ',')
"

km = kaplan_meier(@formula(Srv(time, status) ~ 1), data)
km = kaplan_meier(surv);
confint(km)
