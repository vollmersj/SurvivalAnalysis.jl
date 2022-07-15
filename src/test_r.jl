using RCall
using Random: seed!
using Distributions
using BenchmarkTools
using Survival
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
";

@benchmark R"
    km = survfit(surv ~ 1)
"

@benchmark kaplan(surv)

R"
    paste(round(km$lower, 2),
    round(km$upper, 2), sep = ',')
"

km = kaplan(surv);
confint(km)
