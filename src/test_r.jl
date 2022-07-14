import Survival
using RCall
using Random: seed!
using Distributions
seed!(1)
n = 1000
T = round.(rand(Uniform(1, 10), n));
Δ = rand(Binomial(), n) .== 1;


R"
    library(survival)
    time = $T
    status = $Δ
    surv = Surv(time, status)
    km = survfit(surv ~ 1)
    paste(round(km$lower, 2),
    round(km$upper, 2), sep = ',')
"

surv = Survival.Surv(T, Δ, "right");
km = Survival.kaplan(surv);
km.survival
Survival.confint(km)
