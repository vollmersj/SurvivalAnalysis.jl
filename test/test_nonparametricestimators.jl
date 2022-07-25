seed!(28349)
n = 500
Δ = fill(1, n);
T = randn(n);
data = DataFrame(status = Δ, time = T)

km = kaplan_meier(@formula(Srv(time, status) ~ 1), data)
na = nelson_aalen(@formula(Srv(time, status) ~ 1), data)

km = kaplan_meier(ones((1,1)), Surv(T, Δ, "right"))
km = nelson_aalen(ones((1,1)), Surv(T, Δ, "right"))
