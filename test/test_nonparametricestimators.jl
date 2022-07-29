seed!(28349)
n = 50
Δ = fill(1, n);
T = randn(n) .+ 10;
data = DataFrame(status = Δ, time = T)

km = kaplan_meier(@formula(Srv(time, status) ~ 1), data)
na = nelson_aalen(@formula(Srv(time, status) ~ 1), data)

@test kaplan_meier(Surv(T, Δ, "right")) isa SurvivalAnalysis.KaplanMeier
@test nelson_aalen(Surv(T, Δ, "right")) isa SurvivalAnalysis.NelsonAalen

@test km isa
    StatsModels.TableStatisticalModel{SurvivalAnalysis.KaplanMeier, Matrix{Float64}}
@test na isa
    StatsModels.TableStatisticalModel{SurvivalAnalysis.NelsonAalen, Matrix{Float64}}

R"
    library(survival)
    km = survfit(Surv($T) ~ 1)
    T = km$time
    S = km$surv
    H = km$cumhaz
";

R_T = rcopy(R"T");
R_S = rcopy(R"S");
R_H = rcopy(R"H");

@test R_T ≈ time(km) ≈ time(na)
@test R_S ≈ Sₜ(km)
@test R_H ≈ -log.(Sₜ(na))

confint(km) isa Vector{Tuple{Float64, Float64}}
