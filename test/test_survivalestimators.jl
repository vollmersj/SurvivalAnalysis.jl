seed!(28349)
n = 50
Δ = fill(1, n);
T = randn(n) .+ 10;
data = DataFrame(status = Δ, time = T)

km = kaplan_meier(@formula(Srv(time, status) ~ 1), data)
na = nelson_aalen(@formula(Srv(time, status) ~ 1), data)

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

confint(km) isa Vector{Tuple{Float64, Float64}} # TODO - Improve precision

p = predict(km, DataFrame(x = randn(10)))

@test p isa SurvivalAnalysis.DiscreteSurvivalPrediction{Float64}
@test all(p.lp .=== p.crank .=== p.time .=== fill(NaN, 10))
@test length(unique(p.distr)) == 1
@test p.distr[1] == distribution(km) # FIXME - FIX CONSISTENCY DISTR/DISTRIUTION
@test p.survival_matrix.time == time(km)
@test size(p.survival_matrix.surv) == (10, length(time(km)))
@test unique(p.survival_matrix.surv) == survival(km)  # FIXME - FIX CONSISTENCY SURV/SURVIVAL


# test non-formula interface
na = nelson_aalen(Surv(T, Δ, "right"))
@test kaplan_meier(Surv(T, Δ, "right")) isa SurvivalAnalysis.KaplanMeier
@test na isa SurvivalAnalysis.NelsonAalen

p = predict(na, DataFrame(x = randn(10)))

@test p isa SurvivalAnalysis.DiscreteSurvivalPrediction{Float64}
@test all(p.lp .=== p.crank .=== p.time .=== fill(NaN, 10))
@test length(unique(p.distr)) == 1
@test p.distr[1] == distribution(na) # FIXME - FIX CONSISTENCY DISTR/DISTRIUTION
@test p.survival_matrix.time == time(na)
@test size(p.survival_matrix.surv) == (10, length(time(na)))
@test unique(p.survival_matrix.surv) == survival(na)  # FIXME - FIX CONSISTENCY SURV/SURVIVAL

true
