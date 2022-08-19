seed!(28349)
n = 50
Δ = fill(1, n);
T = randn(n) .+ 10;
data = DataFrame(status = Δ, time = T, X = randn(n))

@testset "Can fit/predict with covariates" begin
    @test predict(kaplan_meier(@formula(Srv(time, status) ~ X), data), data) isa
        SurvivalAnalysis.DiscreteSurvivalPrediction{Float64}
end

km = kaplan_meier(@formula(Srv(time, status) ~ 1), data)
na = nelson_aalen(@formula(Srv(time, status) ~ 1), data)

@testset "Basics" begin
    @test km isa
        StatsModels.TableStatisticalModel{SurvivalAnalysis.KaplanMeier, Matrix{Float64}}
    @test na isa
        StatsModels.TableStatisticalModel{SurvivalAnalysis.NelsonAalen, Matrix{Float64}}

    @test confint(km) isa Vector{Tuple{Float64, Float64}} # TODO - Improve precision & test
    @test confint(km) == confint(km.model)
    @test std(km) isa Vector{Float64}
    @test std(km) === std(km.model)

    @test confint(na) isa Vector{Tuple{Float64, Float64}} # TODO - Improve precision & test
    @test confint(na) == confint(na.model)
    @test std(na) isa Vector{Float64}
    @test std(na) === std(na.model)
end

@testset "Alignment with R" begin
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
end

@testset "Can make predictions" begin
    p = predict(km, DataFrame(x = randn(10)))

    @test p isa SurvivalAnalysis.DiscreteSurvivalPrediction{Float64}
    @test all(p.lp .=== p.crank .=== p.time .=== fill(NaN, 10))
    @test length(unique(p.distr)) == 1
    @test p.distr[1] == distr(km)
    @test p.survival_matrix.time == time(km)
    @test size(p.survival_matrix.survival) == (10, length(time(km)))
    @test unique(p.survival_matrix.survival) == survival(km)
end

@testset "Non-formula interface works" begin
    # test non-formula interface
    na = nelson_aalen(Surv(T, Δ, "right"))
    @test kaplan_meier(Surv(T, Δ, "right")) isa SurvivalAnalysis.KaplanMeier
    @test fit(KaplanMeier, Surv(T, Δ, "right")) isa SurvivalAnalysis.KaplanMeier
    @test na isa SurvivalAnalysis.NelsonAalen

    p = predict(na, DataFrame(x = randn(10)))

    @test p isa SurvivalAnalysis.DiscreteSurvivalPrediction{Float64}
    @test all(p.lp .=== p.crank .=== p.time .=== fill(NaN, 10))
    @test length(unique(p.distr)) == 1
    @test p.distr[1] == distr(na)
    @test p.survival_matrix.time == time(na)
    @test size(p.survival_matrix.survival) == (10, length(time(na)))
    @test unique(p.survival_matrix.survival) == survival(na)
end

true
