seed!(28349)
n = 50
Δ = fill(1, n);
T = randn(n) .+ 10;
data = DataFrame(status = Δ, time = T, X = randn(n))

@testset "expandstats test" begin

    T1 = Float64[3, 1, 5, 9, 7]
    Δ1 = [1, 1, 0, 1, 1]
    T2 = Float64[2, 10, 4, 8, 6, 12]
    Δ2 = [1, 0, 0, 1, 0, 1]
    S1 = Surv(T1, Δ1, :r)
    S2 = Surv(T2, Δ2, :r)
    M = merge(S1, S2)

    stats1 = SurvivalAnalysis.expandstats(S1.stats, unique_outcome_times(M))
    @test isapprox(stats1.time, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12])
    @test isapprox(stats1.nrisk, [5, 4, 4, 3, 3, 2, 2, 1, 1, 0, 0])
    @test isapprox(stats1.ncens, [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0])
    @test isapprox(stats1.nevents, [1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0])
    @test isapprox(stats1.noutcomes, [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0])
end

@testset "logrank test (unstratified)" begin
    R"
    library(survival)
    t = c(3, 3, 5, 9, 7, 2, 10, 4, 6, 6, 10, 13, 15, 15, 14)
    s = c(1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0)
    x = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3)
    d = data.frame(t=t, s=s, x=x)
    r3 = survdiff(Surv(t, s) ~ x, d, rho=0)
    d2 = d[d$x != 3,]
    r2 = survdiff(Surv(t, s) ~ x, d2, rho=0)
    ";

    T1 = Float64[3, 3, 5, 9, 7]
    Δ1 = [1, 1, 0, 1, 1]
    T2 = Float64[2, 10, 4, 6, 6, 10]
    Δ2 = [1, 0, 0, 1, 0, 1]
    T3 = Float64[13, 15, 15, 14]
    Δ3 = [1, 1, 0, 0]
    S1 = Surv(T1, Δ1, :r)
    S2 = Surv(T2, Δ2, :r)
    S3 = Surv(T3, Δ3, :r)
    T = vcat(T1, T2, T3)
    Δ = vcat(Δ1, Δ2, Δ3)
    G = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3]

    # Test logrank statistics against R
    r02j = logrank(S1, S2; wtmethod=:logrank)
    r03j = logrank(S1, S2, S3; wtmethod=:logrank)
    r03jx = logrank(T, Δ, G, zeros(0); wtmethod=:logrank)
    r02rstat = rcopy(R"r2$chisq")
    r03rstat = rcopy(R"r3$chisq")
    @test isapprox(r02rstat, r02j.stat)
    @test isapprox(r02j.dof, 1)
    @test isapprox(r03rstat, r03j.stat)
    @test isapprox(r03j.dof, 2)
    @test isapprox(r03rstat, r03jx.stat)
    @test isapprox(r03jx.dof, 2)

    # Test Wilcoxon statistics against Stata
    r12j = logrank(S1, S2; wtmethod=:wilcoxon)
    r13j = logrank(S1, S2, S3; wtmethod=:wilcoxon)
    @test isapprox(5.2944871, r13j.stat)
    @test isapprox(2, r13j.dof)
    @test isapprox(0.5540201, r12j.stat)
    @test isapprox(1, r12j.dof)

    # Test Tarone-Ware statistics against Stata
    r12j = logrank(S1, S2; wtmethod=:tw)
    r13j = logrank(S1, S2, S3; wtmethod=:tw)
    @test isapprox(6.4063107, r13j.stat)
    @test isapprox(2, r13j.dof)
    @test isapprox(0.88063788, r12j.stat)
    @test isapprox(1, r12j.dof)

    # Test Peto-Peto statistics against Stata
    r12j = logrank(S1, S2; wtmethod=:peto)
    r13j = logrank(S1, S2, S3; wtmethod=:peto)
    @test isapprox(6.029332, r13j.stat)
    @test isapprox(2, r13j.dof)
    @test isapprox(0.61770817, r12j.stat)
    @test isapprox(1, r12j.dof)
end

@testset "logrank test (stratified)" begin
    R"
    library(survival)
    t = c(3, 3, 5, 9, 7, 2, 10, 4, 6, 6, 10, 13, 15, 15, 14)
    s = c(1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0)
    g = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3)
    z = c(1, 1, 1, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 2, 2)
    d = data.frame(t=t, s=s, g=g, z=z)
    r = survdiff(Surv(t, s) ~ g + strata(z), d, rho=0)
    ";

    T = Float64[3, 3, 5, 9, 7, 2, 10, 4, 6, 6, 10, 13, 15, 15, 14]
    Δ = [1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0]
    G = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3]
    Z = [1, 1, 1, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 2, 2]

    # Test logrank statistics against R
    rj = logrank(T, Δ, G, Z; wtmethod=:logrank)
    rstat = rcopy(R"r$chisq")
    @test isapprox(rstat, rj.stat)
    @test isapprox(rj.dof, 2)
end

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
    @test all(confint(km) .=== confint(km.model))
    @test std(km) isa Vector{Float64}
    @test std(km) === std(km.model)

    @test confint(na) isa Vector{Tuple{Float64, Float64}} # TODO - Improve precision & test
    @test all(confint(na) .=== confint(na.model))
    @test std(na) isa Vector{Float64}
    @test std(na) === std(na.model)

    @test show(km) === nothing
    @test show(na.model) === nothing
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
    na = nelson_aalen(Surv(T, Δ, :r))
    @test kaplan_meier(Surv(T, Δ, :r)) isa SurvivalAnalysis.KaplanMeier
    @test fit(KaplanMeier, Surv(T, Δ, :r)) isa SurvivalAnalysis.KaplanMeier
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
