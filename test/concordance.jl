using Random: seed!
using Distributions
using SurvivalAnalysis
using RCall
using Test
seed!(42)
n = 50
ℰ = Exponential()
𝓑 = Bernoulli()
T = rand(ℰ, n);
Δ = rand(𝓑, n);
truth = Surv(T, Δ, :r);
train = Surv(rand(ℰ, n), rand(𝓑, n), :r);
pred = rand(ℰ, n);

@testset "ConcordanceWeights" begin
    @test ConcordanceWeights(1,1,0.5,0.5) isa ConcordanceWeights
    @test_throws ArgumentError ConcordanceWeights(0,0,-1,1)
    @test_throws ArgumentError ConcordanceWeights(0,0,1,2)
end

@testset "Check functions" begin
    @test concordance(truth, pred, :I) isa Number
    @test concordance(truth, pred, :GH) isa Number
end

@testset "Check basic results" begin
    # 0.5 if: all pred same, all times same, no events
    @test concordance(truth, fill(10, n), :I) === 0.5
    @test concordance(Surv(fill(1, n), trues(n), :r), pred, :I) === 0.5
    @test concordance(Surv(randn(n), falses(n), :r), pred, :I) === 0.5

    # perfect prediction
    @test concordance(truth, truth.time, :I) === 1.0
    @test concordance(truth, truth.time, :I, rev=true) === 0.0

    # worst prediction
    @test concordance(truth, 1 .- truth.time, :I) === 0.0
end

@testset "Compare to R" begin
    R"
    library(survival)
    Cng = concordance(Surv($T, $Δ) ~ $pred, timewt = 'n/G')$concordance
    Cng2 = concordance(Surv($T, $Δ) ~ $pred, timewt = 'n/G2')$concordance
    Ci = concordance(Surv($T, $Δ) ~ $pred, timewt = 'I')$concordance

    library(mlr3proba)
    mp_S = mlr3proba:::cindex(Surv($T, $Δ), $pred, weight_meth = 'S', train = Surv($T, $Δ))
    mp_SG = mlr3proba:::cindex(Surv($T, $Δ), $pred, weight_meth = 'SG', train = Surv($T, $Δ))
    mp_GH = mlr3proba:::gonen(sort($pred), 0.5)
    ";

    R_Ci = rcopy(R"Ci");
    R_Cng = rcopy(R"Cng");
    R_Cng2 = rcopy(R"Cng2");

    mp_S = rcopy(R"mp_S")
    mp_SG = rcopy(R"mp_SG")
    mp_GH = rcopy(R"mp_GH");

    J_Ci = concordance(truth, pred, :I, tie_time=0);
    J_Cs = concordance(truth, pred, :S, tie_time=0);
    J_Csg = concordance(truth, pred, :SG, tie_time=0);
    J_Cng = concordance(truth, pred, :G, tie_time=0);
    J_Cng2 = concordance(truth, pred, :G2, tie_time=0);
    J_Cgh = concordance(truth, pred, :GH, tie_time=0);

    @test J_Ci ≈ R_Ci
    @test J_Cng ≈ R_Cng
    @test J_Cng2 ≈ R_Cng2

    # different implementation to (or not implemented in) {survival}
    @test round.(J_Csg, digits=5) ≈ round.(1 - mp_SG, digits=5)
    @test J_Cgh ≈ 1 - mp_GH
    @test J_Cs ≈ 1 - mp_S
end

true
