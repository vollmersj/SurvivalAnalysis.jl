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

@testset "Concordance functions" begin
    @test concordance(truth, pred, :I) isa SurvivalAnalysis.Concordance
    @test concordance(truth, pred, :I).C isa Number
    @test concordance(truth, pred, :GH) isa SurvivalAnalysis.Concordance
end

@testset "Concordance fails as expected" begin
    @test_throws ArgumentError concordance(truth, pred, :I, tied_preds = 2)
    @test_throws ArgumentError concordance(truth, pred, :P)
    @test_throws ArgumentError concordance(truth, [1], :I)
end

@testset "Check fallback results" begin
    # 0.5 if: all pred same, all times same, no events
    @test concordance(truth, fill(10, n), :I).C === 0.5
    @test concordance(truth, fill(10, n), :I).weights.name ===
        "Fallback (all predictions equal)"
    @test concordance(Surv(fill(1, n), trues(n), :r), pred, :I).C === 0.5
    @test concordance(Surv(fill(1, n), trues(n), :r), pred, :I).weights.name ===
        "Fallback (all observed event times equal)"
    @test concordance(Surv(randn(n), falses(n), :r), pred, :I).C === 0.5
    @test concordance(Surv(randn(n), falses(n), :r), pred, :I).weights.name ===
        "Fallback (no observed events)"
end

@testset "Check names" begin
    @test concordance(truth, pred, :I).weights.name === "Harrell's"
    @test concordance(truth, pred, :Harrell).weights.name === "Harrell's"
    @test concordance(truth, pred, :G).weights.name === ""
    @test concordance(truth, pred, :G2).weights.name === "Uno's"
    @test concordance(truth, pred, :Uno).weights.name === "Uno's"
    @test concordance(truth, pred, :GH).weights.name === "Gönen-Heller's"
    @test concordance(truth, pred, :Gonen).weights.name === "Gönen-Heller's"
    @test concordance(truth, pred, :SG).weights.name === "Schemper's"
    @test concordance(truth, pred, :Schemper).weights.name === "Schemper's"
    @test concordance(truth, pred, :S).weights.name === "Peto-Wilcoxon's"
    @test concordance(truth, pred, :Peto).weights.name === "Peto-Wilcoxon's"
end

@testset "Check basic results" begin
    # perfect prediction
    @test concordance(truth, truth.time, :I).C === 1.0
    @test concordance(truth, truth.time, :I, rev=true).C === 0.0

    # worst prediction
    @test concordance(truth, 1 .- truth.time, :I).C === 0.0

    # show returns nothing
    @test concordance(truth, 1 .- truth.time, :I).show === nothing

    # overwrite as expected
    C = concordance(truth, pred, :I, tied_times = 0.1, tied_preds = 0.9)
    @test C.weights.tied_times == Float16(0.1)
    @test C.weights.tied_preds == Float16(0.9)
    @test C.weights.G == 0
    @test C.weights.S == 0

    C = concordance(truth, pred, :I, tied_times = 0.1, tied_preds = 0.9,
        custom_weights = ConcordanceWeights(1, 2, 0.2, 0.8, ""))
    @test C.weights.tied_times == Float16(0.8)
    @test C.weights.tied_preds == Float16(0.2)
    @test C.weights.G == 2
    @test C.weights.S == 1
end

## comparing to {survival} is quick as it's bundled with R but
##  checking to other packages takes a very long time
## (e.g. >5min install for mlr3proba on CI) and has a high failure
## rate -- therefore only test outputs from some indices below --
## the others have been tested locally against mlr3proba (and are equal)
@testset "Compare to R" begin
    R"
    library(survival)
    Cng = concordance(Surv($T, $Δ) ~ $pred, timewt = 'n/G')$concordance
    Cng2 = concordance(Surv($T, $Δ) ~ $pred, timewt = 'n/G2')$concordance
    Ci = concordance(Surv($T, $Δ) ~ $pred, timewt = 'I')$concordance
    ";

    R_Ci = rcopy(R"Ci");
    R_Cng = rcopy(R"Cng");
    R_Cng2 = rcopy(R"Cng2");

    J_Ci = concordance(truth, pred, :I, tied_times=0).C;
    J_Cs = concordance(truth, pred, :S, tied_times=0).C;
    J_Csg = concordance(truth, pred, :SG, tied_times=0).C;
    J_Cng = concordance(truth, pred, :G, tied_times=0).C;
    J_Cng2 = concordance(truth, pred, :G2, tied_times=0).C;
    J_Cgh = concordance(truth, pred, :GH, tied_times=0).C;

    @test J_Ci ≈ R_Ci
    @test J_Cng ≈ R_Cng
    @test J_Cng2 ≈ R_Cng2

    # different implementation to (or not implemented in) {survival}
    # see note above for why tests are reduced
    @test J_Csg isa Number
    @test J_Cgh isa Number
    @test J_Cs isa Number
end

@testset "ties work as expected" begin
    Y = Surv([1,2,3], trues(3), :r)
    # compatible pairs: [1, 2], [1, 3], [2, 3]
    # tied preds
    # when tied_preds=0 only compatible with crank different
    @test concordance(Y, [2, 2, 3], :I; tied_preds = 0).C === 1.0
    # concordant pairs: [0.5, 1, 1] = 2.5
    @test concordance(Y, [2, 2, 3], :I; tied_preds = 0.5).C === 2.5/3
    # concordant pairs: [1, 1, 1] = 3
    @test concordance(Y, [2, 2, 3], :I; tied_preds = 1).C === 1.0

    # tied times
    Y = Surv([1,1,2], trues(3), :r)
    p = [1,2,3]
    # compatible pairs: [1, 2] [1, 2]
    # concordant pairs: [1, 1] = 2
    @test concordance(Y, p, :I; tied_times = 0).C === 1.0
    # compatible pairs: [1, 1] [1, 2] [1, 2]
    # concordant pairs: [0.5, 1, 1] = 2.5
    @test concordance(Y, p, :I; tied_times = 0.5).C === 2.5/3
    # compatible pairs: [1, 1] [1, 2] [1, 2]
    # concordant pairs: [1, 1, 1] = 3
    @test concordance(Y, p, :I; tied_times = 1).C === 1.0
end

@testset "count concordant as expected" begin
    Y = Surv([1, 1, 2, 3, 4, 4], trues(6), :r)
    p = [1, 2, 3, 3, 2.5, 2.5]

    C = concordance(Y, p, :I);
    @test C.pairs == 15
    @test C.tied_times === 1
    @test C.tied_preds === 1
    @test C.tied_both === 1
    @test C.comparable === 13
    @test C.concordant === 8
    @test C.disconcordant === 4

    C = concordance(Y, p, :GH);
    @test C.pairs == 15
    @test C.tied_times === 0
    @test C.tied_preds === 2
    @test C.tied_both === 0
    @test C.comparable === 15
    @test C.concordant === 13
    @test C.disconcordant === 2

    C = concordance(Y, p, :I, tied_preds = 1, tied_times = 1);
    @test C.pairs == 15
    @test C.tied_times === 1
    @test C.tied_preds === 1
    @test C.tied_both === 1
    @test C.comparable === 15
    @test C.concordant === 8
    @test C.disconcordant === 4
end

true
