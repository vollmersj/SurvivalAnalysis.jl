n = 5

@testset "Srv" begin
    data = DataFrame(Δ = trues(n), T = randn(n), X = randn(n))
    @test_throws ArgumentError apply_schema(@formula(Srv(T, Δ, 0) ~ 1), schema(data))
    @test modelcols(apply_schema(@formula(Srv(T, Δ) ~ 1), schema(data)), data)[1] isa
        SurvivalAnalysis.RCSurv
    @test modelcols(apply_schema(@formula(Srv(T, Δ, 1) ~ 1), schema(data)), data)[1] isa
        SurvivalAnalysis.RCSurv
    @test modelcols(apply_schema(@formula(Srv(T, Δ, -1) ~ 1), schema(data)), data)[1] isa
        SurvivalAnalysis.LCSurv
end

@testset "SurvTerm" begin
    st = SurvivalAnalysis.SurvTerm(Term(:X), Term(:Y), ConstantTerm(1))
    @test st isa SurvivalAnalysis.SurvTerm
    @test show(st) === nothing
    @test SurvivalAnalysis.SurvTerm(Term(:X), Term(:Y)) == st
end
