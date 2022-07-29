T = [2, 2, 5, 9, 10]
Δ = [true, true, false, false, true]
T2 = T .+ 1
Δ2 = [false, false, true, true, false]

@testset "OneSidedSurv" begin
    # can construct LCSurv
    @test Surv(T, trues(5), "left") isa SurvivalAnalysis.LCSurv
    # can construct RCSurv without censoring
    srv = Surv(T)
    @test srv isa SurvivalAnalysis.RCSurv
    @test srv.status == fill(1, 5)
    @test srv.stats.ncens == zeros(4)
    @test show(srv) === nothing
    # can construct RCSurv from scalar
    @test Surv(1, 0, "right").time == [1]
    # methods
    srv = Surv(T, Δ, "right")
    @test outcome_times(srv) == T
    @test outcome_status(srv) == Δ
    @test event_times(srv) == [2, 2, 10]
    @test unique_times(srv) == [2, 5, 9, 10]
    @test unique_event_times(srv) == [2, 10]
    @test total_events(srv) == 3
    @test total_censored(srv) == 2
    @test total_outcomes(srv) == 5
    @test total_risk(srv) == 5
    @test total_events(srv, 2) == 2
    @test total_events(srv, 5) == 0
    @test total_censored(srv, 2) == 0
    @test total_censored(srv, 5) == 1
    @test total_outcomes(srv, 2) == 2
    @test total_outcomes(srv, 5) == 1
    @test total_risk(srv, 2) == 5
    @test total_risk(srv, 5) == 3
    srv2 = merge(srv, Surv(T2, Δ2, "right"))
    @test srv2.time == [T..., T2...]
    @test srv2.status == [Δ..., Δ2...]
end

@testset "TwoSidedSurv" begin
    srv = Surv(T, T .+ 10)
    @test srv isa SurvivalAnalysis.IntSurv
    @test srv.start == T
    @test srv.stop == T .+ 10
    @test show(srv) === nothing
    @test outcome_times(srv) == [srv.start, srv.stop]
    @test event_times(srv) == T .+ 10
    srv2 = merge(srv, Surv(T2, T2 .+ 10))
    @test srv2.start == [T..., T2...]
    @test srv2.stop == [(T .+ 10)..., (T2 .+ 10)...]
end

true
