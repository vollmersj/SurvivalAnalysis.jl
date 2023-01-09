function npmle_gendat(n, dist, tdist, rng)
    event, start, stop, ltrunc = Float64[], Float64[], Float64[], Float64[]
    for i = 1:n
        while true
            evt = rand(rng, dist)
            sta = floor.(10 * evt) / 10
            sto = sta .+ 0.1 .- 1e-10

            # Left truncation
            ltr = rand(rng, tdist)
            if sta >= ltr
                push!(event, evt)
                push!(ltrunc, ltr)
                push!(start, sta)
                push!(stop, sto)
                break
            end
        end
    end

    return event, Surv(start, stop, ltrunc, Float64[])
end

@testset "Check the gradient of the NPMLE for interval censored data" begin

    rng = StableRNG(123)

    n = 100
    for k in 1:1
        event, Y = npmle_gendat(n, Exponential(1), Exponential(1), rng)
        ms = HazardNPMLE(Y)
        p = length(ms.support)
        par = 3 .+ rand(rng, p-1)
        sort!(par)
        agrad = zeros(p-1)
        SurvivalAnalysis.score!(ms, agrad, par)
        loglike = par -> SurvivalAnalysis.loglike(ms, par)
        ngrad = grad(central_fdm(8, 1), loglike, par)[1]
        @test isapprox(agrad, ngrad)
    end
end

@testset "Check monotone hazard NPMLE fit for interval censored data" begin

    rng = StableRNG(123)
    n = 250
    nrep = 2

    # The distributions and sample sizes here match Pan and Chappell (1998).
    for (dist, tdist, mv, tol) in [[Exponential(1), Uniform(0, 1.5), 4., 0.1],
                                   [Weibull(4, 1), Uniform(0, 1.5), 2., 0.02]]
        xs = collect(range(0., mv, 20))
        haz = zeros(length(xs), nrep)
        thaz = [hazard(dist, x) for x in xs]
        surv = zeros(length(xs), nrep)
        for k in 1:nrep
            _, Y = npmle_gendat(n, dist, tdist, rng)
            ms = fit(HazardNPMLE, Y)
            @test Optim.converged(ms.optim_rslt) == true
            haz[:, k] = [hazard(ms, x) for x in xs]
            surv[:, k] = [survival(ms, x) for x in xs]
        end

        ehaz = mean(haz, dims=2)[:]
        ii = findall(thaz .< 2)
        irmse = sqrt(mean((ehaz[ii] .- thaz[ii]).^2))
        @test irmse < 0.3

        esurv = mean(surv, dims=2)[:]
        tsurv = [1 - cdf(dist, x) for x in xs]
        irmse = sqrt(mean((esurv .- tsurv).^2))
        @test irmse < 0.1
    end
end

@testset "Check monotone hazard NPMLE fit for interval censored data with weights" begin

    rng = StableRNG(123)
    n = 250
    m = div(n, 2)

    _, Y1 = npmle_gendat(n, Exponential(1), Uniform(0, 1.5), rng)
    w = ones(n)
    w[m:end] .= 2.0
    Y1 = Surv(Y1.start, Y1.stop, Y1.ltrunc, w)

    start2 = vcat(Y1.start, Y1.start[m:end])
    stop2 = vcat(Y1.stop, Y1.stop[m:end])
    ltrunc2 = vcat(Y1.ltrunc, Y1.ltrunc[m:end])
    Y2 = Surv(start2, stop2, ltrunc2, Float64[])
    Y3 = Surv(start2, stop2, ltrunc2, ones(length(start2)))

    opts = Optim.Options(iterations=100000, g_tol=1e-12)
    ms1 = fit(HazardNPMLE, Y1; optim_opts=opts)
    ms2 = fit(HazardNPMLE, Y2; optim_opts=opts)
    ms3 = fit(HazardNPMLE, Y3; optim_opts=opts)

    @test Optim.converged(ms1.optim_rslt) == true
    @test Optim.converged(ms2.optim_rslt) == true
    @test Optim.converged(ms3.optim_rslt) == true

    @test isapprox(SurvivalAnalysis.loglike(ms1, ms1.par), SurvivalAnalysis.loglike(ms2, ms1.par))
    @test isapprox(SurvivalAnalysis.loglike(ms1, ms1.par), SurvivalAnalysis.loglike(ms3, ms1.par))

    ll1 = SurvivalAnalysis.loglike(ms1, ms1.par)
    ll2 = SurvivalAnalysis.loglike(ms2, ms2.par)
    ll3 = SurvivalAnalysis.loglike(ms3, ms3.par)
    ll4 = SurvivalAnalysis.loglike(ms3, ms1.par)
    ll5 = SurvivalAnalysis.loglike(ms1, ms3.par)
    ll = [ll1, ll2, ll3, ll4, ll5]
    @test maximum(ll) - minimum(ll) < 0.5

    G1 = zeros(length(ms1.par))
    G2 = zeros(length(ms1.par))
    G3 = zeros(length(ms1.par))
    SurvivalAnalysis.score!(ms1, G1, ms1.par)
    SurvivalAnalysis.score!(ms2, G2, ms1.par)
    SurvivalAnalysis.score!(ms3, G3, ms1.par)
    @assert isapprox(G1, G2)
    @assert isapprox(G1, G3)

    @assert isapprox(ms1.support, ms2.support)
    @assert isapprox(ms1.support, ms3.support)
    @assert isapprox(ms1.duration, ms2.duration)
    @assert isapprox(ms1.duration, ms3.duration)

    haz1 = clamp.(ms1.hazard, 0, 2)
    haz2 = clamp.(ms2.hazard, 0, 2)
    haz3 = clamp.(ms3.hazard, 0, 2)
    @test mean(abs, haz1 - haz2) < 0.01
    @test mean(abs, haz1 - haz3) < 0.01
end
