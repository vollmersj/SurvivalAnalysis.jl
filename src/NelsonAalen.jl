mutable struct NelsonAalen <: SurvivalAnalysis.NonParametricEstimator
    time::Vector{Float64}
    survival::Vector{Float64}
    sd::Vector{Float64}
    distribution::DiscreteNonParametric

    NelsonAalen() = new()
end

nelson_aalen(args...; kwargs...) = StatsBase.fit(KaplanMeier, args...; kwargs...)

function StatsBase.fit!(obj::NelsonAalen, Y::RCSurv)
    return _fit_npe(
        obj,
        Y,
        (d, n) -> d / n,
        (d, n) -> (d * (n - d)) / (n^3),
        # makes use of NA being a plug-in estimator then converts to surv
        (p) -> exp.(-cumsum(p)),
        (v, p) -> .√v
    )
end

# Calculates pointwise confidence intervals for *survival predictions*
# Unclear if this even needs to be a different method, could just use the confint around
# survival for both NPEs
function StatsBase.confint(na::NelsonAalen, t::Number; α::Float64 = 0.05)
    return _confint_npe(
        na, t, α,
        (E, q, w) -> map(x -> min(1, max(0, exp(-x))), -log(E.survival[w]) ∓ (q * E.sd[w]))
    )
end
