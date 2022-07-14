struct NelsonAalen <: NonParametricEstimator
    time::Vector{Float64}
    survival::Vector{Float64}
    sd::Vector{Float64}
    d::DiscreteNonParametric
end

function nelson(Surv::rcSurv)
    fit_NPE(
        NelsonAalen,
        Surv,
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
    confint_NPE(
        na, t, α,
        (E, q, w) -> map(x -> min(1, max(0, exp(-x))), -log(E.survival[w]) ∓ (q * E.sd[w]))
    )
end
