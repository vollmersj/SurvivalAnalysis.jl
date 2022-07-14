struct KaplanMeier <: NonParametricEstimator
    time::Vector{Float64}
    survival::Vector{Float64}
    sd::Vector{Float64}
    d::DiscreteNonParametric
end

function kaplan(Surv::rcSurv)
    fit_NPE(
        KaplanMeier,
        Surv,
        (d, n) -> 1 - (d / n),
        (d, n) -> d / (n * (n - d)),
        # makes use of KM being a plug-in estimator
        cumprod,
        # Calculates standard error using the method of Kalbfleisch and Prentice (1980)
        (v, p) -> map((qₜ, vₜ) -> √(vₜ / (log(qₜ)^2)), p, v)
    )
end

function StatsBase.confint(km::KaplanMeier, t::Number; α::Float64 = 0.05)
    confint_NPE(
        km, t, α,
        (E, q, w) -> map(x -> exp(-exp(x)), log(-log(E.survival[w])) ∓ (q * E.sd[w]))
    )
end
