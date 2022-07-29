mutable struct KaplanMeier <: NonParametricEstimator
    time::Vector{Float64}
    survival::Vector{Float64}
    std::Vector{Float64}
    distribution::DiscreteNonParametric

    KaplanMeier() = new()
end

kaplan_meier(args...; kwargs...) = StatsBase.fit(KaplanMeier, args...; kwargs...)
kaplan_meier(Y::RCSurv) = StatsBase.fit(KaplanMeier, Y)

function StatsBase.fit!(obj::KaplanMeier, Y::RCSurv)
    return _fit_npe(
        obj,
        Y,
        (d, n) -> 1 - (d / n),
        (d, n) -> d / (n * (n - d)),
        cumprod, # makes use of KM being a plug-in estimator
        # Calculates standard error using the method of Kalbfleisch and Prentice (1980)
        (v, p) -> map((qₜ, vₜ) -> √(vₜ / (log(qₜ)^2)), p, v)
    )
end

function StatsBase.confint(km::KaplanMeier, t::Number; α::Float64 = 0.05)
    return _confint_npe(
        km, t, α,
        (E, q, w) -> map(x -> exp(-exp(x)), log(-log(E.survival[w])) ∓ (q * E.std[w]))
    )
end
