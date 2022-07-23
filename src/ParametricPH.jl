mutable struct ParametricPH{T} <: ParametricSurvival where {T <: ContinuousUnivariateDistribution}
    coefficients::Vector{Float64}
    scale::Float64
    hessian::Matrix
    var_cov::Matrix
    tstats::Vector
    baseline::Union{T, Type{T}}
    routine
end

function ParametricPH(d::Type{T}, X) where {T <: ContinuousUnivariateDistribution}
    nβ = size(X, 2)
    ParametricPH(
        zeros(nβ),
        0.0,
        zeros(nβ, nβ),
        zeros(nβ, nβ),
        zeros(nβ),
        d,
        missing
    )
end

ph(args...; kwargs...) = fit(ParametricPH, args...; kwargs...)

function StatsBase.fit!(obj::ParametricPH, X::AbstractMatrix{<:Real}, Y::rcSurv,
                        init::Number)

    function llik(ζ, x, t, δ, ϕ, β)
        if ζ == Exponential
            l = (δ .* (x*β)) .- (exp.(x*β) .* t)
        else # weibull
            l = (δ .* (log(1/ϕ) .+ (((1/ϕ)-1) .* log.(t)) .+ x*β)) .- (exp.(x*β) .* t.^(1/ϕ))
        end
    end

    function dtrafo(ζ, β, ϕ)
        if ζ == Exponential
            Exponential(exp(-β[1]))
        else # weibull
            Weibull(1/ϕ, exp(-β[1] * ϕ))
        end
    end

    _fitParametricSurvival(
        obj,
        X,
        Y,
        init,
        llik,
        dtrafo
    )
end

function StatsBase.predict(fit::ParametricPH, X::AbstractMatrix{<:Real})
    η = X * fit.coefficients[2:end] # intercept removed
    ζ = ContinuousPHDistribution.(fit.baseline, η)
    SurvivalPrediction(ζ = ζ, η = η, ϕ = η)
end
