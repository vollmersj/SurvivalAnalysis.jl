mutable struct ParametricPH{T} <: FittedParametric where {T <: ContinuousUnivariateDistribution}
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

ph(X, y, args...; kwargs...) = fit(ParametricPH, X, y, args...; kwargs...)

function StatsBase.fit(::Type{ParametricPH}, X::AbstractMatrix{<:Real}, Y::rcSurv,
                d::Type{T}, init::Number = 1) where {T <: ContinuousUnivariateDistribution}
    @assert d in [Weibull, Exponential]
    fit!(ParametricPH(d, X), X, Y, init)
end

function StatsBase.fit!(obj::ParametricPH, X::AbstractMatrix{<:Real}, Y::rcSurv,
                        init::Number)
    # θ[1] = scale, θ[2:end] = βs
    nβ = size(X, 2)
    npar = nβ + 1
    X = Matrix(X)

    init = [init, zeros(nβ)...]

    function loglik(x, t, δ, ϕ, β)
        ϕ <= 0 && return Inf # reject impossible candidates
        if obj.baseline == Exponential
            l = (δ .* (x*β)) .- (exp.(x*β) .* t)
        else
            l = (δ .* (log(1/ϕ) .+ (((1/ϕ)-1) .* log.(t)) .+ x*β)) .- (exp.(x*β) .* t.^(1/ϕ))
        end
        -∑(l)
    end

    func = TwiceDifferentiable(
        θ -> loglik(X, Y.time, Y.status, θ[1], θ[2:end]),
        init; autodiff=:forward);

    opt = optimize(func, init)
    θ̂ = Optim.minimizer(opt)
    obj.routine = opt

    nθ = length(θ̂)
    obj.hessian = hessian!(func, θ̂)
    obj.var_cov = try
        inv(obj.hessian)
    catch
        fill(NaN, (npar, npar))
    end
    obj.scale = θ̂[1]

    β̂ = θ̂[2:nθ]
    obj.tstats = β̂./sqrt.(diag(obj.var_cov)[2:nθ])
    obj.coefficients = β̂

    obj.baseline = obj.baseline == Exponential ?
        Exponential(exp(-β̂[1])) :
        Weibull(1/obj.scale, exp(-β̂[1] * obj.scale))

    obj
end

function StatsBase.predict(fit::ParametricPH, X::AbstractMatrix{<:Real})
    η = X * fit.coefficients
    ζ = ContinuousPHDistribution.(fit.baseline, η)
    SurvivalPrediction(ζ = ζ, η = η, ϕ = η)
end
