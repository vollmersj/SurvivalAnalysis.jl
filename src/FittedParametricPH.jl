mutable struct ParametricPH{T} <: FittedParametric where {T <: ContinuousUnivariateDistribution}
    coefficients::Vector{Float64}
    scale::Float64
    hessian::Matrix
    var_cov::Matrix
    tstats::Vector
    distribution::Union{T, Type{T}}
    routine
end

function ParametricPH(d::Type{T}, X) where {T <: ContinuousUnivariateDistribution}
    nβ = size(X, 2) + 1
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

function StatsBase.fit(::Type{ParametricPH}, X::AbstractMatrix{<:Real}, Y::Survival.rcSurv,
                d::Type{T}, init::Number = 1) where {T <: ContinuousUnivariateDistribution}
    @assert d in [Weibull, Exponential]
    fit!(ParametricPH(d, X), X, Y, init)
end

function StatsBase.fit!(obj::ParametricPH, X::AbstractMatrix{<:Real}, Y::Survival.rcSurv,
                        init::Number)
    # θ[1] = scale, θ[2] = β₀, θ[...] = β...
    nβ = size(X, 2)
    npar = nβ + 2
    X = Matrix(X)

    init = [init, 0, zeros(nβ)...]

    function loglik(x, t, δ, ϕ, β₀, β)
        ϕ <= 0 && return Inf # reject impossible candidates
        if obj.distribution == Exponential
            l = (δ .* (β₀ .+ x*β)) .- (exp.(β₀ .+ x*β) .* t)
        else
            l = (δ .* (log(1/ϕ) .+ (((1/ϕ)-1) .* log.(t)) .+ β₀ .+ x*β)) .- (exp.(β₀ .+ x*β) .* t.^(1/ϕ))
        end
        -Survival.∑(l)
    end

    func = TwiceDifferentiable(
        θ -> loglik(X, Y.time, Y.status, θ[1], θ[2], θ[3:length(θ)]),
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

    obj.distribution = obj.distribution == Exponential ?
        Exponential(exp(-β̂[1])) :
        Weibull(1/obj.scale, exp(-β̂[1] * obj.scale))

    obj
end

function StatsBase.predict(fit::ParametricPH, X::AbstractMatrix{<:Real})
    η = X * fit.coefficients.betas[2:end]
    ζ = ParametricPH.(fit.baseline, η)
    _survPredict(ζ = ζ, η = η, ϕ = η)
end

coef(obj::FittedParametric) = obj.coefficients
scale(object::FittedParametric) = scale(obj.distribution)
shape(object::FittedParametric) = shape(obj.distribution)
