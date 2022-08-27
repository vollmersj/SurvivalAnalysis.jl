#-------------------
# ParametricSurvival
#-------------------
"""
    ParametricSurvival <: SurvivalModel

Abstract type for all fully-parametric survival models implemented in, or extending, this
package. Type 'inherits' [`SurvivalModel`](@ref).
from JuliaStats.StatisticalModel to enable formula fitting and predicting interface.
"""
abstract type ParametricSurvival <: SurvivalModel end

StatsBase.coef(obj::ParametricSurvival) = obj.coefficients

function StatsBase.fit(t::Type{<:ParametricSurvival}, X::AbstractMatrix{<:Real}, Y::RCSurv,
    d::Type{T}; init::Number = 1) where {T <: ContinuousUnivariateDistribution}
    @assert d in [Weibull, Exponential]
    # crude method to force intercept
    X = X[:,1] == ones(size(X, 1)) ? X : hcat(ones(size(X, 1)), X)
    return fit!(t(d, X), X, Y, init)
end

function _fitParametricSurvival(obj, X, Y, init, llik, dtrafo)
    # θ[1] = scale, θ[2:end] = βs
    nβ = size(X, 2)
    npar = nβ + 1
    X = Matrix(X)

    init = [init, zeros(nβ)...]

    func = TwiceDifferentiable(
        θ -> θ[1] <= 0 ? Inf : -∑(llik(obj.baseline, X, Y.time, Y.status, θ[1], θ[2:end])),
        init; autodiff=:forward);

    opt = optimize(func, init)
    θ̂ = Optim.minimizer(opt)
    obj.routine = opt

    obj.hessian = hessian!(func, θ̂)
    ## FIXME #27 - NEED TO FIGURE OUT WHY INVERSION SOMETIMES FAILS
    obj.var_cov = try
        inv(obj.hessian)
    catch
        fill(NaN, (npar, npar))
    end
    obj.scale = θ̂[1]

    obj.coefficients = θ̂[2:end]

    ## FIXME #29 - NEED TO FIGURE OUT WHY VAR_COV SOMETIMES < 0
    obj.tstats = try
        obj.coefficients./sqrt.(diag(obj.var_cov)[2:end])
    catch
        fill(NaN, nβ)
    end

    obj.baseline = dtrafo(obj.baseline, obj.coefficients, obj.scale)

    return obj
end


baseline(mm::StatsModels.TableStatisticalModel{<:ParametricSurvival, <:AbstractMatrix}) =
    mm.model.baseline
baseline(mm::ParametricSurvival) = mm.baseline

Distributions.scale(
    mm::StatsModels.TableStatisticalModel{<:ParametricSurvival, <:AbstractMatrix}) =
    mm.model.scale
Distributions.scale(mm::ParametricSurvival) = mm.scale

function Base.show(io::IO, mm::ParametricSurvival)
    println(io, typeof(mm))
    println(io)
    println(io, "Distr:")
    println(io, baseline(mm))
    println(io)
    println(io,"Coefficients:")

    header = ["(Scale)", ["x$i" for i = 0:(length(coef(mm))-1)]...]
    header[2] = "(Intercept)"
    data = [mm.scale coef(mm)...]
    pretty_table(io, data,  header = header, vlines = :none, hlines = :none)
end

function Base.show(io::IO, mm::StatsModels.TableStatisticalModel{<:ParametricSurvival})
    println(io, typeof(mm))
    println(io)
    println(io, mm.mf.f)
    println(io)
    println(io, "Distr:")
    println(io, baseline(mm))
    println(io)
    println(io, "Coefficients:")

    header = ["(Scale)", coefnames(mm.mf)...]
    data = [scale(mm) coef(mm)...]
    pretty_table(io, data,  header = header, vlines = :none, hlines = :none)
end

#-------------------
# ParametricPH
#-------------------

mutable struct ParametricPH{T} <: ParametricSurvival where
    {T <: ContinuousUnivariateDistribution}
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
    return ParametricPH(
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

function StatsBase.fit!(obj::ParametricPH, X::AbstractMatrix{<:Real}, Y::RCSurv,
                        init::Number)

    function llik(ζ, x, t, δ, ϕ, β)
        if ζ === Exponential
            return (δ .* (x*β)) .- (exp.(x*β) .* t)
        else # weibull
            return (δ .* (log(1/ϕ) .+ (((1/ϕ)-1) .* log.(t)) .+ x*β)) .-
                (exp.(x*β) .* t.^(1/ϕ))
        end
    end

    function dtrafo(ζ, β, ϕ)
        if ζ === Exponential
            return Exponential(exp(-β[1]))
        else # weibull
            return Weibull(1/ϕ, exp(-β[1] * ϕ))
        end
    end

    return _fitParametricSurvival(
        obj,
        X,
        Y,
        init,
        llik,
        dtrafo
    )
end

function StatsBase.predict(fit::ParametricPH, X::AbstractMatrix{<:Real})
    lp = X * fit.coefficients[2:end] # intercept removed
    distr = ContinuousPHDistribution.(fit.baseline, lp)
    return SurvivalPrediction(distr = distr, lp = lp, crank = lp)
end

#-------------------
# ParametricAFT
#-------------------

mutable struct ParametricAFT{T} <: ParametricSurvival where
    {T <: ContinuousUnivariateDistribution}
    coefficients::Vector{Float64}
    scale::Float64
    hessian::Matrix
    var_cov::Matrix
    tstats::Vector
    baseline::Union{T, Type{T}}
    routine
end

function ParametricAFT(d::Type{T}, X) where {T <: ContinuousUnivariateDistribution}
    nβ = size(X, 2)
    return ParametricAFT(
        zeros(nβ),
        0.0,
        zeros(nβ, nβ),
        zeros(nβ, nβ),
        zeros(nβ),
        d,
        missing
    )
end

aft(args...; kwargs...) = fit(ParametricAFT, args...; kwargs...)

function StatsBase.fit!(obj::ParametricAFT, X::AbstractMatrix{<:Real}, Y::RCSurv,
                        init::Number)

    function llik(ζ, x, t, δ, ϕ, β)
        if ζ === Exponential
            return (δ .* (-x*β)) .- (exp.(-x*β) .* t)
        else # weibull
            return (δ .* (log(1/ϕ) .+ (((1/ϕ)-1) .* log.(t)) .- x*β ./ ϕ)) .-
                (exp.((-x*β) ./ ϕ) .* t.^(1/ϕ))
        end
    end

    function dtrafo(ζ, β, ϕ)
        if ζ === Exponential
            return Exponential(exp(β[1]))
        else # weibull
            return Weibull(1/ϕ, exp(β[1]))
        end
    end

    return _fitParametricSurvival(
        obj,
        X,
        Y,
        init,
        llik,
        dtrafo
    )
end

function StatsBase.predict(fit::ParametricAFT, X::AbstractMatrix{<:Real})
    lp = X * fit.coefficients[2:end] # intercept removed
    distr = ContinuousAFTDistribution.(fit.baseline, lp)
    return SurvivalPrediction(distr = distr, lp = lp, crank = -lp)
end
