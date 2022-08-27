#-------------------
# ParametricSurvival
#-------------------
"""
    ParametricSurvival <: SurvivalModel

Abstract type for all fully-parametric survival models implemented in, or extending, this
package. Type 'inherits' [`SurvivalModel`](@ref).

Available methods:

* [`coef`](@ref) - Extract fitted coefficients
* [`fit`](@ref) and [`predict`](@predict) - Fit model and make predictions from fitted model with `@formula` or `matrix` interface, see [Fitting and predicting](@ref)
* [`baseline`](@ref) - Extract fitted baseline distribution, see [ph](@ref) and [aft](@ref) for more
* [`scale`](@ref) - Extract scale parameter of fitted distribution

All distributions are fitted with the Kalbfleisch-Prentice parametrisation and then
converted as required to make use of `Distributions.jl`

Objects inheriting from this should have the following fields:

* `coefficients::Vector{Float64}` - Fitted coefficients
* `scale::Float64` - Fitted scale parameter for baseline distribution *before* transformation
* `hessian::Matrix` - Hessian from `Optim`
* `var_cov::Matrix` - Covariance matrix
* `tstats::Vector` - t-statistics
* `baseline<:ContinuousUnivariateDistribution` - Fitted baseline distribution
* `routine` - Optimisation routine from `Optim`
"""
abstract type ParametricSurvival <: SurvivalModel end

"""
    coef(mm::StatsModels.TableStatisticalModel{<:ParametricSurvival, <:AbstractMatrix})
    coef(mm::ParametricSurvival)

Extract coefficients from fitted [`ParametricSurvival`](@ref) object.
"""
StatsBase.coef(
    mm::StatsModels.TableStatisticalModel{<:ParametricSurvival, <:AbstractMatrix}) =
    mm.model.coefficients
StatsBase.coef(obj::ParametricSurvival) = obj.coefficients

"""
    fit(t::Type{<:ParametricSurvival}, X::AbstractMatrix{<:Real}, Y::RCSurv,
        d::Type{T}; init::Number = 1) where {T <: ContinuousUnivariateDistribution}

Fit a [`ParametricSurvival`](@ref) survival model using matrix interface. It is recommended
to use [`ph`](@ref) or [`aft`](@ref) directly instead.
"""
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

"""
    baseline(mm::StatsModels.TableStatisticalModel{<:ParametricSurvival, <:AbstractMatrix})
    baseline(mm::ParametricSurvival)

Extract baseline distribution from fitted [`ParametricSurvival`](@ref) object.
"""
baseline(mm::StatsModels.TableStatisticalModel{<:ParametricSurvival, <:AbstractMatrix}) =
    mm.model.baseline
baseline(mm::ParametricSurvival) = mm.baseline

"""
    scale(mm::StatsModels.TableStatisticalModel{<:ParametricSurvival, <:AbstractMatrix})
    scale(mm::ParametricSurvival)

Extract estimated scale parameter (*before* transformation) for baseline distribution from
fitted [`ParametricSurvival`](@ref) object.
"""
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
"""
    ParametricPH{<:ContinuousUnivariateDistribution} <: ParametricSurvival

See [`ParametricSurvival`](@ref).
"""
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

"""
    ph(f::@formula, data::DataFrames.DataFrame, d::Type{T}; init::Number = 1)
    ph(X::AbstractMatrix{<:Real}, Y::RCSurv, d::Type{T}; init::Number = 1)
    fit(ParametricPH, X::AbstractMatrix{<:Real}, Y::RCSurv, d::Type{T}; init::Number = 1)

Fit a fully-parametric proportional hazards (PH) model with baseline distribution `d`.
See [Fitting and predicting](@ref) and examples below for fitting interfaces, we recommend
using `ph(::@formula...)`.

Fully-parametric PH models are defined by

```math
h(t) = h₀(t)exp(Xβ)
```

where ``β`` are coefficients to be estimated, ``X`` are covariates, and `h₀` is the hazard
function of an assumed baseline distribution, `d`. Available choices for distributions are:

* Distributions.Exponential
* Distributions.Weibull

No other distributions have the PH assumption, which assumes that the risk of an event
taking place is constant over time.

Future additions:

* Methods for testing if the assumption is valid for your data (see [#11](@ref)).
* Methods to calculate hazards ratios and add to `show` (see [#40](@ref)).

Function returns a [`ParametricPH`](@ref) struct.

# Examples
```jldoctest
julia> Y = [1,1,4,6,8,4,9,4,5,10];

julia> D = [true, false, false, false, true, false, false, true, true, false];

julia> X = [1,9,3,4,20,-4,pi,exp(5),log(8),0];

julia> data = DataFrame(Y = Y, D = D, X = X);

julia> f = ph(@formula(Srv(Y, D) ~ X), data, Exponential)
StatsModels.TableStatisticalModel{ParametricPH{Exponential}, Matrix{Float64}}

(Y,D;+) ~ 1 + X

Distr:
Exponential{Float64}(θ=17.121272620254064)

Coefficients:
 (Scale)  (Intercept)          X
     1.0     -2.84032  0.0101247
```
"""
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

"""
    predict(fit::ParametricPH, X::AbstractMatrix{<:Real})
    predict(fit::StatsModels.TableStatisticalModel{ParametricPH, Matrix{Float64}},
        data::DataFrames.DataFrame)

Make predictions from a fitted [`ParametricPH`](@ref) model. See
[Fitting and predicting](@ref) and examples below for predicting interfaces, we recommend
using `predict(fit, data::DataFrame)`.

Three prediction types can be made from a fitted PH model:

* `lp` - ``Xβ̂``
* `crank` - ``Xβ̂``
* `distr` - ``F(t) = 1 - Ŝ₀(t)^{exp(Xβ̂)}``

where ``β̂`` are estimated coefficients, ``X`` are covariates from the new data, and ``Ŝ₀`` is the estimated baseline distribution survival function. Predicted distributions are returned as `ContinuousPHDistribution <: Distributions.ContinuousUnivariateDistribution`.

Note❗The PH model assumes that a higher linear predictor means a higher risk of event and therefore a lower survival time, i.e., ``βXᵢ > βXⱼ → hᵢ(t) > hⱼ(t)` - hence `crank = lp`. This means when calculating [`concordance`](@ref) you *must* include `rev = true`.

Future updates will add transformation methods for more prediction types (see [#12](@ref)).

Function returns a [`SurvivalPrediction`](@ref) struct.

# Examples
```jldoctest
julia> Y = [1,1,4,6,8,4,9,4,5,10];

julia> D = [true, false, false, false, true, false, false, true, true, false];

julia> X = [1,9,3,4,20,-4,pi,exp(5),log(8),0];

julia> data = DataFrame(Y = Y, D = D, X = X);

julia> f = ph(@formula(Srv(Y, D) ~ X), data, Weibull);

julia> predict(f, DataFrame(X = [2,9,1,1,exp(8),log(9),2^3,5^exp(2),1]));
```
"""
function StatsBase.predict(fit::ParametricPH, X::AbstractMatrix{<:Real})
    lp = X * fit.coefficients[2:end] # intercept removed
    distr = ContinuousPHDistribution.(fit.baseline, lp)
    return SurvivalPrediction(distr = distr, lp = lp, crank = lp)
end

#-------------------
# ParametricAFT
#-------------------
"""
    ParametricAFT{<:ContinuousUnivariateDistribution} <: ParametricSurvival

See [`ParametricSurvival`](@ref).
"""
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

"""
    aft(f::@formula, data::DataFrames.DataFrame, d::Type{T}; init::Number = 1)
    aft(X::AbstractMatrix{<:Real}, Y::RCSurv, d::Type{T}; init::Number = 1)
    fit(ParametricAFT, X::AbstractMatrix{<:Real}, Y::RCSurv, d::Type{T}; init::Number = 1)

Fit a fully-parametric accelerated failure time (AFT) model with baseline distribution `d`. See [Fitting and predicting](@ref) and examples below for fitting interfaces, we recommend using `aft(::@formula...)`.

Fully-parametric AFT models are defined by

```math
h(t) = e^{-Xβ} h₀(t e^{-Xβ})
```

where ``β`` are coefficients to be estimated, ``X`` are covariates, and `h₀` is the hazard function of an assumed baseline distribution, `d`. Available choices for distributions are:

* Distributions.Exponential
* Distributions.Weibull

AFT models assumes that an increase in a covariate results in an acceleration of the event by a constant. This is best explained by example. The above formula can also be expressed as ``S(t) = S₀(exp(-η)t)`` then let ``ηᵢ = log(2)`` and ``ηⱼ = log(1)`` so ``ηᵢ > ηⱼ``. Then Sᵢ(t) = S₀(0.5t) and Sⱼ(t) = S₀(t) and so for all ``t``, ``Sᵢ(t) ≥ Sⱼ(t)`` as S is a decreasing function.

Future additions:

* More baseline distributions (see [#15](@ref))

Function returns a [`ParametricAFT`](@ref) struct.

# Examples
```jldoctest
julia> Y = [1,1,4,6,8,4,9,4,5,10];

julia> D = [true, false, false, false, true, false, false, true, true, false];

julia> X = [1,9,3,4,20,-4,pi,exp(5),log(8),0];

julia> data = DataFrame(Y = Y, D = D, X = X);

julia> f = aft(@formula(Srv(Y, D) ~ X), data, Weibull)
StatsModels.TableStatisticalModel{ParametricAFT{Weibull}, Matrix{Float64}}

(Y,D;+) ~ 1 + X

Distr:
Weibull{Float64}(α=1.7439548025959823, θ=11.754846456706442)

Coefficients:
  (Scale)  (Intercept)            X
 0.573409      2.46427  -0.00741527
```
"""
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

"""
    predict(fit::ParametricAFT, X::AbstractMatrix{<:Real})
    predict(fit::StatsModels.TableStatisticalModel{ParametricAFT, Matrix{Float64}},
        data::DataFrames.DataFrame)

Make predictions from a fitted [`ParametricAFT`](@ref) model. See [Fitting and predicting](@ref) and examples below for predicting interfaces, we recommend using `predict(fit, data::DataFrame)`.

Three prediction types can be made from a fitted AFT model:

* `lp` - ``Xβ̂``
* `crank` - ``-Xβ̂``
* `distr` - ``F(t) = F̂₀(t/exp(Xβ̂))``

where ``β̂`` are estimated coefficients, ``X`` are covariates from the new data, and ``F̂₀`` is the estimated baseline distribution CDF function. Predicted distributions are returned as `ContinuousAFTDistribution <: Distributions.ContinuousUnivariateDistribution`.

Note❗The AFT model assumes that a higher linear predictor means a lower risk of event and therefore a higher survival time, i.e., ``βXᵢ > βXⱼ → hᵢ(t) < hⱼ(t)`` - hence `crank = -lp`.

Future updates will add transformation methods for more prediction types (see [#12](@ref)).

Function returns a [`SurvivalPrediction`](@ref) struct.

# Examples
```jldoctest
julia> Y = [1,1,4,6,8,4,9,4,5,10];

julia> D = [true, false, false, false, true, false, false, true, true, false];

julia> X = [1,9,3,4,20,-4,pi,exp(5),log(8),0];

julia> data = DataFrame(Y = Y, D = D, X = X);

julia> f = aft(@formula(Srv(Y, D) ~ X), data, Exponential);

julia> predict(f, DataFrame(X = [2,9,1,1,exp(8),log(9),2^3,5^exp(2),1]));
```
"""
function StatsBase.predict(fit::ParametricAFT, X::AbstractMatrix{<:Real})
    lp = X * fit.coefficients[2:end] # intercept removed
    distr = ContinuousAFTDistribution.(fit.baseline, lp)
    return SurvivalPrediction(distr = distr, lp = lp, crank = -lp)
end
