#-------------------
# SurvivalEstimator
#-------------------
"""
    SurvivalEstimator <: StatisticalModel

Abstract type for all non-parametric estimators implemented in, or extending, this package. Type 'inherits' from [StatsAPI.StatisticalModel](https://github.com/JuliaStats/StatsAPI.jl/blob/main/src/statisticalmodel.jl) to enable formula fitting and predicting interface. Note❗This may be abstracted further into `ConditionalSurvivalEstimator <: SurvivalEstimator` and `UnconditionalSurvivalEstimator <: SurvivalEstimator`.

Available methods:

* [`fit`](@ref) and [`predict`](@ref) - Fit model and make predictions from fitted model with `@formula` or `matrix` interface, see [Fitting and predicting](@ref)
* [`confint`](@ref) - Calculate confidence intervals around estimates
* [`time`](@ref) - Extract fitted times
* [`survival`](@ref) - Extract estimated survival probabilities
* [`std`](@ref) - Extract computed standard deviation
* [`distr`](@ref) - Extract fitted survival distribution

Objects inheriting from this should have the following fields:

* `time::Vector{Float64}` - Fitted survival times
* `survival::Vector{Float64}` - Estimated survival probabilities
* `std::Vector{Float64}` - Computed standard deviation
* `distr::DiscreteNonParametric` - Fitted survival distribution
* `stats::NamedTuple` - Summary statistics such as numbers at risk, dead, censored
"""
abstract type SurvivalEstimator <: StatisticalModel end

"""
    fit(obj::Type{<:SurvivalEstimator}, X::AbstractMatrix{<:Real}, Y::RCSurv)
    fit(obj::Type{<:SurvivalEstimator}, Y::RCSurv)

Fit a [`SurvivalEstimator`](@ref) survival model using matrix interface. It is recommended
to use [`kaplan_meier`](@ref) or [`nelson_aalen`](@ref) directly instead.
"""
function StatsBase.fit(
    obj::Type{<:SurvivalEstimator}, X::AbstractMatrix{<:Real}, Y::RCSurv)
    # TODO - need to add something here to check if rhs is intercept only or stratified
    #  currently stratified not supported
    obj = obj()
    obj.stats = (n = Y.stats.nrisk[1], ncens = sum(Y.stats.ncens),
        nevents = sum(Y.stats.nevents))
    return fit!(obj, Y)
end

StatsBase.fit(obj::Type{<:SurvivalEstimator}, Y::RCSurv) = fit!(obj(), Y)

"""
    predict(fit::SurvivalEstimator, X::AbstractMatrix{<:Real})
    predict(fit::StatsModels.TableStatisticalModel{SurvivalEstimator, Matrix{Float64}},
        data::DataFrames.DataFrame)

Make predictions from a fitted [`SurvivalEstimator`](@ref) model. See [Fitting and predicting](@ref) and examples below for predicting interfaces, we recommend using `predict(fit, data::DataFrame)`.

Two prediction types can be made from a non-parametric estimator:

* `distr` - See [`kaplan_meier`](@ref) or [`nelson_aalen`](@ref) for formulae
* `survival_matrix` - Matrix of survival probabilities

Predicted distributions are returned as `Distributions.DiscreteNonParametric`.

Function returns a [`SurvivalPrediction`](@ref) struct.

# Examples
```jldoctest
julia> data = DataFrame(Y = [1,1,4,6,8,4,9,4,5,10], D = [true, false, false, false, true, false, false, true, true, false]);

julia> f = km(@formula(Srv(Y, D) ~ 1), data);

julia> predict(f, DataFrame(X = [2,9,1,1,exp(8),log(9),2^3,5^exp(2),1]));
```
"""
StatsBase.predict(fit::SurvivalEstimator, data::DataFrame; kwargs...) =
    predict(fit, Matrix(data); kwargs...)

function StatsBase.predict(fit::SurvivalEstimator, X::AbstractMatrix{<:Real})
    # TODO - need to add something here to check if rhs is intercept only or stratified
    #  currently stratified not supported
    n = size(X, 1)
    return SurvivalPrediction(
        distr = fill(fit.distr, n),
        fit_times = fit.time,
        survival_matrix = repeat(fit.survival', n)
    )
end

function StatsModels.predict(
    mm::StatsModels.TableStatisticalModel{<:SurvivalEstimator, <:AbstractMatrix}, data;
    kwargs...
)
    Tables.istable(data) ||
        throw(ArgumentError("expected data in a Table, got $(typeof(data))"))

    f = mm.mf.f
    if f.rhs isa MatrixTerm{Tuple{InterceptTerm{true}}}
        new_x = reshape(data[!,1], :, 1) # pick an arbitrary column - only care about size
    else
        cols, nonmissings = StatsModels.missing_omit(StatsModels.columntable(data), f.rhs)
        new_x = modelcols(f.rhs, cols)
        new_x = reshape(new_x, size(new_x, 1), :)[:,2:end] # remove intercept
    end
    return predict(mm.model, new_x; kwargs...)
end

function _fit_npe(obj::SurvivalEstimator, Surv::RCSurv, point_est::Function,
    var_est::Function, surv_trafo::Function, std_trafo::Function)
    stats = surv_stats(Surv, events_only = true)
    n = length(stats.time)
    p = zeros(n)
    v = zeros(n)
    for (i, t) in enumerate(stats.time)
        p[i] = point_est(stats.nevents[i], stats.nrisk[i])
        v[i] = var_est(stats.nevents[i], stats.nrisk[i])
    end

    obj.survival = surv_trafo(p)
    obj.time = stats.time
    # for both NPEs variance calculated as plug-in
    obj.std = std_trafo(cumsum(v), obj.survival)
    # calculate pmf and create distr - Set S(0) = 1
    pₓ = [0, abs.(diff(obj.survival))...]
    obj.distr = DiscreteNonParametric([0, stats.time...],  pₓ, check_args = false)

    return obj
end

# Calculates confidence intervals using the method of Kalbfleisch and Prentice (1980)
function _confint_npe(npe, t, level, trafo)
    q = quantile(Normal(), 1 - (1 - level)/2)
    which = searchsortedlast(npe.time, t)
    return trafo(npe, q, which)
end

# syntactic sugar for vectorising over all times
function StatsBase.confint(npe::SurvivalEstimator; level::Float64 = 0.95)
    return confint.(Ref(npe), npe.time; level = level)
end

function StatsBase.confint(
    npe::StatsModels.TableStatisticalModel{<:SurvivalEstimator, <:AbstractMatrix};
    level::Float64 = 0.95
)
    return confint.(Ref(npe.model), npe.model.time; level = level)
end

"""
    time(npe::SurvivalEstimator)
    time(npe::StatsModels.TableStatisticalModel{<:SurvivalEstimator, <:AbstractMatrix})

Return fitted times from a [`SurvivalEstimator`](@ref).

# Examples
```jldoctest
julia> data = DataFrame(Y = [1,1,4,6,8,4,9,4,5,10], D = [true, false, false, false, true, false, false, true, true, false]);

julia> time(km(@formula(Srv(Y, D) ~ 1), data))
4-element Vector{Float64}:
 1.0
 4.0
 5.0
 8.0
```
"""
Base.time(npe::SurvivalEstimator) = npe.time

"""
    survival(npe::SurvivalEstimator)
    survival(npe::StatsModels.TableStatisticalModel{<:SurvivalEstimator, <:AbstractMatrix})

Return fitted survival probabilities from a [`SurvivalEstimator`](@ref).

# Examples
```jldoctest
julia> data = DataFrame(Y = [1,1,4,6,8,4,9,4,5,10], D = [true, false, false, false, true, false, false, true, true, false]);

julia> survival(km(@formula(Srv(Y, D) ~ 1), data))
4-element Vector{Float64}:
 0.9
 0.7875
 0.63
 0.42000000000000004
```
"""
survival(npe::SurvivalEstimator) = npe.survival

"""
    std(npe::SurvivalEstimator)
    std(npe::StatsModels.TableStatisticalModel{<:SurvivalEstimator, <:AbstractMatrix})

Return computed standard deviation from a [`SurvivalEstimator`](@ref).

# Examples
```jldoctest
julia> data = DataFrame(Y = [1,1,4,6,8,4,9,4,5,10], D = [true, false, false, false, true, false, false, true, true, false]);

julia> std(km(@formula(Srv(Y, D) ~ 1), data))
4-element Vector{Float64}:
 1.0004625991132958
 0.712458742539604
 0.6082063644330836
 0.5713145523891875
```
"""
StatsBase.std(npe::SurvivalEstimator) = npe.std

"""
    distr(npe::SurvivalEstimator)
    distr(npe::StatsModels.TableStatisticalModel{<:SurvivalEstimator, <:AbstractMatrix})

Return estimated survival distribution from a [`SurvivalEstimator`](@ref).

# Examples
```jldoctest
julia> data = DataFrame(Y = [1,1,4,6,8,4,9,4,5,10], D = [true, false, false, false, true, false, false, true, true, false]);

julia> distr(km(@formula(Srv(Y, D) ~ 1), data))
DiscreteNonParametric{Float64, Float64, Vector{Float64}, Vector{Float64}}(
support: [0.0, 1.0, 4.0, 5.0, 8.0]
p: [0.0, 0.11250000000000004, 0.15749999999999997, 0.20999999999999996]
)
```
"""
distr(npe::SurvivalEstimator) = npe.distr

function Base.time(
    npe::StatsModels.TableStatisticalModel{<:SurvivalEstimator, <:AbstractMatrix}
)
    return npe.model.time
end
function survival(
    npe::StatsModels.TableStatisticalModel{<:SurvivalEstimator, <:AbstractMatrix}
)
    return npe.model.survival
end
function StatsBase.std(
    npe::StatsModels.TableStatisticalModel{<:SurvivalEstimator, <:AbstractMatrix}
)
    return npe.model.std
end
function distr(
    npe::StatsModels.TableStatisticalModel{<:SurvivalEstimator, <:AbstractMatrix}
)
    return npe.model.distr
end

function Base.show(io::IO, mm::SurvivalEstimator)
    println(io, typeof(mm))
    println(io)
    println(io,"Coefficients:")
    pretty_table(io, [mm.stats.n mm.stats.ncens mm.stats.nevents],
        header = ["n", "ncens", "nevents"], vlines = :none, hlines = :none)
    nothing
end

function Base.show(io::IO, mm::StatsModels.TableStatisticalModel{<:SurvivalEstimator})
    println(io, typeof(mm))
    println(io)
    println(io, mm.mf.f)
    println(io)
    println(io, "Coefficients:")
    pretty_table(io, [mm.model.stats.n mm.model.stats.ncens mm.model.stats.nevents],
        header = ["n", "ncens", "nevents"], vlines = :none, hlines = :none)
    nothing
end

#-------------------
# KaplanMeier
#-------------------
"""
    KaplanMeier <: SurvivalEstimator

See [`SurvivalEstimator`](@ref).
"""
mutable struct KaplanMeier <: SurvivalEstimator
    time::Vector{Float64}
    survival::Vector{Float64}
    std::Vector{Float64}
    distr::DiscreteNonParametric
    stats::NamedTuple

    KaplanMeier() = new()
end

"""
    kaplan_meier(Y::RCSurv)
    kaplan_meier(f::@formula, data::DataFrames.DataFrame)
    kaplan_meier(X::AbstractMatrix{<:Real}, Y::RCSurv)
    fit(KaplanMeier, X::AbstractMatrix{<:Real}, Y::RCSurv)

    Aliases: km, kaplan

Fit a non-parametric Kaplan-Meier estimator on a right-censored survival outcome. See [Fitting and predicting](@ref) and examples below for fitting interfaces, we recommend using `kaplan_meier(::@formula...)` (or aliases).

The Kaplan-Meier Estimator is defined by

```math
Ŝ(τ) = ∏_{i:tᵢ≤τ} (1 - \\frac{dᵢ}{nᵢ})
```

where ``dᵢ`` and ``nᵢ`` are the number of events and nunber at risk at time ``tᵢ`` respectively.

Standard deviation is calculated using the pointwise method of Kalbfleisch and Prentice (1980).

Future additions:

* Support for other censoring types (see [#46](@ref))
* Support for stratified estimator (see [#47](@ref))

Function returns a [`KaplanMeier`](@ref) struct.

# Examples
```jldoctest
julia> data = DataFrame(Y = [1,1,4,6,8,4,9,4,5,10], D = [true, false, false, false, true, false, false, true, true, false]);

julia> f = kaplan_meier(@formula(Srv(Y, D) ~ 1), data)
StatsModels.TableStatisticalModel{KaplanMeier, Matrix{Float64}}

(Y,D;+) ~ 1

Coefficients:
  n  ncens  nevents
 10      6        4
```
"""
kaplan_meier(args...; kwargs...) = StatsBase.fit(KaplanMeier, args...; kwargs...)
kaplan_meier(Y::RCSurv) = StatsBase.fit(KaplanMeier, Y)
const kaplan = kaplan_meier
const km = kaplan_meier

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

#-------------------
# Log-rank tests
#-------------------
function update_wt(wt, nevents, nrisk, wtmethod)
    if wtmethod == :logrank
        return 1
    elseif wtmethod == :wilcoxon
        return nrisk
    elseif wtmethod == :tw
       return sqrt(nrisk)
    elseif wtmethod == :peto
        return wt * (1 - nevents / (nrisk + 1))
    else
        error("Unknown weight")
    end
end

"""
    logrank_test(Y::RCSurv...; wtmethod=:logrank)

Test the null hypothesis that two or more survival functions are identical.
"""
function logrank_test(Y::RCSurv...; wtmethod=:logrank)

    m = length(Y)
    A = merge(Y...)
    ti = unique_outcome_times(A)
    sta = A.stats
    st = [_padstats(y.stats, ti) for y in Y]

    u = zeros(m)
    V = zeros(m, m)
    wt = 1.0
    for i in eachindex(ti)
        d, n = sta.nevents[i], sta.nrisk[i]
        wt = update_wt(wt, d, n, wtmethod)
        if n == 0
            break
        end
        for j in 1:m
            dd, nnj = st[j].nevents[i], st[j].nrisk[i]
            rj = dd / nnj
            fj = nnj / n
            u[j] += wt * (dd - d*fj)
            for k in 1:m
                nnk = st[k].nrisk[i]
                fk = nnk / n
                q = j == k ? 1.0 : 0.0
                if n > 1
                    V[j,k] += wt^2 * (q - fj) * fk * d * (n - d) / (n - 1)
                end
            end
        end
    end

    # Chi-square statistic
    csq = u' * pinv(V) * u

    # Degrees of freedom
    dof = m - 1

    # P-value
    p = 1 - cdf(Chisq(dof), csq)

    return (stat=csq, dof=dof, pvalue=p)
end


"""
    confint(km::KaplanMeier; level::Float64 = 0.95)
    confint(km::KaplanMeier, t::Number; level::Float64 = 0.95)

Calculate the confidence interval (CI) around a fitted Kaplan-Meier estimate to `level`% confidence. If `t` provided then returns CI at that time, otherwise returns CI at all fitted times.
Standard deviation is calculated using the pointwise method of Kalbfleisch and Prentice (1980).

# Examples
```jldoctest
julia> data = DataFrame(Y = [1,1,4,6,8,4,9,4,5,10], D = [true, false, false, false, true, false, false, true, true, false]);

julia> confint(km(@formula(Srv(Y, D) ~ 1), data))
4-element Vector{Tuple{Float64, Float64}}:
 (0.473009271362049, 0.9852813933673431)
 (0.38088152320549545, 0.9425909522237038)
 (0.21830025822743343, 0.8691223292427415)
 (0.0700802713627666, 0.7534316354804488)
```
"""
function StatsBase.confint(km::KaplanMeier, t::Number; level::Float64 = 0.95)
    return _confint_npe(
        km, t, level,
        (E, q, w) -> map(x -> exp(-exp(x)), log(-log(E.survival[w])) ∓ (q * E.std[w]))
    )
end

#-------------------
# NelsonAalen
#-------------------
"""
    NelsonAalen <: SurvivalEstimator

See [`SurvivalEstimator`](@ref).
"""
mutable struct NelsonAalen <: SurvivalEstimator
    time::Vector{Float64}
    survival::Vector{Float64}
    std::Vector{Float64}
    distr::DiscreteNonParametric
    stats::NamedTuple

    NelsonAalen() = new()
end

"""
    nelson_aalen(Y::RCSurv)
    nelson_aalen(f::@formula, data::DataFrames.DataFrame)
    nelson_aalen(X::AbstractMatrix{<:Real}, Y::RCSurv)
    fit(NelsonAalen, X::AbstractMatrix{<:Real}, Y::RCSurv)

    Aliases: na, nelson

Fit a non-parametric Nelson-NelsonAalen estimator on a right-censored survival outcome. See [Fitting and predicting](@ref) and examples below for fitting interfaces, we recommend using `nelson_aalen(::@formula...)` (or aliases).

The Nelson-Aalen Estimator is defined by

```math
Ĥ(τ) = ∑_{i:tᵢ≤τ} \\frac{dᵢ}{nᵢ}
```

where ``dᵢ`` and ``nᵢ`` are the number of events and nunber at risk at time ``tᵢ`` respectively and Ĥ is the estimated cumulative hazard function. The survival function is then Ŝ = exp(-Ĥ).

Future additions:

* Support for other censoring types (see [#46](@ref))
* Support for stratified estimator (see [#47](@ref))

Function returns a [`NelsonAalen`](@ref) struct.

# Examples
```jldoctest
julia> data = DataFrame(Y = [1,1,4,6,8,4,9,4,5,10], D = [true, false, false, false, true, false, false, true, true, false]);

julia> f = nelson_aalen(@formula(Srv(Y, D) ~ 1), data)
StatsModels.TableStatisticalModel{NelsonAalen, Matrix{Float64}}

(Y,D;+) ~ 1

Coefficients:
  n  ncens  nevents
 10      6        4
```
"""
nelson_aalen(args...; kwargs...) = StatsBase.fit(NelsonAalen, args...; kwargs...)
nelson_aalen(Y::RCSurv) = StatsBase.fit(NelsonAalen, Y)
const na = nelson_aalen
const nelson = nelson_aalen

function StatsBase.fit!(obj::NelsonAalen, Y::RCSurv)
    return _fit_npe(
        obj,
        Y,
        (d, n) -> d / n,
        (d, n) -> (d * (n - d)) / (n^3),
        # makes use of NA being a plug-in estimator then converts to survival
        (p) -> exp.(-cumsum(p)),
        (v, p) -> .√v
    )
end

# Calculates pointwise confidence intervals for *survival predictions*
# Unclear if this even needs to be a different method, could just use the confint around
# survival for both NPEs
"""
    confint(na::NelsonAalen; level::Float64 = 0.95)
    confint(na::NelsonAalen, t::Number; level::Float64 = 0.95)


Calculate the confidence interval (CI) around a fitted Nelson-Aalen estimate to `level`% confidence. If `t` provided then returns CI at that time, otherwise returns CI at all fitted times.

# Examples
```jldoctest
julia> data = DataFrame(Y = [1,1,4,6,8,4,9,4,5,10], D = [true, false, false, false, true, false, false, true, true, false]);

julia> confint(na(@formula(Srv(Y, D) ~ 1), data), 10)
(0.23186692958834798, 0.9464141496254767)
```
"""
function StatsBase.confint(na::NelsonAalen, t::Number; level::Float64 = 0.95)
    return _confint_npe(
        na, t, level,
        (E, q, w) -> map(x -> min(1, max(0, exp(-x))), -log(E.survival[w]) ∓ (q * E.std[w]))
    )
end
