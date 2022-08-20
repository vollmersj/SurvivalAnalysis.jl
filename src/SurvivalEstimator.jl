#-------------------
# SurvivalEstimator
#-------------------
abstract type SurvivalEstimator <: StatisticalModel end

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
    obj.std = [0, std_trafo(cumsum(v), obj.survival)...]
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
    return confint.(Ref(npe.model), npe.model.time, level = level)
end

Base.time(npe::SurvivalEstimator) = npe.time
survival(npe::SurvivalEstimator) = npe.survival
StatsBase.std(npe::SurvivalEstimator) = npe.std
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
mutable struct KaplanMeier <: SurvivalEstimator
    time::Vector{Float64}
    survival::Vector{Float64}
    std::Vector{Float64}
    distr::DiscreteNonParametric
    stats::NamedTuple

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

function StatsBase.confint(km::KaplanMeier, t::Number; level::Float64 = 0.95)
    return _confint_npe(
        km, t, level,
        (E, q, w) -> map(x -> exp(-exp(x)), log(-log(E.survival[w])) ∓ (q * E.std[w]))
    )
end

#-------------------
# NelsonAalen
#-------------------
mutable struct NelsonAalen <: SurvivalEstimator
    time::Vector{Float64}
    survival::Vector{Float64}
    std::Vector{Float64}
    distr::DiscreteNonParametric
    stats::NamedTuple

    NelsonAalen() = new()
end

nelson_aalen(args...; kwargs...) = StatsBase.fit(NelsonAalen, args...; kwargs...)
nelson_aalen(Y::RCSurv) = StatsBase.fit(NelsonAalen, Y)

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
function StatsBase.confint(na::NelsonAalen, t::Number; level::Float64 = 0.95)
    return _confint_npe(
        na, t, level,
        (E, q, w) -> map(x -> min(1, max(0, exp(-x))), -log(E.survival[w]) ∓ (q * E.std[w]))
    )
end
