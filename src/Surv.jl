# module Survival
using Distributions
using StatsBase
using RecipesBase
using Plots
using DataFrames

abstract type Surv end
abstract type OneSidedSurv <: Surv end
abstract type TwoSidedSurv <:Surv end

struct rcSurv <: OneSidedSurv
    time::Vector{Float64}
    status::Vector{Bool}
    symbol::Char
    type::String
    rcSurv(time::Union{Vector{Float64}, Vector{Number}}, status::Union{Vector{Bool},BitVector, Vector{Int}}) =
        new(convert(Vector{Float64}, time), convert(BitVector, status), '+', "right")
    rcSurv(time::Number, status::Bool) =
        new([convert(Float64, time)], [status], '+', "right")
    rcSurv(time::Number, status::Int) =
        new([convert(Float64, time)], [status == 1], '+', "right")
end

struct lcSurv <: OneSidedSurv
    time::Vector{Float64}
    status::Vector{Bool}
    symbol::Char
    type::String
    lcSurv(time::Vector{Number}, status::Union{Vector{Bool},BitVector}) =
        new(convert(Vector{Float64}, time), status, '-', "left")
    lcSurv(time::Number, status::Bool) =
        new([convert(Float64, time)], [status], '-', "left")
    lcSurv(time::Number, status::Int) =
        new([convert(Float64, time)], [status == 1], '-', "left")
end

struct intSurv <: TwoSidedSurv ## FIXME - VECTORISE LIKE OneSidedSurv
    start::Float64
    stop::Float64
    type::String
    intSurv(start, stop) = new(start, stop, "interval")
end

function Surv(start::AbstractFloat, stop::AbstractFloat)
    intSurv(start, stop)
end

function Surv(time::Vector{Float64}, status::BitVector, type::String)
    @assert type in ["left", "right"]
    if type == "right"
        rcSurv(time, status)
    elseif type == "left"
        lcSurv(time, status)
    end
end

Base.show(io::IO, oss::OneSidedSurv) =
    print(io, map((t, δ) -> string(t, δ ? "" : "+"), oss.time, oss.status))
Base.show(io::IO, oss::TwoSidedSurv) =
    print(io, "(", oss.start, ", ", oss.stop, "]")

outcomeTimes(v::Union{rcSurv,lcSurv}) = v.time
outcomeTimes(v::Vector{intSurv}) = map(x -> [x.start, x.stop], v)

eventTimes(v::Union{rcSurv,lcSurv}) = v.time[v.status]
eventTimes(v::Vector{intSurv}) = map(x -> x.stop, v)

outcomeStatus(v::Union{rcSurv, lcSurv}) = v.status

uniqueTimes(v::Union{rcSurv, lcSurv}) = sort(unique(outcomeTimes(v)))
uniqueEventTimes(v::Union{rcSurv, lcSurv}) = sort(unique(eventTimes(v)))

riskSet(v::Union{rcSurv,lcSurv}, t::Number) = map(x -> x >= t, v.time)

totalEvents(v::Union{rcSurv,lcSurv}) = length(v.status)
totalEvents(v::Union{rcSurv,lcSurv}, t::Number) = sum(map(x -> x == t, v.time))

totalDeaths(v::Union{rcSurv,lcSurv}) = sum(v.outcome.status)
totalDeaths(v::Union{rcSurv,lcSurv}, t::Number) = sum(map((τ, δ) -> τ == t &&  δ, v.time, v.status))

totalRisk(v::Union{rcSurv,lcSurv}, t::Number) = sum(riskSet(v, t))

abstract type NonParametricEstimator end
struct KaplanMeier <: NonParametricEstimator
    times::Vector{Float64}
    survs::Vector{Float64}
    sd::Vector{Float64}
end

struct NelsonAalen <: NonParametricEstimator
    times::Vector{Float64}
    survs::Vector{Float64}
    sd::Vector{Float64}
end

function fit_NPE(class, Surv::rcSurv, point_est::Function, var_est::Function,
                surv_trafo::Function, sd_trafo::Function)
    ut = uniqueEventTimes(Surv)
    p = zeros(length(ut))
    v = zeros(length(ut))
    for (i, t) in enumerate(ut)
        d = totalDeaths(Surv, t)
        n = sum(riskSet(Surv, t))
        p[i] = point_est(d, n)
        v[i] = var_est(d, n)
    end
    variance = cumsum(v)
    surv = surv_trafo(p)
    sd = sd_trafo(variance, surv)
    class(ut, surv, sd)
end

function kaplan(Surv::rcSurv)
    fit_NPE(
        KaplanMeier,
        Surv,
        (d, n) -> 1 - (d / n),
        (d, n) -> d / (n * (n - d)),
        cumprod,
        # Calculates standard error using the method of Kalbfleisch and Prentice (1980)
        (v, p) -> map((qₜ, vₜ) -> √(vₜ / (log(qₜ)^2)), p, v)
    )
end

function nelson(Surv::rcSurv)
    fit_NPE(
        NelsonAalen,
        Surv,
        (d, n) -> d / n,
        (d, n) -> (d * (n - d)) / (n^3),
        (p) -> exp.(-cumsum(p)),
        (v, p) -> .√v
    )
end

±(x, y) = (x - y, x + y)
∓(x, y) = (x + y, x - y)

# Calculates confidence intervals using the method of Kalbfleisch and Prentice (1980)
function confint_NPE(npe, t, α, trafo)
    q = quantile(Normal(), 1 - α/2)
    which = searchsortedlast(npe.times, t)
    trafo(npe, q, which)
end

function StatsBase.confint(km::KaplanMeier, t::Number, α::Float64 = 0.05)
    confint_NPE(
        km, t, α,
        (E, q, w) -> map(x -> exp(-exp(x)), log(-log(E.survs[w])) ∓ (q * E.sd[w]))
    )
end

# Calculates pointwise confidence intervals for *survival predictions*
function StatsBase.confint(na::NelsonAalen, t::Number, α::Float64 = 0.05)
    confint_NPE(
        na, t, α,
        (E, q, w) -> map(x -> min(1, max(0, exp(-x))), -log(E.survs[w]) ∓ (q * E.sd[w]))
    )
end

function survival(npe::NonParametricEstimator, t)
    npe.survs[searchsortedlast(npe.times, t)]
end
function chf(npe::NonParametricEstimator, t)
    -log(survival(npe, t))
end
function Distributions.cdf(npe::NonParametricEstimator, t)
    1 - survival(npe, t)
end

@recipe function f(npe::NonParametricEstimator, plot_confint::Bool = true, α::Float64 = 0.05)
    seriestype := :steppost
    ylims := (0, 1)
    legend := false
    @series begin
        linecolor   --> :black
        npe.times, npe.survs
    end
    if plot_confint
        linecolor   --> :blue
        cis = confint.(Ref(npe), npe.times, α)
        lb = map(x -> x[1], cis)
        ub = map(x -> x[2], cis)
        @series begin
            npe.times, lb
        end
        @series begin
            npe.times, ub
        end
    end
end
