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

# Calculates standard error using the method of Kalbfleisch and Prentice (1980)
function kaplan(v::rcSurv)
    ut = uniqueEventTimes(v)
    surv = []
    # V̂ = []
    V̂ = ones(length(ut))
    for tmax in ut
        p = 1
        # Vᵢ = 0
        for t in ut
            if t <= tmax
                d = totalDeaths(v, t)
                n = sum(riskSet(v, t))
                p *= (1 - (d / n))
                # Vᵢ += d / (n * (n - d))
            end
        end
        push!(surv, p)
        # push!(V̂, Vᵢ / (log(p)^2))
    end
    KaplanMeier(ut, surv, .√V̂)
end

function kaplan2(v::rcSurv)
    ut = uniqueEventTimes(v)
    q = zeros(length(ut))
    V̂ = ones(length(ut))
    for (i, t) in enumerate(ut)
        p = 1
        # Vᵢ = 0
        d = totalDeaths(v, t)
        n = sum(riskSet(v, t))
        q[i] = 1 - (d / n)
        # push!(V̂, Vᵢ / (log(p)^2))
    end
    KaplanMeier(ut, cumprod(q), .√V̂)
end

function nelson(v::Vector{rcSurv})
    ut = uniqueEventTimes(v)
    chf = []
    sd = []
    for tmax in ut
        p = 0
        sdᵢ = 0
        for t in ut
            if t <= tmax
                d = totalDeaths(v, t)
                n = sum(riskSet(v, t))
                p += (d / n)
                sdᵢ += (d * (n - d)) / (n^3)
            end
        end
        push!(chf, p)
        push!(sd, √(sdᵢ))
    end
    NelsonAalen(ut, exp.(-chf), sd)
end

±(x, y) = (x - y, x + y)
∓(x, y) = (x + y, x - y)

# Calculates confidence intervals using the method of Kalbfleisch and Prentice (1980)
function StatsBase.confint(km::KaplanMeier, t::Number, α::Float64 = 0.05)
    q = quantile(Normal(), 1 - α/2)
    which = searchsortedlast(km.times, t)
    map(x -> exp(-exp(x)), log(-log(km.survs[which])) ∓ (q * km.sd[which]))
end

# Calculates pointwise confidence intervals for *survival predictions*
function StatsBase.confint(na::NelsonAalen, t::Number, α::Float64 = 0.05)
    q = quantile(Normal(), 1 - α/2)
    which = searchsortedlast(na.times, t)
    map(x -> min(1, max(0, exp(-x))),
        -log(na.survs[which]) ∓ (q * na.sd[which]))
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
