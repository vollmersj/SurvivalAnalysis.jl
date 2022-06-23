abstract type Surv end
abstract type OneSidedSurv <: Surv end
abstract type TwoSidedSurv <:Surv end

struct rcSurv <: OneSidedSurv
    time::Float16
    status::Bool
    symbol::Char
    type::String
    rcSurv(time, status) = new(time, status, '+', "right")
end

struct lcSurv <: OneSidedSurv
    time::Float16
    status::Bool
    symbol::Char
    type::String
    lcSurv(time, status) = new(time, status, '-', "left")
end

struct intSurv <: TwoSidedSurv
    start::Float16
    stop::Float16
    type::String
    intSurv(start, stop) = new(start, stop, "interval")
end

function Surv(start::AbstractFloat, stop::AbstractFloat)
    intSurv(start, stop)
end

function Surv(time::AbstractFloat, status::Bool, type::String)
    @assert type in ["left", "right"]
    if type == "right"
        rcSurv(time, status)
    elseif type == "left"
        lcSurv(time, status)
    end
end

Base.show(io::IO, oss::OneSidedSurv) =
    print(io, oss.time, oss.status ? "" : oss.symbol)

Base.show(io::IO, oss::TwoSidedSurv) =
    print(io, "(", oss.start, ", ", oss.stop, "]")

outcomeTimes(v::Vector{<:Union{rcSurv,lcSurv}}) = map(x -> x.time, v)
function eventTimes(v::Vector{<:Union{rcSurv,lcSurv}})
    res = [];
    foreach(x -> x.status && push!(res, x.time), v)
    res
end
outcomeTimes(v::Vector{intSurv}) = map(x -> [x.start, x.stop], v)
eventTimes(v::Vector{intSurv}) = map(x -> x.stop, v)

outcomeStatus(v::Vector{<:Union{rcSurv, lcSurv}}) = map(x -> x.status, v)
uniqueTimes(v::Vector) = sort(unique(outcomeTimes(v)))
uniqueEventTimes(v::Vector) = sort(unique(eventTimes(v)))
riskSet(v::Vector{<:Union{rcSurv,lcSurv}}, t::Number) = map(x -> x.time >= t, v)
totalDeaths(v::Vector{<:Union{rcSurv,lcSurv}}, t::Number) = sum(map(x -> x.status && x.time == t, v))
totalEvents(v::Vector{<:Union{rcSurv,lcSurv}}, t::Number) = sum(map(x -> x.time == t, v))

abstract type NonParametricEstimator end
struct KaplanMeier <: NonParametricEstimator
    times::Vector{Float16}
    survs::Vector{Float16}
end

struct NelsonAalen <: NonParametricEstimator
    times::Vector{Float16}
    survs::Vector{Float16}
end

kaplan = function (v::Vector{rcSurv})
    ut = uniqueTimes(survs)
    surv = []
    for tmax in ut
        p = 1
        sd = 0
        for t in ut
            if t <= tmax
                d = totalDeaths(v, t)
                n = totalEvents(v, t)
                p *= 1 - (d / n)
                sd += √(d / (n * (n - d)))
            end
        end
        push!(surv, p)
    end
    KaplanMeier(ut, surv)
end

nelson = function (v::Vector{rcSurv})
    ut = uniqueTimes(survs)
    chf = []
    for tmax in ut
        p = 0
        for t in ut
            if t <= tmax
                p += (totalDeaths(v, t) / totalEvents(v, t))
            end
        end
        push!(chf, p)
    end
    NelsonAalen(ut, exp.(-chf))
end

function survival(npe::NonParametricEstimator, t)
    npe.survs[searchsortedlast(npe.times, t)]
end
function chf(npe::NonParametricEstimator, t)
    -log(survival(npe, t))
end
function cdf(npe::NonParametricEstimator, t)
    1 - survival(npe, t)
end


survival.(Ref(km), [1,2, 3, 8])
