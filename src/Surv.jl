abstract type Surv end
abstract type OneSidedSurv <: Surv end
abstract type TwoSidedSurv <:Surv end

struct rcSurv <: OneSidedSurv
    time::Vector{Float64}
    status::Vector{Bool}
    symbol::Char
    type::String
    stats::NamedTuple{(:time, :nrisk, :ncens, :nevents, :noutcomes), Tuple{Vector{Float64}, Vector{Int64}, Vector{Int64}, Vector{Int64}, Vector{Int64}}}
    rcSurv(time::Vector{Float64}, status::BitVector) =
        new(time, status, '+', "right", _tabulateSurv(time, status))
end

struct lcSurv <: OneSidedSurv
    time::Vector{Float64}
    status::Vector{Bool}
    symbol::Char
    type::String
    stats::NamedTuple{(:time, :nrisk, :ncens, :nevents, :noutcomes), Tuple{Vector{Float64}, Vector{Int64}, Vector{Int64}, Vector{Int64}, Vector{Int64}}}
    lcSurv(time::Vector{Float64}, status::BitVector) =
        new(time, status, '-', "left", _tabulateSurv(time, status))
end

struct intSurv <: TwoSidedSurv
    start::Vector{Float64}
    stop::Vector{Float64}
    type::String
    intSurv(start::Vector{Float64}, stop::Vector{Float64}) = new(start, stop, "interval")
end

function Surv(start::Union{Vector{T}, T} where T <: Number,
            stop::Union{Vector{T}, T} where T <: Number)
    start = start isa Vector ? convert(Vector{Float64}, start) : convert(Vector{Float64}, [start])
    stop = stop isa Vector ? convert(Vector{Float64}, stop) : convert(Vector{Float64}, [stop])
    intSurv(start, stop)
end

function Surv(time::Union{Vector{T}, T} where T <: Number,
            status::Union{BitVector, Vector{Bool}, Bool, Int, Vector{Int}},
            type::String)
    @assert type in ["left", "right"]

    time = time isa Vector ? convert(Vector{Float64}, time) : convert(Vector{Float64}, [time])
    status = (status isa Bool || status isa Int) ? convert(BitVector, [status]) :
        convert(BitVector, status)

    type == "right" ? rcSurv(time, status) : lcSurv(time, status)
end

Base.show(io::IO, oss::OneSidedSurv) =
    print(io, map((t, δ) -> string(t, δ ? "" : "+"), oss.time, oss.status))
Base.show(io::IO, oss::TwoSidedSurv) =
    print(io, "(", oss.start, ", ", oss.stop, "]")

outcomeTimes(v::OneSidedSurv) = v.time
outcomeTimes(v::TwoSidedSurv) = [v.start, v.stop]

eventTimes(v::OneSidedSurv) = v.time[v.status]
eventTimes(v::TwoSidedSurv) = v.stop

outcomeStatus(v::OneSidedSurv) = v.status

uniqueTimes(v::OneSidedSurv) = v.stats.time
uniqueEventTimes(v::OneSidedSurv) = unique(eventTimes(v))

totalEvents(v::OneSidedSurv) = sum(v.status)
totalCensored(v::OneSidedSurv) = sum(!v.status)
totalOutcomes(v::OneSidedSurv) = length(v.status)
totalRisk(v::OneSidedSurv) = length(v.status)

# unsure if should return 0, NA, or error if > observed times
#  also unclear about behaviour between times
function _totalOutcome(v::OneSidedSurv, t::Number, var::String)
    if t == 0
        return var == "nrisk" ? v.stats.time[1] : 0
    elseif t > maximum(v.stats.time)
        return 0
    else
        getproperty(v.stats, Symbol(var))[searchsortedlast(v.stats.time, t)]
    end
end

totalEvents(v::OneSidedSurv, t::Number) = _totalOutcome(v, t, "nevents")
totalCensored(v::OneSidedSurv, t::Number) = _totalOutcome(v, t, "ncens")
totalOutcomes(v::OneSidedSurv, t::Number) = _totalOutcome(v, t, "noutcomes")
totalRisk(v::OneSidedSurv, t::Number) = _totalOutcome(v, t, "nrisk")

survStats(v::OneSidedSurv; events_only = false) =
    events_only ? (w = map(x -> x > 0, v.stats.nevents); map(s -> s[w], v.stats)) : v.stats

function _tabulateSurv(T, Δ)
    ut = sort(unique(T))
    n = length(ut)
    nrisk = zeros(n)
    ncens = zeros(n)
    nevents = zeros(n)
    noutcomes = zeros(n)
    for (i, t) in enumerate(ut)
        nrisk[i] = sum(map(x -> x >= t, T))
        ncens[i] = sum(map((τ, δ) -> τ == t &&  !δ, T, Δ))
        nevents[i] = sum(map((τ, δ) -> τ == t &&  δ, T, Δ))
        noutcomes[i] = ncens[i] + nevents[i]
    end
    (time = ut, nrisk = nrisk, ncens = ncens, nevents = nevents, noutcomes = noutcomes)
end

function Base.merge(A::OneSidedSurv...)
    type = map(v -> v.type, A)
    @assert length(unique(type)) == 1
    T = zeros(0)
    Δ = falses(0)
    foreach(v -> push!(T, v.time...), A)
    foreach(v -> push!(Δ, v.status...), A)
    Surv(T, Δ, type[1])
end

function Base.merge(A::TwoSidedSurv...)
    start = zeros(0)
    stop = zeros(0)
    foreach(v -> push!(start, v.start...), A)
    foreach(v -> push!(stop, v.stop...), A)
    Surv(start, stop)
end
