abstract type Surv end
abstract type OneSidedSurv <: Surv end
abstract type TwoSidedSurv <:Surv end

struct RCSurv <: OneSidedSurv
    time::Vector{Float64}
    status::Vector{Bool}
    symbol::Char
    type::String
    stats::NamedTuple{(:time, :nrisk, :ncens, :nevents, :noutcomes), Tuple{Vector{Float64}, Vector{Int64}, Vector{Int64}, Vector{Int64}, Vector{Int64}}}
    RCSurv(time::Vector{Float64}, status::BitVector) =
        new(time, status, '+', "right", _tabulate_surv(time, status))
end

struct LCSurv <: OneSidedSurv
    time::Vector{Float64}
    status::Vector{Bool}
    symbol::Char
    type::String
    stats::NamedTuple{(:time, :nrisk, :ncens, :nevents, :noutcomes), Tuple{Vector{Float64}, Vector{Int64}, Vector{Int64}, Vector{Int64}, Vector{Int64}}}
    LCSurv(time::Vector{Float64}, status::BitVector) =
        new(time, status, '-', "left", _tabulate_surv(time, status))
end

struct IntSurv <: TwoSidedSurv
    start::Vector{Float64}
    stop::Vector{Float64}
    type::String
    IntSurv(start::Vector{Float64}, stop::Vector{Float64}) = new(start, stop, "interval")
end

function Surv(start::Union{Vector{T}, T} where T <: Number,
            stop::Union{Vector{T}, T} where T <: Number)
    start = start isa Vector ? convert(Vector{Float64}, start) : convert(Vector{Float64}, [start])
    stop = stop isa Vector ? convert(Vector{Float64}, stop) : convert(Vector{Float64}, [stop])
    IntSurv(start, stop)
end

function Surv(time::Union{Vector{T}, T} where T <: Number,
            status::Union{BitVector, Vector{Bool}, Bool, Int, Vector{Int}},
            type::String)
    @assert type in ["left", "right"]

    time = time isa Vector ? convert(Vector{Float64}, time) : convert(Vector{Float64}, [time])
    status = (status isa Bool || status isa Int) ? convert(BitVector, [status]) :
        convert(BitVector, status)

    type == "right" ? RCSurv(time, status) : LCSurv(time, status)
end

Base.show(io::IO, oss::OneSidedSurv) =
    print(io, map((t, δ) -> string(t, δ ? "" : "+"), oss.time, oss.status))
Base.show(io::IO, oss::TwoSidedSurv) =
    print(io, "(", oss.start, ", ", oss.stop, "]")

outcome_times(v::OneSidedSurv) = v.time
outcome_times(v::TwoSidedSurv) = [v.start, v.stop]

event_times(v::OneSidedSurv) = v.time[v.status]
event_times(v::TwoSidedSurv) = v.stop

outcome_status(v::OneSidedSurv) = v.status

unique_times(v::OneSidedSurv) = v.stats.time
unique_event_times(v::OneSidedSurv) = unique(event_times(v))

total_events(v::OneSidedSurv) = sum(v.status)
total_censored(v::OneSidedSurv) = sum(!v.status)
total_outcomes(v::OneSidedSurv) = length(v.status)
total_risk(v::OneSidedSurv) = length(v.status)

# unsure if should return 0, NA, or error if > observed times
#  also unclear about behaviour between times
function _total_outcome(v::OneSidedSurv, t::Number, var::String)
    if t == 0
        return var == "nrisk" ? v.stats.time[1] : 0
    elseif t > maximum(v.stats.time)
        return 0
    else
        getproperty(v.stats, Symbol(var))[searchsortedlast(v.stats.time, t)]
    end
end

total_events(v::OneSidedSurv, t::Number) = _total_outcome(v, t, "nevents")
total_censored(v::OneSidedSurv, t::Number) = _total_outcome(v, t, "ncens")
total_outcomes(v::OneSidedSurv, t::Number) = _total_outcome(v, t, "noutcomes")
total_risk(v::OneSidedSurv, t::Number) = _total_outcome(v, t, "nrisk")

surv_stats(v::OneSidedSurv; events_only = false) =
    events_only ? (w = map(x -> x > 0, v.stats.nevents); map(s -> s[w], v.stats)) : v.stats

function _tabulate_surv(T, Δ)
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
