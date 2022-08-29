"""
    Surv
    OneSidedSurv <: Surv
    TwoSidedSurv <:Surv
    RCSurv <: OneSidedSurv
    LCSurv <: OneSidedSurv
    IntSurv <: TwoSidedSurv

`Surv` is the abstract type for all survival outcome representations in this package. See [`Surv`](@ref) for constructing a `Surv` object.

Available methods:

* [`length`](@ref) - Get length of a survival object
* [`outcome_times`](@ref) - Get times at which an outcome (event or censoring) takes place
* [`event_times`](@ref) - Get times at which an event (not censoring) takes place
* [`outcome_status`](@ref) - Get vector of survival indicators
* [`unique_outcome_times`](@ref) - Get unique times at which an outcome takes place
* [`unique_event_times`](@ref) - Get unique times at which an event takes place
* [`total_events`](@ref) - Get total number of events (optionally at a given time)
* [`total_censored`](@ref) - Get total number censored (optionally at a given time)
* [`total_outcomes`](@ref) - Get total number of outcomes (optionally at a given time)
* [`total_risk`](@ref) -  Get total number at risk (optionally at a given time)
* [`surv_stats`](@ref) - Get set of useful summary statistics
* [`threshold_risk`](@ref) - Get the time at which a given proportion of observations are no longer at risk
* [`merge`](@ref) - Merge survival objects
* [`reverse`](@ref) - Reverse survival outcome (not in-place)
"""
abstract type Surv end
abstract type OneSidedSurv <: Surv end
abstract type TwoSidedSurv <:Surv end

struct RCSurv <: OneSidedSurv
    time::Vector{Float64}
    status::Vector{Bool}
    symbol::Char
    stats::NamedTuple{(:time, :nrisk, :ncens, :nevents, :noutcomes),
        Tuple{Vector{Float64}, Vector{Int64}, Vector{Int64}, Vector{Int64}, Vector{Int64}}}
    RCSurv(time::Vector{Float64}, status::BitVector) =
        new(time, status, '+', _tabulate_surv(time, status))
end

struct LCSurv <: OneSidedSurv
    time::Vector{Float64}
    status::Vector{Bool}
    symbol::Char
    stats::NamedTuple{(:time, :nrisk, :ncens, :nevents, :noutcomes),
        Tuple{Vector{Float64}, Vector{Int64}, Vector{Int64}, Vector{Int64}, Vector{Int64}}}
    LCSurv(time::Vector{Float64}, status::BitVector) =
        new(time, status, '-', _tabulate_surv(time, status))
end

struct IntSurv <: TwoSidedSurv
    start::Vector{Float64}
    stop::Vector{Float64}
    IntSurv(start::Vector{Float64}, stop::Vector{Float64}) = new(start, stop)
end

"""
    Surv(start, stop)
    Surv(time)
    Surv(time, status, type)

Survival analysis is dependent on representing both known and unknown survival time outcomes for an observation. When the true survival time is unknown we say an observation is *censored* and we instead record their censoring time. Let ``(T, Δ)`` be a one-sided survival outcome, then ``T`` is the outcome time (the time at which an outcome, event or censoring, is observed) and ``Δ`` is the survival indicator (1 if the outcome is an event and 0 if censoring). For example, if a patient drops out of a study at time 5, then they are recorded as ``(T=5, Δ=0)``. Two-sided survival outcomes are more simply recorded as ``(T₁,T₂)`` which means the true event time was somewhere between ``T₁`` and ``T₂``.

There are three censoring types:

* Right (`Surv(time, status, :r)` - The true outcome time occurs at some time *after* the observed censoring time
* Left (`Surv(time, status, :l)`) - The true outcome time occurs at some time *before* the observed censoring time
* Interval (`Surv(start, stop)`) - The true outcome time occurs at some time *within* the observed censoring times

If no `status` vector is passed to the function then it is assumed no-one is censored - this is rarely useful in practice as in this case one could simply use regression models.

Note❗ Whilst this package supports functionality for all censoring types, currently only methods for right censoring are included.

# Examples
```jldoctest
julia> Surv([1, 2, 3], [true, false, true], :r) # right-censoring
["1.0", "2.0+", "3.0"]

julia> Surv([1, 2, 3], [1, 0, 1], :l) # left-censoring
["1.0", "2.0-", "3.0"]

julia> Surv([1, 2, 3], [5, 6, 7]) # interval-censoring
[(1.0, 5.0), (2.0, 6.0), (3.0, 7.0)]
```
"""
function Surv(start::Union{Vector{T}, T} where T <: Number,
            stop::Union{Vector{T}, T} where T <: Number)
    start = start isa Vector ? convert(Vector{Float64}, start) :
        convert(Vector{Float64}, [start])
    stop = stop isa Vector ? convert(Vector{Float64}, stop) :
        convert(Vector{Float64}, [stop])
    return IntSurv(start, stop)
end

Surv(time::Union{Vector{T}, T} where T <: Number) = Surv(time, trues(length(time)), :right)


function Surv(time::Union{Vector{T}, T} where T <: Number,
            status::Union{BitVector, Vector{Bool}, Bool, Int, Vector{Int}},
            type::Symbol)

    type in [:left, :right, :l, :r] ||
        throw(ArgumentError("`type` must be one of `:right`, `:r`, `:left`, `:l`"))

    time = time isa Vector ? convert(Vector{Float64}, time) :
        convert(Vector{Float64}, [time])
    status = (status isa Bool || status isa Int) ? convert(BitVector, [status]) :
        convert(BitVector, status)

    return (type === :right || type === :r) ? RCSurv(time, status) : LCSurv(time, status)
end

function Base.show(io::IO, srv::OneSidedSurv)
    print(io,
        map((t, δ) -> string(round(t; digits=3), δ ? "" : srv.symbol), srv.time, srv.status)
    )
end

function Base.show(io::IO, srv::TwoSidedSurv)
    print(io,
        map((t0, t1) -> (round(t0; digits=3), round(t1; digits=3)), srv.start, srv.stop)
    )
end

"""
    length(srv::OneSidedSurv)
    length(srv::TwoSidedSurv)

Get length of a survival object.

# Examples
```jldoctest
julia> length(Surv([1, 2, 3], [true, false, true], :r))
3

julia> length(Surv([1, 2, 3], [5, 6, 7]))
3
```
"""
Base.length(srv::OneSidedSurv) = length(srv.time)
Base.length(srv::TwoSidedSurv) = length(srv.start)

"""
    outcome_times(srv::OneSidedSurv)
    outcome_times(srv::TwoSidedSurv)

    Get times at which an outcome (event or censoring) takes place.

# Examples
```jldoctest
julia> outcome_times(Surv([1, 2, 3], [true, false, true], :r))
3-element Vector{Float64}:
 1.0
 2.0
 3.0

julia> outcome_times(Surv([1, 2, 3], [5, 6, 7]))
2-element Vector{Vector{Float64}}:
 [1.0, 2.0, 3.0]
 [5.0, 6.0, 7.0]
```
"""
outcome_times(srv::OneSidedSurv) = srv.time
outcome_times(srv::TwoSidedSurv) = [srv.start, srv.stop]

"""
    outcome_times(srv::OneSidedSurv)

Get times at which an event (not censoring) takes place.

# Examples
```jldoctest
julia> event_times(Surv([1, 2, 3], [true, false, true], :r))
2-element Vector{Float64}:
 1.0
 3.0
```
"""
event_times(srv::OneSidedSurv) = srv.time[srv.status]

"""
    outcome_status(srv::OneSidedSurv)

    Get vector of survival indicators.

# Examples
```jldoctest
julia> outcome_status(Surv([1, 2, 3], [true, false, true], :r))
3-element Vector{Bool}:
 1
 0
 1
```
"""
outcome_status(srv::OneSidedSurv) = srv.status

"""
    unique_outcome_times(srv::OneSidedSurv)

Get unique times at which an outcome takes place

# Examples
```jldoctest
julia> unique_outcome_times(Surv([1, 2, 1], [true, false, true], :r))
2-element Vector{Float64}:
 1.0
 2.0
```
"""
unique_outcome_times(srv::OneSidedSurv) = srv.stats.time

"""
    unique_event_times(srv::OneSidedSurv)

Get unique times at which an event takes place.

# Examples
```jldoctest
julia> unique_event_times(Surv([1, 2, 1], [true, false, true], :r))
1-element Vector{Float64}:
 1.0
```
"""
unique_event_times(srv::OneSidedSurv) = unique(event_times(srv))

"""
    total_events(srv::OneSidedSurv)
    total_events(srv::OneSidedSurv, t::Number)

Get total number of events (optionally at a given time).

# Examples
```jldoctest
julia> total_events(Surv([1, 2, 3], [true, false, true], :r))
2

julia> total_events(Surv([1, 2, 3], [true, false, true], :r), 3)
1
```
"""
total_events(srv::OneSidedSurv) = sum(srv.status)

"""
    total_censored(srv::OneSidedSurv)
    total_censored(srv::OneSidedSurv, t::Number)

Get total number censored (optionally at a given time).

# Examples
```jldoctest
julia> total_censored(Surv([1, 2, 3], [false, true, false], :r))
2

julia> total_censored(Surv([1, 2, 3], [true, false, true], :r), 3)
1
```
"""
total_censored(srv::OneSidedSurv) = length(srv.status) - sum(srv.status)

"""
    total_outcomes(srv::OneSidedSurv)
    total_outcomes(srv::OneSidedSurv, t::Number)

Get total number of outcomes (optionally at a given time).

# Examples
```jldoctest
julia> total_outcomes(Surv([1, 2, 3], [false, true, false], :r))
3

julia> total_outcomes(Surv([1, 2, 3], [true, false, true], :r), 3)
1
```
"""
total_outcomes(srv::OneSidedSurv) = length(srv.status)

"""
    total_risk(srv::OneSidedSurv)
    total_risk(srv::OneSidedSurv, t::Number)

Get total number at risk (optionally at a given time).

# Examples
```jldoctest
julia> total_risk(Surv([1, 2, 3], [false, true, false], :r))
3

julia> total_risk(Surv([1, 2, 3], [true, false, true], :r), 2)
2
```
"""
total_risk(srv::OneSidedSurv) = length(srv.status)

# TODO - unsure if should return 0, NA, or error if > observed times
#  also unclear about behaviour between times
function _total_outcome(srv::OneSidedSurv, t::Number, var::String)
    if t < minimum(srv.stats.time)
        return var == "nrisk" ? srv.stats.nrisk[1] : 0
    elseif t > maximum(srv.stats.time)
        return 0
    else
        return getproperty(srv.stats, Symbol(var))[searchsortedlast(srv.stats.time, t)]
    end
end

total_events(srv::OneSidedSurv, t::Number) = _total_outcome(srv, t, "nevents")
total_censored(srv::OneSidedSurv, t::Number) = _total_outcome(srv, t, "ncens")
total_outcomes(srv::OneSidedSurv, t::Number) = _total_outcome(srv, t, "noutcomes")
total_risk(srv::OneSidedSurv, t::Number) = _total_outcome(srv, t, "nrisk")

"""
    surv_stats(srv::OneSidedSurv; events_only = false)

Get set of useful summary statistics over time - if `events_only = false` (default) then returned for all outcomes, otherwise only return the times at which an event occurred (useful for computing non-parametric estimators).

# Examples
```jldoctest
julia> surv_stats(Surv([1, 2, 3], [false, true, false], :r))
(time = [1.0, 2.0, 3.0], nrisk = [3, 2, 1], ncens = [1, 0, 1], nevents = [0, 1, 0], noutcomes = [1, 1, 1])

julia> surv_stats(Surv([1, 2, 3], [false, true, false], :r); events_only = true)
(time = [2.0], nrisk = [2], ncens = [0], nevents = [1], noutcomes = [1])
```
"""
function surv_stats(srv::OneSidedSurv; events_only = false)
    events_only ? (w = map(x -> x > 0, srv.stats.nevents); map(s -> s[w], srv.stats)) :
        srv.stats
end

"""
    threshold_risk(srv::OneSidedSurv, p::Number)

Get the time at which a given proportion of observations are no longer at risk. Useful when computing metrics and calculations become increasingly unstable as number of observations at risk decreases over time.

# Examples
```jldoctest
# time at which 80% of observations have experienced the event or been censored
julia> threshold_risk(Surv([1, 2, 3, 4, 5, 6], [false, true, false, false, true, true], :r), 0.8)
6.0
```
"""
function threshold_risk(srv::OneSidedSurv, p::Number)
    test_proportion(p) || throw(ArgumentError("Expected 0≤`p`≤1, got $(p)"))
    p == 1 && return srv.stats.time[end]
    srv.stats.time[searchsortedfirst(1 .- srv.stats.nrisk / length(srv), p)]
end

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
    return (time = ut, nrisk = nrisk, ncens = ncens, nevents = nevents,
            noutcomes = noutcomes)
end

"""
    merge(A::OneSidedSurv...)
    merge(A::TwoSidedSurv...)

    Merge survival objects.

# Examples
```jldoctest
julia> srv1 = Surv([1, 2, 3], [false, true, false], :r)
["1.0+", "2.0", "3.0+"]

julia> srv2 = Surv([4, 5, 6], [false, true, true], :r)
["4.0+", "5.0", "6.0"]

julia> merge(srv1, srv2) # OneSidedSurv
["1.0+", "2.0", "3.0+", "4.0+", "5.0", "6.0"]

julia> srv1 = Surv([1, 2, 3], [4, 5, 6])
[(1.0, 4.0), (2.0, 5.0), (3.0, 6.0)]

julia> srv2 = Surv([7, 8, 9], [10, 11, 12])
[(7.0, 10.0), (8.0, 11.0), (9.0, 12.0)]

julia> merge(srv1, srv2) # TwoSidedSurv
[(1.0, 4.0), (2.0, 5.0), (3.0, 6.0), (7.0, 10.0), (8.0, 11.0), (9.0, 12.0)]
```
"""
function Base.merge(A::OneSidedSurv...)
    length(unique(map(typeof, A))) > 1 &&
        throw(ArgumentError("Objects to merge must either be all `RCSurv` or `LCSurv`"))
    T = zeros(0)
    Δ = BitVector()
    foreach(srv -> push!(T, srv.time...), A)
    foreach(srv -> push!(Δ, srv.status...), A)
    return typeof(A[1])(T, Δ)
end

function Base.merge(A::TwoSidedSurv...)
    start = zeros(0)
    stop = zeros(0)
    foreach(srv -> push!(start, srv.start...), A)
    foreach(srv -> push!(stop, srv.stop...), A)
    return Surv(start, stop)
end

"""
    reverse(srv::OneSidedSurv)

Reverse survival outcome (not in-place). Useful for computing non-parametric estimators of the censoring distribution instead of the survival distribution.

# Examples
```jldoctest
julia> srv = Surv([1, 2, 3, 4, 5, 6], [false, true, false, false, true, true], :r)
["1.0+", "2.0", "3.0+", "4.0+", "5.0", "6.0"]

julia> reverse(srv)
["1.0", "2.0+", "3.0", "4.0", "5.0+", "6.0+"]
```
"""
Base.reverse(srv::OneSidedSurv) = typeof(srv)(srv.time, .!srv.status)
