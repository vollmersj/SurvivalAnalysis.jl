abstract type OneSidedSurv end
abstract type TwoSidedSurv end

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

struct OutcomeTimes
    times::Vector{OneSidedSurv}
end

outcomeTimes(ot::OutcomeTimes) = map(x -> x.time, ot.times)
outcomeStatus(ot::OutcomeTimes) = map(x -> x.status, ot.times)
uniqueTimes(ot::OutcomeTimes) = unique(map(x -> x.time, ot.times))

outcomeTimes(ot)
outcomeStatus(ot)
uniqueTimes(ot)
