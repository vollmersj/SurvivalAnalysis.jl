"""
     HazardNPMLE <: SurvivalEstimator

See [`SurvivalEstimator`](@ref).
"""
mutable struct HazardNPMLE <: SurvivalEstimator

    # The data
    Y::IntSurv

    # The NPMLE takes the form of a piecewise-constant hazard function
    # with support on the intervals [support[i], support[i+1]).
    support::Vector{Float64}

    # The duration of each interval in `support`.  The length of
    # `duration` is one less than the length of `support`.
    duration::Vector{Float64}

    # A n x p matrix (n=sample size, p=number of support intervals),
    # such that D[i,j] = 1 iff the j^th support interval is contained
    # inside the interval [ltrunc[i], start[i]).
    D::AbstractMatrix

    # A n x p matrix (n=sample size, p=number of support intervals),
    # such that Dstar[i,j] = 1 iff the j^th support interval is
    # contained inside the interval [start[i], stop[i]).
    Dstar::AbstractMatrix

    # The untransformed parameter vector, directly from the optimizer.
    par::Vector{Float64}

    # The estimated hazard values on the intervals defined by support.
    hazard::Vector{Float64}

    # The estimated cumulative hazard function.
    cum_hazard::Vector{Float64}

    # The shape constraint imposed on the hazard function estimate.
    shape::Symbol

    # An indicator that the NPMLE converged.
    optim_rslt::Any
end

# Build design matrices D and Dstar which respectively indicate which
# truncation intervals and which observation intervals contain each
# support interval.
function design(spt, Y)

    start = Y.start
    stop = Y.stop
    ltrunc = length(Y.ltrunc) > 0 ? Y.ltrunc : spzeros(length(start))

    n = length(Y)
    p = length(spt)
    D = spzeros(n, p - 1)
    Dstar = spzeros(n, p - 1)
    itrunc0 = Interval{:Closed,:Open}(0.0, 0.0)

    for i = 1:n

        # A truncation interval
        itrunc = Interval{:Closed,:Open}(ltrunc[i], start[i])

        # An observation interval
        iobs = Interval{:Closed,:Open}(start[i], stop[i])

        # Loop over the support intervals
        for j = 1:p-1
            ivl = Interval{:Closed,:Open}(spt[j], spt[j+1])
            if ivl ⊆ itrunc
                D[i, j] = 1
            end
            if ivl ⊆ iobs
                Dstar[i, j] = 1
            end
        end
    end

    return D, Dstar
end

# Calculate the support intervals for the NPMLE.
function get_support(Y)
    spt = vcat(Y.start, Y.stop, Y.ltrunc)
    spt = unique(spt)

    sort!(spt)
    return spt
end

function HazardNPMLE(Y::IntSurv)
    spt = get_support(Y)
    D, Dstar = design(spt, Y)
    return HazardNPMLE(Y, spt, diff(spt), D, Dstar, [], [], [], Symbol(), nothing)
end

"""
    loglike(ms::HazardNPMLE, par::Vector{Float64})

Evaluate and return the nonparametric log-likelihood function for the
data in `ms` at the parameter vector `par`.
"""
function loglike(ms::HazardNPMLE, par::Vector{Float64})

    Y = ms.Y
    support = ms.support
    duration = ms.duration
    Dstar = ms.Dstar

    n = length(Y)
    p = length(support)
    dpar = cumsum(exp.(par))
    ddpar = dpar .* duration
    lps = Dstar * ddpar

    wgt = length(Y.weight) == 0 ? ones(n) : Y.weight

    v = log.(1 .- exp.(-lps)) - ms.D * ddpar
    return dot(wgt, v)
end

"""
    score!(ms::HazardNPMLE, G::Vector{Float64}, par::Vector{Float64})

Calculate the score function for the nonparametric log-likelihood
function for the data in `ms` at the parameter vector `par`.  The
score vector is copied into the array `G` which must have the same
length as `par`.
"""
function score!(ms::HazardNPMLE, G::Vector{Float64}, par::Vector{Float64})

    Y = ms.Y
    support = ms.support
    duration = ms.duration
    Dstar = ms.Dstar

    length(G) == length(par) || throw(ArgumentError("G and par must have the same length"))

    n = length(Y)
    p = length(support)
    dpar = cumsum(exp.(par))
    ddpar = dpar .* duration

    jac = [i >= j ? 1.0 : 0.0 for i = 1:p-1, j = 1:p-1]

    lps = Dstar * ddpar

    wgt = length(Y.weight) == 0 ? ones(n) : Y.weight

    G .= 0
    ee = exp.(-lps)
    G .= Dstar' * ((wgt .* ee) ./ (1 .- ee)) - ms.D' * wgt
    G .*= duration
    G .= Diagonal(exp.(par)) * jac' * G
end

"""
     fit(::Type{HazardNPMLE}, Y::IntSurv; shape=:nondecreasing_hazard)

Estimate a distribution based on interval-censored data subjet to left
truncation.  It is known that the nonparametric MLE may not exist or
may fail to be unique unless constraints are imposed.  Currently the
only impemented option is to constrain the hazard function to be
non-decreasing, in which case the approach of Pan and Chappell
(Biometrics, 1998) is used.
"""
function StatsBase.fit(::Type{HazardNPMLE}, Y::IntSurv; shape = :nondecreasing_hazard,
                       optim_opts=Optim.Options(iterations = 100000))

    shape == :nondecreasing_hazard ||
        throw(ArgumentError("Only shape=:nondecreasing_hazard is currently allowed"))

    if any(Y.start .== Y.stop)
        throw(ArgumentError("Observations may not have identical start and stop times"))
    end

    ms = HazardNPMLE(Y)
    return _fit_nondecreasing_hazard(ms, optim_opts)
end

function _optfuns(ms)

    # The MLE may have some parameters equal to -infinity
    # corresponding to flat regions of the hazard function.
    # To make the optimization have a finite solution, introduce
    # a barrier function that blows up at negative infinity.
    ee = 1e-6

    f = x -> -loglike(ms, x) + ee*sum(exp.(-x))

    g! = (G, x) -> begin
        score!(ms, G, x)
        G .*= -1
        G .-= ee*exp.(-x)
    end

    return f, g!
end

function _fit_nondecreasing_hazard(ms::HazardNPMLE, opts)
    Y = ms.Y
    support = ms.support
    duration = ms.duration
    Dstar = ms.Dstar
    D = ms.D

    p = length(support)

    # Starting values
    par = -5*ones(p - 1)

    f, g! = _optfuns(ms)

    # More iterations than the default is needed since
    # the optimization is high-dimensional.
    rr = Optim.optimize(f, g!, par, LBFGS(), opts)
    ms.optim_rslt = rr
    if !Optim.converged(rr)
        @warn("monotone hazard estimation did not converge")
    end
    ms.par = Optim.minimizer(rr)
    ms.hazard = cumsum(exp.(ms.par))
    ms.cum_hazard = cumsum(ms.duration .* ms.hazard)

    return ms
end

"""
     hazard(m::HazardNPMLE, x::Real)

Return the estimated hazard function at `x`.
"""
function hazard(ms::HazardNPMLE, x::Real)
    ii = searchsortedfirst(ms.support[1:end-1], x) - 1
    return ii <= 0 ? 0.0 : ms.hazard[ii]
end

"""
     cum_hazard(m::HazardNPMLE, x::Real)

Return the estimated cumulative hazard function at `x`.
"""
function cum_hazard(ms::HazardNPMLE, x::Real)
    ii = searchsortedfirst(ms.support[1:end-1], x) - 1
    return ii <= 0 ? 0.0 : ms.cum_hazard[ii]
end

"""
     survival(m::HazardNPMLE, x::Real)

Return the estimated survival function at `x`.
"""
function survival(ms::HazardNPMLE, x::Real)
    ii = searchsortedfirst(ms.support[1:end-1], x) - 1
    return ii <= 0 ? 1.0 : exp(-ms.cum_hazard[ii])
end
