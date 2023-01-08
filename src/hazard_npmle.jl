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
    converged::Bool
end

# Build design matrices D and Dstar which respectively indicate which
# truncation intervals and which observation intervals contain each
# support interval.
function design(spt, Y)

    start = Y.start
    stop = Y.stop
    ltrunc = Y.ltrunc

    n = length(Y)
    p = length(spt)
    D = spzeros(n, p - 1)
    Dstar = spzeros(n, p - 1)

    for i = 1:n
        itrunc = Interval{:Closed,:Open}(ltrunc[i], start[i])
        iobs = Interval{:Closed,:Open}(start[i], stop[i])
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
    spt = Float64[]
    for i in eachindex(Y.start)
        push!(spt, Y.start[i])
        push!(spt, Y.stop[i])
        push!(spt, Y.ltrunc[i])
    end
    spt = unique(spt)
    sort!(spt)

    return spt
end

function HazardNPMLE(Y::IntSurv)
    spt = get_support(Y)
    D, Dstar = design(spt, Y)
    return HazardNPMLE(Y, spt, diff(spt), D, Dstar, [], [], [], Symbol(), false)
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
    D = ms.D

    n = length(Y)
    p = length(support)
    dpar = cumsum(exp.(par))
    ddpar = dpar .* duration

    lp = D * ddpar
    lps = Dstar * ddpar

    wgt = length(ms.Y.weight) == 0 ? ones(n) : ms.Y.weight

    ll = dot(wgt, log.(1 .- exp.(-lps)) - lp)
    return ll
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
    D = ms.D

    length(G) == length(par) || throw(ArgumentError("G and par must have the same length"))

    n = length(Y)
    p = length(support)
    dpar = cumsum(exp.(par))
    ddpar = dpar .* duration

    jac = [i >= j ? 1.0 : 0.0 for i = 1:p-1, j = 1:p-1]

    lp = D * ddpar
    lps = Dstar * ddpar

    wgt = length(ms.Y.weight) == 0 ? ones(n) : ms.Y.weight

    G .= 0
    ee = exp.(-lps)
    G .= Dstar' * ((wgt .* ee) ./ (1 .- ee)) - D' * wgt
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

    ms = HazardNPMLE(Y)
    return _fit_nondecreasing_hazard(ms, optim_opts)
end

function _fit_nondecreasing_hazard(ms::HazardNPMLE, opts)
    Y = ms.Y
    support = ms.support
    duration = ms.duration
    Dstar = ms.Dstar
    D = ms.D

    p = length(support)

    # Starting values
    par = ones(p - 1) / (p - 1)

    f = x -> -loglike(ms, x)
    g! = (G, x) -> begin
        score!(ms, G, x)
        G .*= -1
    end

    # More iterations than the default is needed since
    # the optimization is high-dimensional.
    rr = Optim.optimize(f, g!, par, LBFGS(), opts)
    ms.converged = Optim.converged(rr)
    if !ms.converged
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
