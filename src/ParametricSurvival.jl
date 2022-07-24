abstract type ParametricSurvival <: SurvivalModel end

StatsBase.coef(obj::ParametricSurvival) = obj.coefficients
Distributions.scale(obj::ParametricSurvival) = scale(obj.distribution)
Distributions.shape(obj::ParametricSurvival) = shape(obj.distribution)

function StatsBase.fit(t::Type{<:ParametricSurvival}, X::AbstractMatrix{<:Real}, Y::RCSurv,
    d::Type{T}; init::Number = 1) where {T <: ContinuousUnivariateDistribution}
    @assert d in [Weibull, Exponential]
    # crude method to force intercept
    X = X[:,1] == ones(size(X, 1)) ? X : hcat(ones(size(X, 1)), X)
    fit!(t(d, X), X, Y, init)
end

function _fitParametricSurvival(obj, X, Y, init, llik, dtrafo)
    # θ[1] = scale, θ[2:end] = βs
    nβ = size(X, 2)
    npar = nβ + 1
    X = Matrix(X)

    init = [init, zeros(nβ)...]

    func = TwiceDifferentiable(
        θ -> θ[1] <= 0 ? Inf : -∑(llik(obj.baseline, X, Y.time, Y.status, θ[1], θ[2:end])),
        init; autodiff=:forward);

    opt = optimize(func, init)
    θ̂ = Optim.minimizer(opt)
    obj.routine = opt

    obj.hessian = hessian!(func, θ̂)
    ## FIXME - NEED TO FIGURE OUT WHY INVERSION SOMETIMES FAILS
    obj.var_cov = try
        inv(obj.hessian)
    catch
        fill(NaN, (npar, npar))
    end
    obj.scale = θ̂[1]

    obj.coefficients = θ̂[2:end]

    ## FIXME - NEED TO FIGURE OUT WHY VAR_COV SOMETIMES < 0
    obj.tstats = try
        obj.coefficients./sqrt.(diag(obj.var_cov)[2:end])
    catch
        fill(NaN, nβ)
    end

    obj.baseline = dtrafo(obj.baseline, obj.coefficients, obj.scale)

    obj
end
