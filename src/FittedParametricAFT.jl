struct FittedParametricAFT <: FittedParametric
    coefficients::NamedTuple{(:betas, :scale), Tuple{Vector{Float64}, Float64}}
    hessian::Matrix
    var_cov::Matrix
    t::Vector
    baseline::ContinuousUnivariateDistribution
    routine
end

function fit_AFT(X::DataFrame, Y::Survival.rcSurv, d::String, init::Number = 1)
    @assert d in ["Weibull", "Exponential"]

    # θ[1] = scale, θ[2] = β₀, θ[...] = β...
    nβ = ncol(X)
    npar = nβ + 2
    X = Matrix(X)

    init = [init, 0, zeros(nβ)...]

    function loglik(x, t, δ, ϕ, β₀, β)
        ϕ <= 0 && return Inf # reject impossible candidates
        if d == "Exponential"
            l = (δ .* (-β₀ .- x*β)) .- (exp.(-β₀ .- x*β) .* t)
        elseif d == "Weibull"
            l = (δ .* (log(1/ϕ) .+ (((1/ϕ)-1) .* log.(t)) .- (β₀ .+ x*β) ./ ϕ)) .- (exp.((-β₀ .- x*β) ./ ϕ) .* t.^(1/ϕ))
        end
        -Survival.∑(l)
    end

    func = TwiceDifferentiable(
        θ -> loglik(X, Y.time, Y.status, θ[1], θ[2], θ[3:length(θ)]),
        init; autodiff=:forward);

    opt = optimize(func, init)
    θ̂ = Optim.minimizer(opt)
    nθ = length(θ̂)
    H = hessian!(func, θ̂)
    V = try
        inv(H)
    catch
        fill(NaN, (npar, npar))
    end
    ϕ̂ = θ̂[1]
    β̂ = θ̂[2:nθ]
    t = β̂./sqrt.(diag(V)[2:nθ])

    ζ̂ = d == "Exponential" ? Exponential(exp(β̂[1])) : Weibull(1/ϕ̂, exp(β̂[1]))
    FittedParametricAFT((betas = β̂, scale = ϕ̂), H, V, t, ζ̂, opt)
end

function predict_Parametric(fit::FittedParametricAFT, X::DataFrame)
    η = Matrix(X) * fit.coefficients.betas[2:end]
    ζ = ContinuousAFTDistribution.(fit.baseline, η)
    _survPredict(ζ = ζ, η = η, ϕ = -η)
end
