using Distributions
using Optim, NLSolversBase, Random
using LinearAlgebra: diag
using Survival
using DataFrames

# L = ∏ λ(t|x)^δ S(t|x)
# l = ∑(Δ .* log.(hₜ.(ζ, T).*exp.(x'β)) .- exp.(x'β) * Hₜ.(ζ, T))

struct FittedParametricPH
    coefficients::NamedTuple{(:betas, :pars), NTuple{2, Vector{Float64}}}
    hessian::Matrix
    var_cov::Matrix
    t::Vector
    baseline::ContinuousUnivariateDistribution
    routine
end

function fit_PH(X::DataFrame, Y::Survival.rcSurv, d::String, init::Number = 0.1)
    @assert d in ["Weibull", "Exponential"]

    nβ = ncol(X)
    nϕ = d == "Weibull" ? 2 : 1
    npar = nβ + nϕ
    X = Matrix(X)

    init = [zeros(nβ)..., (ones(nϕ) .* init)...]

    function loglik(x, t, δ, β, ϕ)
        any(i -> i <= 0, ϕ) && return Inf # reject impossible candidates
        ζ = getfield(Distributions, Symbol(d))(ϕ...)
        -Survival.∑(δ .* log.(Survival.hₜ.(ζ, t).*exp.(x*β)) .- exp.(x*β) .* Survival.Hₜ.(ζ, t))
    end

    func = TwiceDifferentiable(
        θ -> loglik(X, Y.time, Y.status, θ[1:nβ], θ[nβ+1:npar]),
        init; autodiff=:forward);

    opt = optimize(func, init)
    θ̂ = Optim.minimizer(opt)
    H = hessian!(func, θ̂)
    V = inv(H)
    β̂ = θ̂[1:nβ]
    t = β̂./sqrt.(diag(V)[1:nβ])
    ϕ̂ = θ̂[nβ+1:npar]

    ζ̂ = getfield(Distributions, Symbol(d))(ϕ̂...)

    FittedParametricPH((betas = β̂, pars = ϕ̂), H, V, t, ζ̂, opt)
end

function predict_PH(fit::FittedParametricPH, X::DataFrame)
    η = exp.(Matrix(X) * fit.coefficients.betas)
    ζ = ParametricPH.(fit.baseline, η)
    (lp = η, distr = ζ)
end

n = 100;
Xs = Uniform();
X = DataFrame("X1" => rand(Xs, n), "X2" => rand(Xs, n),
        "X3" => rand(Xs, n));
T = X.X1 * 2 + X.X2 * -0.1 + X.X3 * 3.1;
Δ = rand(Binomial(), n) .== 1;
Y = Surv(T, Δ, "right");

#w = fit_PH(X, Y, "Weibull");
e = fit_PH(X, Y, "Exponential");
p = predict_PH(e, X);
hazard.(p.dist, 1)

# using RCall

# R"
#     library(survival)
#     time = $T
#     status = $Δ
#     df = $X
#     surv = Surv(time, status)
#     survreg(surv ~ ., data = df, dist = 'exponential')
# "
