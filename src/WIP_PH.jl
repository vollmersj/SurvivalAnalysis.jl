using Distributions
using LineSearches, Optim, NLSolversBase
using LinearAlgebra: diag
using Survival
using DataFrames
using RCall

n = 100;
Xs = Uniform();
X = DataFrame("X1" => rand(Xs, n), "X2" => rand(Xs, n),
        "X3" => rand(Xs, n));
T = X.X1 * 2 + X.X2 * -0.1 + X.X3 * 3.1;
Δ = rand(Binomial(), n) .== 1;
Y = Surv(T, Δ, "right");

struct FittedParametricPH
    coefficients::NamedTuple{(:betas, :scale), Tuple{Vector{Float64}, Float64}}
    hessian::Matrix
    var_cov::Matrix
    t::Vector
    baseline::ContinuousUnivariateDistribution
    routine
end

function fit_PH(X::DataFrame, Y::Survival.rcSurv, d::String, init::Number = 1)
    @assert d in ["Weibull", "Exponential"]

    # θ[1] = scale, θ[2] = β₀, θ[...] = β...
    nβ = ncol(X)
    npar = nβ + 2
    X = Matrix(X)

    init = [init, 0, zeros(nβ)...]

    function loglik(x, t, δ, ϕ, β₀, β)
        ϕ <= 0 && return Inf # reject impossible candidates
        if d == "Exponential"
            l = (δ .* (log(ϕ) .+ β₀ .+ x*β)) .- (exp.(β₀ .+ x*β) .* ϕ .* t)
        else
            l = (δ .* (log(ϕ) .+ ((ϕ-1) .* log.(t)) .+ β₀ .+ x*β)) .- (exp.(β₀ .+ x*β) .* t.^ϕ)
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

#    ζ̂ = getfield(Distributions, Symbol(d))(ϕ̂...)
    ζ̂ = Exponential() # dummy argument
    FittedParametricPH((betas = β̂, scale = ϕ̂), H, V, t, ζ̂, opt)
end

function predict_PH(fit::FittedParametricPH, X::DataFrame)
    η = exp.(Matrix(X) * fit.coefficients.betas)
    ζ = ParametricPH.(fit.baseline, η)
    (lp = η, distr = ζ)
end


e = fit_PH(X, Y, "Exponential");

R"
library(survival)
surv = Surv($T, $Δ)
survreg(surv ~ ., data = $X, dist = 'exponential')
"

w = fit_PH(X, Y, "Weibull");
R"
library(survival)
surv = Surv($T, $Δ)
survreg(surv ~ ., data = $X, dist = 'weibull')
"

R"
library(eha)
w = weibreg(surv ~ ., data = $X)
w$coefficients / exp(w$shape)
"
