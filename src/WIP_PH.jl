using Distributions
using Optim, NLSolversBase, Random
using LinearAlgebra: diag
using Survival

function fit_PH(d::String, init::Vector{T} = [0.1]) where {T <: Real}
    @assert d in ["Weibull", "Exponential"]

    init = d == "Weibull" ? ones(2) * 0.1 : ones(1) * 0.1

    function loglik(T, Δ, θ)
        any(i -> i <= 0, θ) && return Inf # reject impossible candidates
        ζ = getfield(Distributions, Symbol(d))(θ...)
        -Survival.∑(Δ .* log.(Survival.hₜ.(ζ, T)) .- Survival.Hₜ.(ζ, T))
    end

    func = TwiceDifferentiable(
        θ -> loglik(surv.time, surv.status, θ), init; autodiff=:forward);

    opt = optimize(func, init)
    θ̂ = Optim.minimizer(opt)
    H = hessian!(func, θ̂)
    Var = inv(H)
    t = θ̂./sqrt.(diag(Var))

    (θ̂ = θ̂, H = H, Var = Var, t = t, routine = opt)
end

n = 100
T = round.(rand(Uniform(1, 10), n));
Δ = rand(Binomial(), n) .== 1;
surv = Survival.Surv(T, Δ, "right");

w = fit_PH("Weibull");
e = fit_PH("Exponential");
