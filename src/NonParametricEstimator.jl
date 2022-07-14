abstract type NonParametricEstimator end

function fit_NPE(class, Surv::rcSurv, point_est::Function, var_est::Function,
    surv_trafo::Function, sd_trafo::Function)
    stats = survStats(Surv, events_only = true)
    n = length(stats.time)
    p = zeros(n)
    v = zeros(n)
    for (i, t) in enumerate(stats.time)
    p[i] = point_est(stats.nevents[i], stats.nrisk[i])
    v[i] = var_est(stats.nevents[i], stats.nrisk[i])
    end
    # for both NPEs variance calculated as plug-in
    variance = cumsum(v)
    surv = surv_trafo(p)
    sd = sd_trafo(variance, surv)
    # Set S(0) = 1
    T = [0, stats.time...]
    Sₜ = [1, surv...]
    # calculate pmf and create distribution
    pₓ = [0, abs.(diff(Sₜ))...]
    d = DiscreteNonParametric(T,  pₓ, check_args = false)
    # return class
    class(T, Sₜ, [0, sd...], d)
end

# Calculates confidence intervals using the method of Kalbfleisch and Prentice (1980)
function confint_NPE(npe, t, α, trafo)
    q = quantile(Normal(), 1 - α/2)
    which = searchsortedlast(npe.time, t)
    trafo(npe, q, which)
end

# syntactic sugar for vectorising over all times
function StatsBase.confint(npe::NonParametricEstimator; α::Float64 = 0.05)
    confint.(Ref(npe), npe.time, α = α)
end

@recipe function f(npe::NonParametricEstimator, plot_confint::Bool = true, α::Float64 = 0.05)
    seriestype := :steppost
    ylims := (0, 1)
    legend := false
    @series begin
        linecolor   --> :black
        npe.time, npe.survival
    end
    if plot_confint
        linecolor   --> :blue
        cis = confint.(Ref(npe), npe.time, α = α)
        lb = map(x -> x[1], cis)
        ub = map(x -> x[2], cis)
        @series begin
            npe.time, lb
        end
        @series begin
            npe.time, ub
        end
    end
end
