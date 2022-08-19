mutable struct ContinuousPHDistribution <: ContinuousUnivariateDistribution
    distr::ContinuousUnivariateDistribution ## Baseline
    lp::Float64 # linear predictor
end

#@distr_support ContinuousPHDistribution 0.0 Inf

Distributions.params(d::ContinuousPHDistribution) = (d.distr, d.lp)
Distributions.partype(::ContinuousPHDistribution) = Float64

#### Evaluation

Distributions.pdf(d::ContinuousPHDistribution, x::Real) =
    hazard(d.distr, x) * exp(d.lp) * (ccdf(d.distr, x)^exp(d.lp))
Distributions.cdf(d::ContinuousPHDistribution, x::Real) = 1 - (ccdf(d.distr, x)^exp(d.lp))
Distributions.quantile(d::ContinuousPHDistribution, q::Real) =
    quantile(d.distr, (1 - exp(log(1 - q) / exp(d.lp))))
Distributions.ccdf(d::ContinuousPHDistribution, x::Real) = ccdf(d.distr, x)^exp(d.lp)
hazard(d::ContinuousPHDistribution, x::Real) = hazard(d.distr, x) * exp(d.lp)

Base.show(io::IO, d::ContinuousPHDistribution) =
    print(io, "PH(distr=$(d.distr), lp=$(d.lp))")
