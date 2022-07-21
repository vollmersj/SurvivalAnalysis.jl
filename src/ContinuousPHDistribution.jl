mutable struct ContinuousPHDistribution{T} <: ContinuousUnivariateDistribution where {T <: ContinuousUnivariateDistribution}
    ζ::T ## Baseline
    η::Float64 # linear predictor
end

#@distr_support ContinuousPHDistribution 0.0 Inf

params(d::ContinuousPHDistribution) = ()
partype(::ContinuousPHDistribution) = Float64

#### Evaluation

Distributions.pdf(d::ContinuousPHDistribution, x::Real) = hazard(d.ζ, x) * exp(d.η) * (ccdf(d.ζ, x)^exp(d.η))
Distributions.cdf(d::ContinuousPHDistribution, x::Real) = 1 - (ccdf(d.ζ, x)^exp(d.η))
Distributions.ccdf(d::ContinuousPHDistribution, x::Real) = ccdf(d.ζ, x)^exp(d.η)
hazard(d::ContinuousPHDistribution, x::Real) = hazard(d.ζ, x) * exp(d.η)
