mutable struct ContinuousPHDistribution <: ContinuousUnivariateDistribution
    ζ::ContinuousUnivariateDistribution ## Baseline
    η::Float64 # linear predictor
end

ContinuousPHDistribution(d) = ContinuousPHDistribution(d, 0.0)

#@distr_support ContinuousPHDistribution 0.0 Inf

Distributions.params(d::ContinuousPHDistribution) = (d.ζ, d.η)
Distributions.partype(::ContinuousPHDistribution) = Float64

#### Evaluation

Distributions.pdf(d::ContinuousPHDistribution, x::Real) = hazard(d.ζ, x) * exp(d.η) * (ccdf(d.ζ, x)^exp(d.η))
Distributions.cdf(d::ContinuousPHDistribution, x::Real) = 1 - (ccdf(d.ζ, x)^exp(d.η))
Distributions.ccdf(d::ContinuousPHDistribution, x::Real) = ccdf(d.ζ, x)^exp(d.η)
hazard(d::ContinuousPHDistribution, x::Real) = hazard(d.ζ, x) * exp(d.η)

Base.show(io::IO, d::ContinuousAFTDistribution) = print(io, "PH(ζ=$(d.ζ), η=$(d.η))")
