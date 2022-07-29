struct ContinuousAFTDistribution <: ContinuousUnivariateDistribution
    ζ::ContinuousUnivariateDistribution ## Baseline
    η::Float64 # linear predictor
end

#@distr_support ContinuousAFTDistribution 0.0 Inf

Distributions.params(d::ContinuousAFTDistribution) = (d.ζ, d.η)
Distributions.partype(::ContinuousAFTDistribution) = Float64

#### Evaluation

Distributions.pdf(d::ContinuousAFTDistribution, x::Real) =
    exp(-d.η) * hazard(d.ζ, x / exp(d.η)) * ccdf(d.ζ, x / exp(d.η))
Distributions.cdf(d::ContinuousAFTDistribution, x::Real) = cdf(d.ζ, x / exp(d.η))
Distributions.quantile(d::ContinuousAFTDistribution, q::Real) = quantile(d.ζ, q) * exp(d.η)
Distributions.ccdf(d::ContinuousAFTDistribution, x::Real) = ccdf(d.ζ, x / exp(d.η))
hazard(d::ContinuousAFTDistribution, x::Real) = exp(-d.η) * hazard(d.ζ, x / exp(d.η))

Base.show(io::IO, d::ContinuousAFTDistribution) = print(io, "AFT(ζ=$(d.ζ), η=$(d.η))")
