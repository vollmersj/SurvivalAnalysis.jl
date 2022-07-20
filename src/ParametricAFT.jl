struct ParametricAFT <: ContinuousUnivariateDistribution
    ζ::ContinuousUnivariateDistribution ## Baseline
    η::Float64 # linear predictor
end

#@distr_support ParametricAFT 0.0 Inf

params(d::ParametricAFT) = ()
partype(::ParametricAFT) = Float64

#### Evaluation

# need to confirm sign of all the below
Distributions.pdf(d::ParametricAFT, x::Real) = exp(-d.η) * hazard(d.ζ, x / exp(d.η)) * ccdf(d.ζ, x / exp(d.η))
Distributions.cdf(d::ParametricAFT, x::Real) = cdf(d.ζ, x / exp(d.η))
Distributions.quantile(d::ParametricAFT, q::Real) = quantile(d.ζ, q) * exp(d.η)
Distributions.ccdf(d::ParametricAFT, x::Real) = ccdf(d.ζ, x / exp(d.η))
hazard(d::ParametricAFT, x::Real) = exp(-d.η) * hazard(d.ζ, x / exp(d.η))
