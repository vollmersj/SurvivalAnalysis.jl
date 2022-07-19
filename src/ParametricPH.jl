struct ParametricPH <: ContinuousUnivariateDistribution
    ζ::ContinuousUnivariateDistribution ## Baseline
    η::Float64 # linear predictor
end

#@distr_support ParametricPH 0.0 Inf

params(d::ParametricPH) = ()
partype(::ParametricPH) = Float64

#### Evaluation

Distributions.pdf(d::ParametricPH, x::Real) = hazard(d, x) * ccdf(d, x)
Distributions.ccdf(d::ParametricPH, x::Real) = ccdf(d.ζ, x)^exp(d.η)
hazard(d::ParametricPH, x::Real) = hazard(d.ζ, x)*exp(d.η)
