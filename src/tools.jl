## tools = exported utils
"""
    hazard(d::UnivariateDistribution, x::Real)
    hₜ(d::UnivariateDistribution, x::Real)

Compute the hazard function of distribution `d` at point `x`.

The hazard function for random variable ``t`` is defined as

```math
h_t(x) = \\frac{f_t(x)}{S_t(x)}
```

where ``fₜ`` is the pdf of ``t`` and ``Sₜ`` is the survival function of ``t``.

# Examples
```jldoctest
julia> using Distributions

julia> hazard(Binomial(5, 0.5), 3)
1.6666666666666679
```
"""
hazard(d::UnivariateDistribution, x::Real) = pdf(d, x) / survival(d, x)

"""
    cum_hazard(d::UnivariateDistribution, x::Real)
    Hₜ(d::UnivariateDistribution, x::Real)

Compute the cumulative hazard function of distribution `d` at point `x`.

The cumulative hazard function for random variable ``t`` is defined as

```math
H_t(x) = \\int^x_0 h_t(u) du = -log(S_t(x))
```

where ``hₜ`` is the hazard function of ``t`` and ``Sₜ`` is the survival function of ``t``.

# Examples
```jldoctest
julia> using Distributions

julia> cum_hazard(Binomial(5, 0.5), 3)
1.6739764335716711
```
"""
cum_hazard(d::UnivariateDistribution, x::Real) = -log(survival(d, x))

# aliases auto-documented
const Fₜ = cdf
const fₜ = pdf
const pₜ = pdf
const survival = ccdf
const Sₜ = ccdf
const Hₜ = cum_hazard
const hₜ = hazard
