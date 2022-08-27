#-----------------
# exported tools
#-----------------
"""
    hazard(d::UnivariateDistribution, x::Real)

    Aliases: hₜ

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

    Aliases: Hₜ

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

#-------------------
# non-exported tools
#-------------------
∓(x, y) = (x + y, x - y)
∑(A) = sum(A)
∑(A...) = sum(A)
∏(A) = prod(A)
∏(A...) = prod(A)
test_proportion(x) = x <= 1 && x>= 0

function pstring(s, p)
    s = p < 0 ? string("/", s) : s
    str = string(abs(p))

    str == "0" && return "1"
    str == "1" && return s

    str = replace(str, "1" => Char('\U000B9'))
    str = replace(str, "2" => Char('\U000B2'))
    str = replace(str, "3" => Char('\U000B3'))
    str = replace(str, "4" => Char('\U02074'))
    str = replace(str, "5" => Char('\U02075'))
    str = replace(str, "6" => Char('\U02076'))
    str = replace(str, "7" => Char('\U02077'))
    str = replace(str, "8" => Char('\U02078'))
    str = replace(str, "9" => Char('\U02079'))
    str = replace(str, "0" => Char('\U02070'))

    return string(s, str)
end

function c_pstring(s1, s2)
    s1 == "1" && s2 == "1" && return "1"
    # cleanup
    if first(s1) == '/' && first(s2) != '/'
        tmp = s2
        s2 = s1
        s1 = tmp
    end

    return first(s1) == '/' ? string("1", s1, s2) : string(s1, s2)
end
