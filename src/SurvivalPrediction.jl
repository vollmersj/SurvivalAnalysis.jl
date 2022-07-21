abstract type SurvivalPrediction end

struct DeterministicSurvivalPrediction{T<:Float64} <: SurvivalPrediction
    lp::Vector{T}
    crank::Vector{T}
    survivaltime::Vector{T}
end

struct DiscreteSurvivalPrediction{T<:Float64} <: SurvivalPrediction
    distr::Vector{DiscreteNonParametric}
    lp::Vector{T}
    crank::Vector{T}
    survivaltime::Vector{T}
    times::Vector{T}
    mat::Matrix{T}
end

struct ContinuousSurvivalPrediction{T<:Float64} <: SurvivalPrediction
    distr::Vector{ContinuousUnivariateDistribution}
    lp::Vector{T}
    crank::Vector{T}
    survivaltime::Vector{T}
end

function _survPredict(;ζ::Vector{Distribution}, η::Vector{T}, ϕ::Vector{T},
                        T̂::Vector{T}, Ts::Vector{T}, Ŝ::Matrix{T}) where {T<:Float64}
    n = []
    !missing(ζ) && push!(n, length(ζ))
    !missing(η) && push!(n, length(η))
    !missing(ϕ) && push!(n, length(ϕ))
    !missing(T̂) && push!(n, length(T̂))
    !missing(Ŝ) && push!(n, size(Ŝ, 1))
    n = unique(n)

    @assert length n == 1 error("Supplied parameters of different lengths")

    # construct distribution from matrix if available
    if ismissing(Ts) + ismissing(Ŝ) == 1
        error("Either both 'Ts' should be provided 'Ŝ' or neither")
    elseif !ismissing(Ts) && !ismissing(Ŝ) && ismissing(ζ)
        @assert length(Ts) == size(Ŝ, 2)
        if Ts[1] != 0
            # Set S(0) = 1
            Ŝ = hcat(ones(n), Ŝ)
            Ts = [0, Ts...]
        end
        # calculate pmf and create distribution
        p = mapreduce(x -> abs.(diff(x)), hcat, eachrow(Ŝ))'
        ζ = DiscreteNonParametric(Ts, p, check_args = false)
    end

    # no transformations assumed
    η = ismissing(η) ? fill(NaN, n) : η
    ϕ = ismissing(ϕ) ? fill(NaN, n) : ϕ
    T̂ = ismissing(T̂) ? fill(NaN, n) : T̂

    if ismissing(ζ)
        DeterministicSurvivalPrediction(η, ϕ, T̂)
    elseif ζ[1] isa ContinuousUnivariateDistribution
        ContinuousSurvivalPrediction(ζ, η, ϕ, T̂)
    elseif ζ[1] isa DiscreteNonParametric
        if ismissing(Ts)
            Ts = unique(support.(ζ))
            @assert length(Ts) == 1 "Predicted distributions (ζ) must have same survival times"
            Ŝ = mapreduce(s -> 1 .- cumsum(probs(s)), hcat, ζ)'
        end
        DiscreteSurvivalPrediction(ζ, η, ϕ, T̂, Ts, Ŝ)
    end
end
