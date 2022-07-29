abstract type SurvivalPrediction end

struct DeterministicSurvivalPrediction{T<:Float64} <: SurvivalPrediction
    lp::Vector{T}
    crank::Vector{T}
    time::Vector{T}
end

struct DiscreteSurvivalPrediction{T<:Float64} <: SurvivalPrediction
    distr::Vector{DiscreteNonParametric}
    lp::Vector{T}
    crank::Vector{T}
    time::Vector{T}
    survival_matrix::NamedTuple{(:time, :surv), Tuple{Vector{T}, Matrix{T}}}
end

struct ContinuousSurvivalPrediction{T<:Float64} <: SurvivalPrediction
    distr::Vector{<:ContinuousUnivariateDistribution}
    lp::Vector{T}
    crank::Vector{T}
    time::Vector{T}
end

function SurvivalPrediction(;ζ::Vector{<:Distribution} = nothing,
    η::Union{Vector{T}, Nothing}  = nothing, ϕ::Union{Vector{T}, Nothing} = nothing,
    T̂::Union{Vector{T}, Nothing} = nothing, Ts::Union{Vector{T}, Nothing} = nothing,
    Ŝ::Union{Vector{T}, Nothing} = nothing) where {T<:Float64}

    n = []
    ζ ≠ nothing && push!(n, length(ζ))
    η ≠ nothing && push!(n, length(η))
    ϕ ≠ nothing && push!(n, length(ϕ))
    T̂ ≠ nothing && push!(n, length(T̂))
    Ŝ ≠ nothing && push!(n, size(Ŝ, 1))
    n = unique(n)

    length(n) === 1 || throw(ArgumentError("Supplied parameters of different lengths"))
    n = n[1]

    # construct distribution from matrix if available
    if (Ts === nothing) + (Ŝ === nothing) === 1
        error("Either both 'Ts' should be provided 'Ŝ' or neither")
    elseif Ts ≠ nothing && Ŝ ≠ nothing && ζ ≠ nothing
        @assert length(Ts) == size(Ŝ, 2)
        if Ts[1] ≠ 0
            # Set S(0) = 1
            Ŝ = hcat(ones(n), Ŝ)
            Ts = [0, Ts...]
        end
        # calculate pmf and create distribution
        p = mapreduce(x -> abs.(diff(x)), hcat, eachrow(Ŝ))'
        ζ = DiscreteNonParametric(Ts, p, check_args = false)
    end

    # no transformations assumed
    η = η === nothing ? fill(NaN, n) : η
    ϕ = ϕ === nothing ? fill(NaN, n) : ϕ
    T̂ = T̂ === nothing ? fill(NaN, n) : T̂

    if ζ === nothing
        return DeterministicSurvivalPrediction(η, ϕ, T̂)
    elseif ζ[1] isa ContinuousUnivariateDistribution
        return ContinuousSurvivalPrediction(ζ, η, ϕ, T̂)
    elseif ζ[1] isa DiscreteNonParametric
        if Ts === nothing
            Ts = unique(support.(ζ))
            @assert length(Ts) == 1 "Predicted distributions (ζ) must have same survival
times"
            Ŝ = mapreduce(s -> 1 .- cumsum(probs(s)), hcat, ζ)'
        end
        return DiscreteSurvivalPrediction(ζ, η, ϕ, T̂, (time = Ts, surv = Ŝ))
    end
end

Base.show(io::IO, sp::T where {T <: DeterministicSurvivalPrediction}) =
    print(io, DataFrame("lp" => sp.lp, "crank" => sp.crank, "time" => sp.time))
Base.show(io::IO, sp::T where {T <: DiscreteSurvivalPrediction}) =
    print(io, DataFrame("lp" => sp.lp, "crank" => sp.crank, "time" => sp.time,
                        "distr" => sp.distr,
                        "mat" => sp.survival_matrix))
Base.show(io::IO, sp::T where {T <: ContinuousSurvivalPrediction}) =
    print(io, DataFrame("lp" => sp.lp, "crank" => sp.crank, "time" => sp.time,
                        "distr" => sp.distr))
