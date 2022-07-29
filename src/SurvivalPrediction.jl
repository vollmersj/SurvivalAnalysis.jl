abstract type SurvivalPrediction end

struct DeterministicSurvivalPrediction{T<:Float64} <: SurvivalPrediction
    lp::Vector{T}
    crank::Vector{T}
    time::Vector{T}
end

struct DiscreteSurvivalPrediction{T<:Float64} <: SurvivalPrediction
    distr::Vector{<:DiscreteNonParametric}
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

function SurvivalPrediction(;
    ζ::Union{Nothing, Vector{<:Distribution}} = nothing,
    η::Union{Nothing, Vector{T}}  = nothing,
    ϕ::Union{Nothing, Vector{T}} = nothing,
    T̂::Union{Nothing, Vector{T}} = nothing,
    Ts::Union{Nothing, Vector{T}} = nothing,
    Ŝ::Union{Matrix{T}, Nothing} = nothing
) where {T<:Number}

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
        throw(ArgumentError(("Either both 'Ts' should be provided 'Ŝ' or neither")))
    elseif Ts ≠ nothing && Ŝ ≠ nothing && ζ === nothing
        if length(Ts) != size(Ŝ, 2)
            throw(ArgumentError("Length of 'Ts' must equal number of columns of 'Ŝ'"))
        end
        ζₛ = Ŝ
        ζₜ = Ts
        if ζₜ[1] ≠ 0
            # Set S(0) = 1
            ζₛ = hcat(ones(n), ζₛ)
            ζₜ = [0, ζₜ...]
        end
        # calculate pmf and create distribution
        ζ = map(x -> DiscreteNonParametric(ζₜ, [1 - x[1], abs.(diff(x))...];
                check_args=false), eachrow(ζₛ))
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
            length(Ts) == 1 || throw(ArgumentError(
                "Predicted distributions (ζ) must have same survival times"))
            Ts = Ts[1]
            Ŝ = Matrix(mapreduce(s -> 1 .- cumsum(probs(s)), hcat, ζ)')
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
