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
    survival_matrix::NamedTuple{(:time, :survival), Tuple{Vector{T}, Matrix{T}}}
end

struct ContinuousSurvivalPrediction{T<:Float64} <: SurvivalPrediction
    distr::Vector{<:ContinuousUnivariateDistribution}
    lp::Vector{T}
    crank::Vector{T}
    time::Vector{T}
end

function SurvivalPrediction(;
    distr::Union{Nothing, Vector{<:Distribution}} = nothing,
    lp::Union{Nothing, Vector{T}}  = nothing,
    crank::Union{Nothing, Vector{T}} = nothing,
    time::Union{Nothing, Vector{T}} = nothing,
    fit_times::Union{Nothing, Vector{T}} = nothing,
    survival_matrix::Union{Matrix{T}, Nothing} = nothing
) where {T<:Number}

    n = []
    distr ≠ nothing && push!(n, length(distr))
    lp ≠ nothing && push!(n, length(lp))
    crank ≠ nothing && push!(n, length(crank))
    time ≠ nothing && push!(n, length(time))
    survival_matrix ≠ nothing && push!(n, size(survival_matrix, 1))
    n = unique(n)

    length(n) === 1 || throw(ArgumentError("Supplied parameters of different lengths"))
    n = n[1]

    # construct distr from matrix if available
    if (fit_times === nothing) + (survival_matrix === nothing) === 1
        throw(ArgumentError(("Either both 'fit_times' should be provided 'survival_matrix' or neither")))
    elseif fit_times ≠ nothing && survival_matrix ≠ nothing && distr === nothing
        if length(fit_times) != size(survival_matrix, 2)
            throw(ArgumentError("Length of 'fit_times' must equal number of columns of 'survival_matrix'"))
        end
        ζₛ = survival_matrix
        ζₜ = fit_times
        if ζₜ[1] ≠ 0
            # Set S(0) = 1
            ζₛ = hcat(ones(n), ζₛ)
            ζₜ = [0, ζₜ...]
        end
        # calculate pmf and create distr
        distr = map(x -> DiscreteNonParametric(ζₜ, [1 - x[1], abs.(diff(x))...];
                check_args=false), eachrow(ζₛ))
    end

    # no transformations assumed
    lp = lp === nothing ? fill(NaN, n) : lp
    crank = crank === nothing ? fill(NaN, n) : crank
    time = time === nothing ? fill(NaN, n) : time

    if distr === nothing
        return DeterministicSurvivalPrediction(lp, crank, time)
    elseif distr[1] isa ContinuousUnivariateDistribution
        return ContinuousSurvivalPrediction(distr, lp, crank, time)
    elseif distr[1] isa DiscreteNonParametric
        if fit_times === nothing
            fit_times = unique(support.(distr))
            length(fit_times) == 1 || throw(ArgumentError(
                "Predicted distributions (distr) must have same survival times"))
            fit_times = fit_times[1]
            survival_matrix = Matrix(mapreduce(s -> 1 .- cumsum(probs(s)), hcat, distr)')
        end
        return DiscreteSurvivalPrediction(distr, lp, crank, time,
                                        (time = fit_times, survival = survival_matrix))
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
