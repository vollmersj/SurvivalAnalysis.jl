abstract type SurvivalModel <: StatisticalModel end

# Stripped down version of https://github.com/JuliaStats/StatsModels.jl/blob/5299399a903cdb8ea4b9492a9c09db72f703c7b7/src/statsmodel.jl#L172
#  for SurvivalModels (licensed under MIT)
StatsBase.predict(fit::SurvivalModel, data::DataFrame; kwargs...) =
    predict(fit, Matrix(data); kwargs...)

function StatsModels.predict(
    mm::StatsModels.TableStatisticalModel{<:SurvivalModel, <:AbstractMatrix}, data;
    kwargs...
)
    Tables.istable(data) ||
        throw(ArgumentError("expected data in a Table, got $(typeof(data))"))

    f = mm.mf.f
    cols, nonmissings = StatsModels.missing_omit(StatsModels.columntable(data), f.rhs)
    new_x = modelcols(f.rhs, cols)
    new_x = reshape(new_x, size(new_x, 1), :)[:,2:end] # remove intercept
    return StatsModels.predict(mm.model, new_x; kwargs...)
end

StatsBase.coef(mm::StatsModels.TableStatisticalModel{<:SurvivalModel, <:AbstractMatrix}) =
    coef(mm.model)

baseline(mm::StatsModels.TableStatisticalModel{<:SurvivalModel, <:AbstractMatrix}) =
    mm.model.baseline
baseline(mm::SurvivalModel) = mm.baseline

Distributions.scale(mm::StatsModels.TableStatisticalModel{<:SurvivalModel, <:AbstractMatrix}) =
    mm.model.scale
Distributions.scale(mm::SurvivalModel) = mm.scale

function StatsModels.coeftable(mm::SurvivalModel)
    header = ["(Scale)", ["x$i" for i = 0:(length(coef(mm))-1)]...]
    header[2] = "(Intercept)"
    data = [mm.scale coef(mm)...]
    pretty_table(data,  header = header, vlines = :none, hlines = :none)
end

function StatsModels.coeftable(mm::StatsModels.TableStatisticalModel{<:SurvivalModel})
    header = ["(Scale)", coefnames(mm.mf)...]
    data = [scale(mm) coef(mm)...]
    pretty_table(data,  header = header, vlines = :none, hlines = :none)
end

function Base.show(io::IO, mm::SurvivalModel)
    println(io, typeof(mm))
    println(io)
    println(io, "Distr:")
    println(io, baseline(mm))
    println(io)
    println(io,"Coefficients:")
    coeftable(mm)
end

function Base.show(io::IO, mm::StatsModels.TableStatisticalModel{<:SurvivalModel})
    println(io, typeof(mm))
    println(io)
    println(io, mm.mf.f)
    println(io)
    println(io, "Distr:")
    println(io, baseline(mm))
    println(io)
    println(io, "Coefficients:")
    coeftable(mm)
end
