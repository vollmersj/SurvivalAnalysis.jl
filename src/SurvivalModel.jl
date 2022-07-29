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
