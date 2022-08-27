"""
    plot(npe::SurvivalEstimator, plot_confint::Bool = true; level = 0.95)
    plot(npe::StatsModels.TableStatisticalModel{SurvivalEstimator, Matrix{Float64}},
        plot_confint::Bool = true; level = 0.95)

Recipe for plotting fitted non-parametric estimators, `npe`. If `plot_confint` then
confidence intervals also plotted at a `level`% confidence level.

# Examples
```jldoctest
julia> using Plots

julia> data = DataFrame(t = randn(10), d = [trues(5)..., falses(5)...]);

julia> plot(kaplan_meier(@formula(Srv(t, d) ~ 1), data));

julia> plot(nelson_aalen(@formula(Srv(t, d) ~ 1), data).model);
```
"""
@recipe function f(npe::SurvivalEstimator, plot_confint::Bool = true; level = 0.95)
    test_proportion(level) || throw(ArgumentError("level must be a number in [0, 1]"))
    seriestype := :steppost
    ylims := (0, 1)
    legend := false
    @series begin
        linecolor   --> :black
        npe.time, npe.survival
    end
    if plot_confint
        linecolor   --> :blue
        cis = confint.(Ref(npe), npe.time; level = level)
        lb = map(x -> x[1], cis)
        ub = map(x -> x[2], cis)
        @series begin
            npe.time, lb
        end
        @series begin
            npe.time, ub
        end
    end
    return nothing
end

@recipe function f(
        npe::StatsModels.TableStatisticalModel{<:SurvivalEstimator, Matrix{Float64}},
        plot_confint::Bool = true; level = 0.95
        )
        test_proportion(level) || throw(ArgumentError("level must be a number in [0, 1]"))
        npe = npe.model
        seriestype := :steppost
        ylims := (0, 1)
        legend := false
        @series begin
            linecolor   --> :black
            npe.time, npe.survival
        end
        if plot_confint
            linecolor   --> :blue
            cis = confint.(Ref(npe), npe.time; level = level)
            lb = map(x -> x[1], cis)
            ub = map(x -> x[2], cis)
            @series begin
                npe.time, lb
            end
            @series begin
                npe.time, ub
            end
        end
        return nothing
end
