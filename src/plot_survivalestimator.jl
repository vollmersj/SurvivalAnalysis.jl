@recipe function f(
    npe::SurvivalEstimator, plot_confint::Bool = true;
    α::Float64 = 0.05
)
    seriestype := :steppost
    ylims := (0, 1)
    legend := false
    @series begin
        linecolor   --> :black
        npe.time, npe.survival
    end
    if plot_confint
        linecolor   --> :blue
        cis = confint.(Ref(npe), npe.time, α = α)
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
