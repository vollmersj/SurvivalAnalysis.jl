"""
    SurvivalTimeMeasure

Abstract type for all survival time measures implemented in, or extending, this package. A survival time measure is any measure of the form ``L(t, t̂)`` where `t` is the observed event time for an observation and `t̂` is a prediction of the survival time.

Note❗Survival time measures are poorly researched in the literature and no good studies demonstrate their utility. Nor is there concensus on the 'right' way to implement these. Here the measures remove all censored observations and then compare predicted survival times for uncensored observations. Alternatives could include: using an IPCW weighting method, or retaining censoring observations and applying an uncertainty penalty.

# Fields
Structs inheriting from this type include the following fields

* `losses` - Pointwise losses for the measure, i.e. ``L(tᵢ, t̂ᵢ)``
* `mean` - Aggregation of losses, usually the sample mean
* `se` - Standard error of the mean

Note❗ [`RMSE`](@ref) is an aggregated measure which means `mean` is not simply the sample
mean over the losses and instead `mean` is just used as a shorthand.
"""
abstract type SurvivalTimeMeasure end

Base.show(io::IO, t::SurvivalTimeMeasure) =
    println(io, "$(typeof(t)) = $(round(t.mean; digits=4)) (σ = $(round(t.se; digits=4)))")

struct MSE <: SurvivalTimeMeasure
    losses::Vector{Float64}
    mean::Float64
    se::Float64
end

struct MAE <: SurvivalTimeMeasure
    losses::Vector{Float64}
    mean::Float64
    se::Float64
end

struct RMSE <: SurvivalTimeMeasure
    losses::Vector{Float64}
    mean::Float64
    se::Float64
end

"""
    MSE(truth::OneSidedSurv, prediction::Vector{<:Number})

Calculates the Mean Squared Error (MSE), defined by

```math
MSEᵢ = δᵢ(tᵢ - t̂ᵢ)^2 \\
mean = \\frac{1}{m} ∑ᵢ MSEᵢ \\
se = \\frac{σ}{√m}
```
where `t` is the true survival time, `t̂` is the predicted survival time, `δ` is the censoring indicator, ``m = ∑ᵢ δᵢ``, and `σ` is the sample standard deviation of `MSE`.

Returns a `MSE` struct inheriting from [`SurvivalTimeMeasure`](@ref).

Note❗Censored observations in the test set are ignored.

# Examples
```jldoctest
julia> truth = Surv([1, 4, 2, 2, 8], [true, false, false, true, false], :r)

julia> pred = [pi, sqrt(3), 1, 8, 9]

julia> MSE(truth, pred)
MSE = 20.2932 (σ = 15.7068)
```
"""
function MSE(truth::OneSidedSurv, prediction::Vector{<:Number})
    Δ = outcome_status(truth)
    mse = (outcome_times(truth)[Δ] - prediction[Δ]).^2
    MSE(mse, mean(mse), std(mse) / sqrt(length(mse)))
end

"""
    MAE(truth::OneSidedSurv, prediction::Vector{<:Number})

Calculates the Mean Absolute Error (MAE), defined by

```math
MAEᵢ = δᵢ|tᵢ - t̂ᵢ| \\
mean = \\frac{1}{m} ∑ᵢ MAEᵢ \\
se = \\frac{σ}{√m}
```
where `t` is the true survival time, `t̂` is the predicted survival time, `δ` is the censoring indicator, ``m = ∑ᵢ δᵢ``, and `σ` is the sample standard deviation of `MAE`.

Returns a `MAE` struct inheriting from [`SurvivalTimeMeasure`](@ref).

Note❗Censored observations in the test set are ignored.

# Examples
```jldoctest
julia> truth = Surv([1, 4, 2, 2, 8], [true, false, false, true, false], :r)

julia> pred = [pi, sqrt(3), 1, 8, 9]

julia> MAE(truth, pred)
MAE = 4.0708 (σ = 1.9292)
```
"""
function MAE(truth::OneSidedSurv, prediction::Vector{<:Number})
    Δ = outcome_status(truth)
    mae = abs.(outcome_times(truth)[Δ] - prediction[Δ])
    MAE(mae, mean(mae), std(mae) / sqrt(length(mae)))
end

"""
    RMSE(truth::OneSidedSurv, prediction::Vector{<:Number})

Calculates the Root Mean Squared Error (RMSE), defined by

```math
MSEᵢ = δᵢ(tᵢ - t̂ᵢ)^2 \\
RMSE = \\sqrt{\\frac{1}{m} ∑ᵢ MSEᵢ} \\
se = σ/(2√m * rmse)
```
where `t` is the true survival time, `t̂` is the predicted survival time, `δ` is the censoring indicator, ``m = ∑ᵢ δᵢ``, and `σ` is the sample standard deviation of `MSE`.

Returns a `RMSE` struct inheriting from [`SurvivalTimeMeasure`](@ref).

Note❗Censored observations in the test set are ignored.
Note❗RMSE is an aggregated measure which means `losses` actually correspond to `MSEᵢ`.

# Examples
```jldoctest
julia> truth = Surv([1, 4, 2, 2, 8], [true, false, false, true, false], :r)

julia> pred = [pi, sqrt(3), 1, 8, 9]

julia> RMSE(truth, pred)
RMSE = 4.5048 (σ = 1.7433)
```
"""
function RMSE(truth::OneSidedSurv, prediction::Vector{<:Number})
    mse = MSE(truth, prediction)
    rmse = sqrt(mse.mean)
    RMSE(mse.losses, rmse, mse.se / (2 * rmse))
end
