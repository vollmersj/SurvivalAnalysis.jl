"""
    ConcordanceWeights(S::Int8, G::Int8, tied_preds::Float64 tied_times::Float64)

Weights used in the [`concordance`](@ref) function. `S` and `G` reflect the power
applied to the Kaplan-Meier estimates of the survival and censoring distributions of the
fitted data used to apply weighting at a given time. For example
`ConcordanceWeights(1, -2, 0.5, 0.5)` will tell the `concordance` function to multiply
concordant pairs at time `t` by `S(t)/G(t)²`. `tied_preds` and `tied_times` determine how to
handle ties in predictions and observed times respectively.

See [`concordance`](@ref) for full examples.

Note❗It is strongly recommended that `S ≥ 0` and `G ≤ 0`.
"""
struct ConcordanceWeights
    S::Int8
    G::Int8
    tied_preds::Float64
    tied_times::Float64
    name::String

    ConcordanceWeights(S, G, tied_preds, tied_times) =
        ConcordanceWeights(S, G, tied_preds, tied_times, "")

    function ConcordanceWeights(S::Int, G::Int, tied_preds::Number, tied_times::Number,
        name::String)
        test_proportion(tied_preds) ||
            throw(ArgumentError("Expected 0≤`tied_preds`≤1, got $(tied_preds)"))
        test_proportion(tied_times) ||
            throw(ArgumentError("Expected 0≤`tied_times`≤1, got $(tied_times)"))
        new(convert(Int8, S), convert(Int8, G), convert(Float64, tied_preds),
        convert(Float64, tied_times), name)
    end
end

struct Concordance <: SurvivalMeasure
    C::Float64
    numerator::Float64
    denominator::Float64
    weights::ConcordanceWeights
    tied_times::Int64
    tied_preds::Int64
    tied_both::Int64
    pairs::Int64
    comparable::Int64
    concordant::Int64
    disconcordant::Int64
    reversed::Bool
    cutoff::Float64
end

function Base.show(io::IO, C::Concordance)
    println(io, typeof(C))
    println(io)
    println(io, lstrip("$(C.weights.name) C = $(C.C)"), " (reversed)"^C.reversed)
    println(io)
    println(io, "Cutoff: T ≤ $(C.cutoff)")
    println(io, "Counts:")
    pretty_table(io, [C.pairs C.comparable C.concordant C.disconcordant],
        header = ["Pairs", "Comparable", "Concordant", "Disconcordant"], vlines = :none,
        hlines = :none)
    println(io, "Weights:")
    str = c_pstring(pstring("S", C.weights.S), pstring("G", C.weights.G))
    pretty_table(io, [str C.weights.tied_preds C.weights.tied_times],
        header = ["IPCW", "Tied preds", "Tied times"], vlines = :none, hlines = :none)
    println(io, "Ties:")
    pretty_table(io, [C.tied_times C.tied_preds C.tied_both],
        header = ["Times", "Preds", "Both"], vlines = :none, hlines = :none)
    println(io, "Weighted calculation:")
    pretty_table(io, [C.numerator C.denominator C.C],
        header = ["Numerator", "Denominator", "C"], vlines = :none, hlines = :none)
    return nothing
end

"""
    concordance(
        truth::OneSidedSurv, prediction::Vector{<:Number}, weights::Union{Symbol, ConcordanceWeights};
        tied_preds=0.5, tied_times=0, cutoff=nothing, train::OneSidedSurv=nothing, rev=false
    )

    Aliases: cindex

Generic function to call any concordance index method. Concordance is a measure of
discrimination which evaluates if a prediction is concordant with the truth, i.e.
let Tᵢ,Tⱼ be the true survival times, `truth`, for observations `i` and `j` and let  ϕᵢ, ϕⱼ
be predictions, `prediction`, then these are concordant if ``ϕᵢ > ϕⱼ ⟺ Tᵢ > Tⱼ``.

In survival analysis the generic C-index is defined as follows

```math
C = \\frac{\\sum_{i≠j} W(tᵢ)I(tᵢ < tⱼ, yᵢ < yⱼ, tᵢ < τ)δᵢ}{W(tᵢ)I(tᵢ < tⱼ, tᵢ < τ)δᵢ}
```

where `tᵢ,tⱼ` are true survival times, `yᵢ,yⱼ` are predictions, `τ` is a cutoff time used to
ensure stability even when censoring is high, `δᵢ` is the censoring indicator and `W` is a
weighting function determined by `weights` as follows:

* `:I` or `:Harrell` - Harrell's C - ``W(tᵢ) = 1``
* `:G2` or `:Uno` - Uno's C - ``W(tᵢ) = 1/G(tᵢ)``
* `:SG` or `:Schemper` - Schemper's C - ``W(tᵢ) = S(tᵢ)/G(tᵢ)``
* `:S` or `:Peto` - Peto-Wilcoxon's C - ``W(tᵢ) = S(tᵢ)``

where `S(tᵢ)` and `G(tᵢ)` are respectively the Kaplan-Meier estimates of the survival and
censoring distributions of the training data at time `tᵢ`. For any other combination of
weights pass a [`ConcordanceWeights`](@ref) object to `weights`.
Note❗If training data is not provided to `train` then `truth` is used to estimate `S` and
`G` but it is strongly recommended to provide the training data if possibe.

We also include an implementation of Gönen-Heller's C with `weights = :GH` or `:Gonen`.
Note❗This is actually a very different method from the others and calculates concordance
for predictions from a Cox PH model only. We may move this to its own function in the future
to avoid misuse.

There is open debate about how to handle ties when calculating the concordance. Ties can
occur in predictions and in the observed survival times. The defaults here `tied_preds=0.5`
and `tied_times=0` are set as these seem to be the most common but they can be changed.
Note❗If you pass a [`ConcordanceWeights`](@ref) object then the tied
weights specified in this will take priority.

Note❗For predictions from PH models or any model where the prediction represents a relative
risk then a higher value of ϕ implies a higher risk of event which will result in a *lower*
survival time. In this case a prediction is concordant with the survival time if
``ϕᵢ < ϕⱼ ⟺ Tᵢ > Tⱼ``. To do this within the function just set `rev=true`.

# Examples
```jldoctest
julia> T = [1.0,0.1,pi,0.9,0.4,20,1,5,9,2.5];

julia> ϕ = [0.1,0.2,0.1,0.9,0.25,exp(2),8,log(9),2,8];

julia> Δ = [true,true,false,true,true,false,true,false,true,true];

julia> Y = Surv(T, Δ, :r);

julia> train = Surv([3,3,9,12,2,1,8,4,10,8], Δ, :r);

julia> concordance(Y, ϕ, cutoff = threshold_risk(Y, 0.8)) # Harrell's C cutoff when 80% data is censored or dead
SurvivalAnalysis.Concordance

Harrell's C = 0.6153846153846154

Cutoff: T ≤ 9.0
Counts:
 Pairs  Comparable  Concordant  Disconcordant
    90          39          23             14
Weights:
 IPCW  Tied preds  Tied times
    1         0.5         0.0
Ties:
 Times  Preds  Both
     1      2     0
Weighted calculation:
 Numerator  Denominator         C
      24.0         39.0  0.615385

julia> concordance(Y, ϕ, :Uno, train=train) # Uno's C
SurvivalAnalysis.Concordance

Uno's C = 0.6221374045801525

Cutoff: T ≤ 20.0
Counts:
 Pairs  Comparable  Concordant  Disconcordant
    90          39          23             14
Weights:
 IPCW  Tied preds  Tied times
 1/G²         0.5         0.0
Ties:
 Times  Preds  Both
     1      2     0
Weighted calculation:
 Numerator  Denominator         C
   28.1728       45.284  0.622137

julia> cindex(Y, ϕ, ConcordanceWeights(5, -3, 0.5, 0.5, "Silly Weights"); train=train) # Custom weights
SurvivalAnalysis.Concordance

Silly Weights C = 0.6070334788405888

Cutoff: T ≤ 20.0
Counts:
 Pairs  Comparable  Concordant  Disconcordant
    90          40          23             14
Weights:
  IPCW  Tied preds  Tied times
 S⁵/G³         0.5         0.5
Ties:
 Times  Preds  Both
     1      2     0
Weighted calculation:
 Numerator  Denominator         C
   25.6265       42.216  0.607033
```
"""
function concordance(truth::OneSidedSurv, prediction::Vector{<:Number},
    weights::Union{Symbol, ConcordanceWeights}=:I;
    tied_preds=0.5, tied_times=0, cutoff=nothing,
    train::Union{Nothing, OneSidedSurv}=nothing, rev=false)

    n_truth = length(truth)
    n_pred = length(prediction)
    n_pairs = binomial(n_pred, 2)*2
    n_truth != n_pred &&
        throw(ArgumentError(
            "`truth` (n=$(n_truth)) and `prediction` (n=$(n_pred)) unequal lengths"))
    cutoff = cutoff === nothing ? maximum(truth.time) : cutoff

    # if everything is tied or no events then return 0.5
    msg = nothing
    if length(unique(prediction)) == 1
        msg = "all predictions equal"
    elseif length(unique_event_times(truth)) == 1
        msg = "all observed event times equal"
    elseif total_events(truth) == 0
        msg = "no observed events"
    end

    if msg !== nothing
        weights = ConcordanceWeights(0, 0, 0, 0, string("Fallback ($(msg))"))
        out = (numerator=0, denominator=0, tied_times=0, tied_preds=0, tied_both=0,
            comparable=0, concordant=0, disconcordant=0)
        C = 0.5
    else
        if typeof(weights) === Symbol
            if weights === :I || weights === :Harrell
                weights = ConcordanceWeights(0, 0, tied_preds, tied_times, "Harrell's")
            elseif weights === :G2 || weights === :Uno
                weights = ConcordanceWeights(0, -2, tied_preds, tied_times, "Uno's")
            elseif weights === :GH || weights === :Gonen
                weights = ConcordanceWeights(0, 0, 0, 0, "Gönen-Heller's")
                out = gonen(sort(prediction), tied_preds)
            elseif weights === :SG || weights === :Schemper
                weights = ConcordanceWeights(1, -1, tied_preds, tied_times, "Schemper's")
            elseif weights === :S || weights === :Peto
                weights = ConcordanceWeights(1, 0, tied_preds, tied_times, "Peto-Wilcoxon's")
            else
                throw(ArgumentError("`weights` must be one of `:I`, `:G2`, :`SG:, `:S`,
`:Uno`, `:Harrell`, `:GH`, `:Gonen`, `:Schemper`, `:Peto`, or be a `ConcordanceWeights`
object"))
            end
        end

        if weights.name != "Gönen-Heller's"
            ord = sortperm(truth.time)
            time = truth.time[ord]

            prediction = prediction[ord]
            idx = findall(!iszero, truth.status[ord])

            train = train === nothing ? truth : train
            surv = weights.S != 0 ? fit(KaplanMeier, zeros(0,0), train) : nothing
            cens = weights.G != 0 ? fit(KaplanMeier, zeros(0,0), reverse(train)) : nothing



            out = _concordance(time, idx, prediction, cutoff, weights, cens, surv)
        end

        C = rev ? 1 - out.numerator/out.denominator : out.numerator/out.denominator
    end

    return Concordance(
        C,
        out.numerator,
        out.denominator,
        weights,
        out.tied_times,
        out.tied_preds,
        out.tied_both,
        n_pairs,
        out.comparable,
        out.concordant,
        out.disconcordant,
        rev,
        cutoff
    )
end

function _concordance(time, idx, pred, cutoff, wts, cens, surv)
    num = 0 # weighted numerator
    den = 0 # weighted denominator
    Nₜ = 0 # numbered tied time
    Nₚ = 0 # number tied pred
    Nₜₚ = 0 # number tied time and pred
    N₊ = 0 # concordant
    N₋ = 0 # disconcordant
    N₌ = 0 # comparable

    for i in idx # only calculate when earlier time is an event
        time[i] > cutoff && break
        for j in (i+1):length(time)
        if (time[i] == time[j] && wts.tied_times == 0) || (pred[i] == pred[j] && wts.tied_preds == 0)
            if time[i] == time[j] && pred[i] == pred[j]
            Nₜₚ += 1
            elseif time[i] == time[j]
            Nₜ += 1
            elseif pred[i] == pred[j]
            Nₚ += 1
            end
            continue
        else
            weight_G = (wts.G == 0 || time[i] < cens.time[1]) ? 1 :
            cens.survival[searchsortedlast(cens.time, time[i])]^wts.G
            weight_S = (wts.S == 0 || time[i] < surv.time[1]) ? 1 :
            surv.survival[searchsortedlast(surv.time, time[i])]^wts.S
            weight = weight_G * weight_S

            den += weight
            N₌ += 1

            if time[i] == time[j]
                num += weight * wts.tied_times
                if pred[i] == pred[j]
                Nₜₚ += 1
                else
                Nₜ += 1
                end
            elseif time[i] < time[j]
            if pred[i] < pred[j]
                num += weight
                N₊ += 1
            elseif pred[i] == pred[j]
                num += weight * wts.tied_preds
                Nₚ += 1
            else
                N₋ += 1
            end
            end
        end
        end
    end

    return (numerator=num, denominator=den, tied_times=Nₜ, tied_preds=Nₚ, tied_both=Nₜₚ,
        comparable=N₌, concordant=N₊, disconcordant=N₋)
end

function gonen(pred, tie)
    # we assume pred to be sorted
    n = length(pred)
    ghci = 0.0
    Nₚ = 0 # tied preds
    N₊ = 0 # concordant
    N₌ = binomial(n, 2) # comparable

    for i in 1:(n-1)
        ci = pred[i]
        for j in (i+1):n
            cj = pred[j]
            if ci == cj
            Nₚ += 1
            ghci += tie / (1 + exp(ci - cj))
            else
            N₊ += 1
            ghci += 1 / (1 + exp(ci - cj))
            end
        end
    end

    num = (n * (n - 1)) - (2 * ghci) / (n * (n - 1))
    den = (n * (n - 1))

    N₋ = N₌ - N₊

    return (numerator=num, denominator=den, tied_times=0, tied_preds=Nₚ, tied_both=0,
        comparable=N₌, concordant=N₊, disconcordant=N₋)
end

const cindex = concordance
