struct ConcordanceWeights
  S::Int8
  G::Int8
  tie_pred::Float16
  tie_time::Float16
  name::String

  ConcordanceWeights(S, G, tie_pred, tie_time) =
    ConcordanceWeights(S, G, tie_pred, tie_time, "")

  function ConcordanceWeights(S::Int, G::Int, tie_pred::Number, tie_time::Number,
    name::String)
    test_proportion(tie_pred) ||
      throw(ArgumentError("Expected 0≤`tie_pred`≤1, got $(tie_pred)"))
    test_proportion(tie_time) ||
      throw(ArgumentError("Expected 0≤`tie_time`≤1, got $(tie_time)"))
    new(convert(Int8, S), convert(Int8, G), convert(Float16, tie_pred),
      convert(Float16, tie_time), name)
  end
end

struct Concordance
  C::Float64
  weights::ConcordanceWeights
  tied_times::Int64
  tied_pred::Int64
  tied_both::Int64
  comparable::Int64
  concordant::Int64
  disconcordant::Int64
  reversed::Bool
end

function Base.show(io::IO, C::Concordance)
  println(io, typeof(C))
  println(io)
  println(io, lstrip("$(C.weights.name) C = $(C.C)"), " (reversed)"^C.reversed)
  println(io)
  println(io, "Ties:")
  pretty_table(io, [C.tied_times C.tied_pred C.tied_both],
      header = ["Times", "Preds", "Both"], vlines = :none, hlines = :none)
  println(io, "Counts:")
  pretty_table(io, [C.comparable C.concordant C.disconcordant],
    header = ["Comparable", "Concordant", "Disconcordant"], vlines = :none, hlines = :none)
  println(io, "Weights:")
  pretty_table(io, [c_pstring(pstring("S", C.weights.S), pstring("G", C.weights.G)) C.weights.tie_pred C.weights.tie_time],
    header = ["IPCW", "Tied preds", "Tied times"], vlines = :none, hlines = :none)
  nothing
end

function concordance(truth::OneSidedSurv, prediction::Vector{<:Number}, weights::Symbol;
    tie_pred=0.5, tie_time=0, cutoff=nothing, train=nothing, rev=false,
    custom_weights::Union{Nothing, ConcordanceWeights}=nothing)

    n_truth = length(truth)
    n_pred = length(prediction)
    n_truth != n_pred &&
      throw(ArgumentError(
        "`truth` (n=$(n_truth)) and `prediction` (n=$(n_pred)) unequal lengths"))

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
      msg = string("Fallback ($(msg))")
      return Concordance(0.5, ConcordanceWeights(0, 0, 0, 0, msg), 0, 0, 0, 0, 0, 0, false)
    end

    if custom_weights !== nothing
      weights = custom_weights
    elseif weights === :I || weights === :Harrell
      weights = ConcordanceWeights(0, 0, tie_pred, tie_time, "Harrell's")
    elseif weights === :G
      weights = ConcordanceWeights(0, -1, tie_pred, tie_time)
    elseif weights === :G2 || weights === :Uno
      weights = ConcordanceWeights(0, -2, tie_pred, tie_time, "Uno's")
    elseif weights === :GH || weights === :Gonen
      return gonen(sort(prediction), tie_pred, rev)
    elseif weights === :SG || weights === :Schemper
      weights = ConcordanceWeights(1, -1, tie_pred, tie_time, "Schemper's")
    elseif weights === :S || weights === :Peto
      weights = ConcordanceWeights(1, 0, tie_pred, tie_time, "Peto-Wilcoxon's")
    else
      throw(ArgumentError("`weights` must be one of `:I`, `:G`, `:G2`, :`SG:, `:S`, `:Uno`,
`:Harrell`, `:GH`, `:Gonen`, `:Schemper`, `:Peto`"))
    end

    ord = sortperm(truth.time)
    time = truth.time[ord]
    prediction = prediction[ord]
    idx = findall(!iszero, truth.status[ord])

    train = train === nothing ? truth : train
    surv = weights.S != 0 ? fit(KaplanMeier, zeros(0,0), train) : nothing
    cens = weights.G != 0 ? fit(KaplanMeier, zeros(0,0), reverse(train)) : nothing

    cutoff = cutoff === nothing ? maximum(time) + 1 : cutoff

    return _concordance(time, idx, prediction, cutoff, weights, cens, surv,
                tie_pred, tie_time, rev)
end

function _concordance(time, idx, pred, cutoff, weights, cens, surv, tie_pred, tie_time, rev)
  num = 0 # weighted numerator
  den = 0 # weighted denominator
  Nₜ = 0 # numbered tied time
  Nₚ = 0 # number tied pred
  Nₜₚ = 0 # number tied time and pred
  N₊ = 0 # concordant
  N₌ = 0 # comparable

  for i in idx # only calculate when earlier time is an event
    time[i] > cutoff && break
    for j in (i+1):length(time)
      if (time[i] == time[j] && tie_time == 0) || (pred[i] == pred[j] && tie_pred == 0)
        if time[i] == time[j] && pred[i] == pred[j]
          Nₜₚ += 1
        elseif time[i] == time[j]
          Nₜ += 1
        elseif pred[i] == pred[j]
          Nₚ += 1
        end
        continue
      else
        weight_G = (weights.G == 0 || time[i] < cens.time[1]) ? 1 :
          cens.survival[searchsortedlast(cens.time, time[i])]^weights.G
        weight_S = (weights.S == 0 || time[i] < surv.time[1]) ? 1 :
          surv.survival[searchsortedlast(surv.time, time[i])]^weights.S
        weight = weight_G * weight_S

        den += weight
        N₌ += 1

        if time[i] == time[j]
            num += weight * tie_time
            Nₜ += 1
        elseif time[i] < time[j]
          if pred[i] < pred[j]
            num += weight
            N₊ += 1
          elseif pred[i] == pred[j]
            num += weight * tie_pred
            Nₚ += 1
          end
        end
      end
    end
  end

  C = rev ? 1 - num/den : num/den
  N₋ = N₌ - N₊ - Nₚ - Nₜ - Nₜₚ # disconcordant
  return Concordance(C, weights, Nₜ, Nₚ, Nₜₚ, N₌, N₊, N₋, rev)
end

function gonen(pred, tie, rev)
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

    C = rev ? (2 * ghci) / (n * (n - 1)) : 1 - (2 * ghci) / (n * (n - 1))
    N₋ = N₌ - N₊ - Nₚ # disconcordant
    return Concordance(C, ConcordanceWeights(0, 0, 0, 0, "Gönen-Heller's"),
                      0, Nₚ, 0, N₌, N₊, N₋, rev)
end

const cindex = concordance
