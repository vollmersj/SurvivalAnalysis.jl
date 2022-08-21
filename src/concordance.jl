struct ConcordanceWeights
  S::Int8
  G::Int8
  tie_pred::Float16
  tie_time::Float16

  function ConcordanceWeights(S::Int, G::Int, tie_pred::Number, tie_time::Number)
    test_proportion(tie_pred) ||
      throw(ArgumentError("Expected 0≤`tie_pred`≤1, got $(tie_pred)"))
    test_proportion(tie_time) ||
      throw(ArgumentError("Expected 0≤`tie_time`≤1, got $(tie_time)"))
    new(convert(Int8, S), convert(Int8, G), convert(Float16, tie_pred),
      convert(Float16, tie_time))
  end
end

function concordance(truth::OneSidedSurv, prediction::Vector{<:Number}, weights::Symbol;
    tie_pred=0.5, tie_time=0, cutoff=nothing, train=nothing, rev=false,
    custom_weights::Union{Nothing, ConcordanceWeights}=nothing)

    # if everything is tied or no events then return 0.5
    ((length(unique(prediction)) == 1) || (length(unique_times(truth)) == 1) ||
      (total_events(truth) == 0)) && return 0.5

    if custom_weights !== nothing
      weights = custom_weights
    elseif weights === :I || weights === :Harrell
      weights = ConcordanceWeights(0, 0, tie_pred, tie_time)
    elseif weights === :G
      weights = ConcordanceWeights(0, -1, tie_pred, tie_time)
    elseif weights === :G2 || weights === :Uno
      weights = ConcordanceWeights(0, -2, tie_pred, tie_time)
    elseif weights === :GH || weights === :Gonen
      return gonen(sort(prediction), tie_pred)
    elseif weights === :SG || weights === :Schemper
      weights = ConcordanceWeights(1, -1, tie_pred, tie_time)
    elseif weights === :S || weights === :Peto
      weights = ConcordanceWeights(1, 0, tie_pred, tie_time)
    elseif weights == :G
      weights = ConcordanceWeights(0, -1, tie_pred, tie_time)
    else
      throw(AssertionError("`weights` must be one of `:I`, `:G`, `:G2`, :`SG:, `:S`, `:Uno`,
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

    out = _concordance(time, idx, prediction, cutoff, weights, cens, surv,
                tie_pred, tie_time)
    out = rev ? 1 - out : out
    return out
end

function _concordance(time, idx, crank, cutoff, weights, cens, surv, tie_pred, tie_time)
  num = 0
  den = 0

  for i in idx # only calculate when earlier time is an event
    time[i] > cutoff && break
    for j in (i+1):length(time)
      if time[i] == time[j] && tie_time == 0
        continue
      else
        weight_G = (weights.G == 0 || time[i] < cens.time[1]) ? 1 :
          cens.survival[searchsortedlast(cens.time, time[i])]^weights.G
        weight_S = (weights.S == 0 || time[i] < surv.time[1]) ? 1 :
          surv.survival[searchsortedlast(surv.time, time[i])]^weights.S
        weight = weight_G * weight_S

        den += weight

        if time[i] == time[j]
            num += weight * tie_time
        elseif time[i] < time[j]
          if crank[i] < crank[j]
            num += weight
          elseif crank[i] == crank[j]
            num += weight * tie_pred
          end
        end
      end
    end
  end

  return num/den
end

function gonen(crank, tie)
    # we assume crank to be sorted
    n = length(crank)
    ghci = 0.0

    for i in 1:(n-1)
      ci = crank[i]
      for j in (i+1):n
        cj = crank[j]
        ghci += ((ci < cj) ? 1 : tie) / (1 + exp(ci - cj))
      end
    end

    return 1 - (2 * ghci) / (n * (n - 1))
end

const cindex = concordance
