abstract type SurvivalDistribution <: UnivariateDistribution end
abstract type SurvivalDistribution <: ContinuousUnivariateDistribution end
abstract type DiscreteSurvivalDistribution <: ContinuousUnivariateDistribution end

Base.show(io::IO, oss::OneSidedSurv) =
    print(io, map((t, δ) -> string(t, δ ? "" : "+"), oss.time, oss.status))
