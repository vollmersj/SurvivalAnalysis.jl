# struct for behavior
struct SurvTerm{X, Y} <: AbstractTerm
    T::X
    Δ::Y
    type::ConstantTerm
end
SurvTerm(T::X where {X<:AbstractTerm}, Δ::Y where {Y<:AbstractTerm}) =
    SurvTerm(T, Δ, term(1))

Base.show(io::IO, t::SurvTerm) =
    print(io, string("(", t.T, ",",  t.Δ, t.type.n == 1 ? ";+" : ";-", ")"))

Srv(T::Symbol, Δ::Symbol, type::Int = 1) = SurvTerm(term(T), term(Δ), term(type))
Srv(T::Symbol) = SurvTerm(term(T), term(trues(length(T))), term(1))

function StatsModels.apply_schema(t::FunctionTerm{typeof(Srv)},
                                    sch::StatsModels.Schema,
                                    Mod::Type{<:Any})
    return apply_schema(SurvTerm(t.args_parsed...), sch, Mod)
end

function StatsModels.apply_schema(t::SurvTerm,
                                    sch::StatsModels.Schema,
                                    Mod::Type{<:Any})
    T = apply_schema(t.T, sch, Mod)
    isa(T, ContinuousTerm)  ||
        throw(ArgumentError("SurvTerm only works with continuous terms (got $T)"))
    Δ = apply_schema(t.Δ, sch, Mod)
    isa(Δ, ContinuousTerm) ||
        throw(ArgumentError("SurvTerm only works with discrete terms (got $Δ)"))
    @assert abs(t.type.n) == 1 "Srv type must be 1 (right) or -1 (left) censoring"
    return SurvTerm(T, Δ, t.type)
end

function StatsModels.modelcols(t::SurvTerm, d::NamedTuple)
    T = modelcols(t.T, d)
    Δ = modelcols(t.Δ, d)
    return Surv(T, Δ, t.type.n == 1 ? "right" : "left")
end

StatsModels.terms(t::SurvTerm) = [terms(t.T), terms(t.Δ), terms(t.type)]
StatsModels.termvars(t::SurvTerm) = [StatsModels.termvars(t.T), StatsModels.termvars(t.Δ),
                                        StatsModels.termvars(t.type)]
StatsModels.width(t::SurvTerm) = 1
StatsBase.coefnames(t::SurvTerm) = [coefnames(t.T), coefnames(t.Δ), "(type)"]
