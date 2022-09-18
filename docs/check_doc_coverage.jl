using SurvivalAnalysis

documented = map(x -> replace(string(x), "SurvivalAnalysis." => ""),
    collect(keys(Docs.meta(SurvivalAnalysis))));

nms = names(SurvivalAnalysis);
expected = [];
map(x ->
    x != :SurvivalAnalysis && ## remove package name
    Symbol(eval(x)) == x && # remove aliases
    parentmodule(eval(x)) == SurvivalAnalysis && ## check if this module is parent module
    push!(expected, strip(String(x))), nms); ## push if all true

missings = setdiff(expected, documented); ## under-documented, critical error
extras = setdiff(documented, expected); ## over-documented, warning

if length(missings) > 0
    error("Forgot to document: ", missings)
elseif length(extras) > 0
    @warn "Extra objects documented:" extras
else
    @info "All objects documented as expected:" expected
end
