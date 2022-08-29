using Gumbo
using Cascadia
using SurvivalAnalysis

h = parsehtml(read("docs/build/api/index.html", String)).root[2];
txt = eachmatch(sel".docstring-binding code", h);
documented = map(x -> replace(string(x[1]), "SurvivalAnalysis." => ""), txt);

f = readchomp("src/SurvivalAnalysis.jl");
exp = f[(findfirst("# documented", f)[end]+1):(findlast("# documented", f)[1]-1)];
exp = filter(x -> x != "", strip.(split(exp, r"export|\n|,")));

ed_diff = setdiff(exp, documented);
de_diff = setdiff(documented, exp)

if length(ed_diff) > 0
    @error "Forgot to document:" ed_diff
elseif length(de_diff) > 0
    @warn "Extra objects documented:" de_diff
else
    @info "All objects documented as expected:" exp
end
