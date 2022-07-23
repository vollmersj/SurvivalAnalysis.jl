using SurvivalAnalysis
using Distributions
using DataFrames

n = 20
Δ = fill(1, n);
T_Exp = rand(Exponential(5), n);
Y = Surv(T_Exp, Δ, "right");
data = DataFrame(
    "X" => randn(n),
    "T" => T_Exp,
    "δ" => Δ
)

# Method 1 - ph(::formula)
ph(@formula(Srv(T, δ) ~ X), data, Exponential)
# ph(@formula(Srv(T, δ) ~ X), data, Exponential, init = 0.1)

# # Method 1 - ph(::matrix)
# f2=ph(hcat(ones(20), reshape(data.X, 20, 1)), Y, Exponential)

# predict(fit, data)

# # # Method 1 - ph(::formula)
# # fit = ph(@formula(Srv(T, δ) ~ X), data, Exponential)
