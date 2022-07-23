seed!(42)
n = 500
Δ = fill(1, n);
T_Exp = rand(Exponential(5), n);
T_Weib = rand(Weibull(2, 5), n);
Xs = Uniform();
data = DataFrame(
    "X1" => rand(Xs, n),
    "X2" => rand(Xs, n),
    "X3" => rand(Xs, n),
    "T_Exp" => T_Exp,
    "T_Weib" => T_Weib,
    "δ" => Δ);

weib_ph = ph(@formula(Srv(T_Weib, δ) ~ 1), data, Weibull)
weib_aft = aft(@formula(Srv(T_Weib, δ) ~ 1), data, Weibull)
exp_ph = ph(@formula(Srv(T_Exp, δ) ~ 1), data, Exponential)
exp_aft = aft(@formula(Srv(T_Exp, δ) ~ 1), data, Exponential)

@test weib_ph isa StatsModels.TableStatisticalModel
@test weib_aft isa StatsModels.TableStatisticalModel
@test exp_ph isa StatsModels.TableStatisticalModel
@test exp_aft isa StatsModels.TableStatisticalModel

@test weib_ph.model isa ParametricPH{Weibull}
@test weib_aft.model isa ParametricAFT{Weibull}
@test exp_ph.model isa ParametricPH{Exponential}
@test exp_aft.model isa ParametricAFT{Exponential}

# test parameters against theoretical
@test (2.0, 5.0) === round.(params(weib_aft.model.baseline)) === round.(params(weib_ph.model.baseline))
@test (5.0, ) === round.(params(exp_aft.model.baseline)) === round.(params(exp_ph.model.baseline))

# test coefficients against R
weib_ph = ph(@formula(Srv(T_Weib, δ) ~ X1 + X2 + X3), data, Weibull)
weib_aft = aft(@formula(Srv(T_Weib, δ) ~ X1 + X2 + X3), data, Weibull)
exp_ph = ph(@formula(Srv(T_Exp, δ) ~ X1 + X2 + X3), data, Exponential)
exp_aft = aft(@formula(Srv(T_Exp, δ) ~ X1 + X2 + X3), data, Exponential)

R"
    library(survival)
    weib = survreg(Surv($T_Weib) ~ X1 + X2 + X3, data = $data, dist = 'weibull')
    exp = survreg(Surv($T_Exp) ~ X1 + X2 + X3, data = $data, dist = 'exponential')
    R_coeff_exp = coefficients(exp)
    R_coeff_weib = coefficients(weib)
    R_scale_weib = weib$scale
";

βₑ = rcopy(R"R_coeff_exp");
βᵥ = rcopy(R"R_coeff_weib");
αᵥ = rcopy(R"R_scale_weib");
@test round.(βᵥ, digits = 6) == round.(coef(weib_aft), digits = 6)
@test round.(βᵥ, digits = 6) == round.(-coef(weib_ph) * αᵥ, digits = 6)
@test round.(βₑ, digits = 6) == round.(coef(exp_aft), digits = 6)
@test round.(βₑ, digits = 6) == -round.(coef(exp_ph), digits = 6)

# test parameters
# test shape
@test round(shape(baseline(weib_aft)), digits = 6) ===
    round(shape(baseline(weib_ph)), digits = 6) ===
    round(1 / αᵥ, digits = 6)

# test scale
@test round(scale(baseline(weib_aft)), digits = 6) ===
    round(scale(baseline(weib_ph)), digits = 6) ===
    round(exp(βᵥ[1]) , digits = 6)

@test round(scale(baseline(exp_aft)), digits = 6) ===
    round(scale(baseline(exp_ph)), digits = 6) ===
    round(exp(βₑ[1]) , digits = 6)

# test predictions

new_X = DataFrame("X1" => rand(Xs, n), "X2" => rand(Xs, n),
    "X3" => rand(Xs, n));
R"
library(survival)
R_lp_w = predict(weib, $new_X, type = 'lp')
R_lp_e = predict(exp, $new_X, type = 'lp')
R_q_w = predict(weib, $new_X, type = 'quantile', p = 0.2)
R_q_e = predict(exp, $new_X, type = 'quantile', p = 0.2)
";

pred_e = predict(exp_aft, new_X);
pred_w = predict(weib_aft, new_X);

# note we return lp not including intercept
@test all(round.(rcopy(R"R_lp_e"), digits = 6) .==
    round.(pred_e.lp .+ coef(exp_aft)[1], digits = 6))
@test all(round.(rcopy(R"R_lp_w"), digits = 6) .==
    round.(pred_w.lp .+ coef(weib_aft)[1], digits = 6))

@test all(round.(cdf.(pred_w.distr, rcopy(R"R_q_w")), digits = 6) .== 0.2)
@test all(round.(cdf.(pred_e.distr, rcopy(R"R_q_e")), digits = 6) .== 0.2)

# FIXME - ADD TESTS FOR MATRIX INTERFACE, FIT INTERFACE, PREDICT ON INTERCEPT ONLY
true
