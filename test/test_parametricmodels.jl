n = 100
Δ = fill(1, n);
X = DataFrame("X" => fill(1, n));
T_Weib = rand(Weibull(2, 5), n);
T_Exp = rand(Exponential(5), n);

Y = Surv(T_Weib, Δ, "right");
weib_ph = fit_PH(X, Y, "Weibull");

Y = Surv(T_Exp, Δ, "right");
exp_ph = fit_PH(X, Y, "Exponential");

Y = Surv(T_Weib, Δ, "right");
weib_aft = fit_AFT(X, Y, "Weibull");

Y = Surv(T_Exp, Δ, "right");
exp_aft = fit_AFT(X, Y, "Exponential");

# test parameters against theoretical
@assert (2.0, 5.0) === round.(params(weib_aft.baseline)) === round.(params(weib_ph.baseline))
@assert (5.0, ) === round.(params(exp_aft.baseline)) === round.(params(exp_ph.baseline))

# test coefficients against R

Xs = Uniform();
X = DataFrame("X1" => rand(Xs, n), "X2" => rand(Xs, n),
    "X3" => rand(Xs, n));

Y = Surv(T_Weib, Δ, "right");
weib_ph = fit_PH(X, Y, "Weibull");

Y = Surv(T_Exp, Δ, "right");
exp_ph = fit_PH(X, Y, "Exponential");

Y = Surv(T_Weib, Δ, "right");
weib_aft = fit_AFT(X, Y, "Weibull");

Y = Surv(T_Exp, Δ, "right");
exp_aft = fit_AFT(X, Y, "Exponential");

R"
library(survival)
weib = survreg(Surv($T_Weib) ~ ., data = $X, dist = 'weibull')
exp = survreg(Surv($T_Exp) ~ ., data = $X, dist = 'exponential')
R_coeff_exp = coefficients(exp)
R_coeff_weib = coefficients(weib)
R_scale_weib = weib$scale
";

βₑ = rcopy(R"R_coeff_exp");
βᵥ = rcopy(R"R_coeff_weib");
αᵥ = rcopy(R"R_scale_weib");
@assert round.(βᵥ, digits = 6) == round.(weib_aft.coefficients.betas, digits = 6)
@assert round.(βᵥ, digits = 6) == round.(-weib_ph.coefficients.betas * αᵥ, digits = 6)
@assert round.(βₑ, digits = 6) == round.(exp_aft.coefficients.betas, digits = 6)
@assert round.(βₑ, digits = 6) == -round.(exp_ph.coefficients.betas, digits = 6)

# test parameters
# test shape
@assert round(shape(weib_aft.baseline), digits = 6) ===
    round(shape(weib_ph.baseline), digits = 6) ===
    round(1 / αᵥ, digits = 6)

# test scale
@assert round(scale(weib_aft.baseline), digits = 6) ===
    round(scale(weib_ph.baseline), digits = 6) ===
    round(exp(βᵥ[1]) , digits = 6)

@assert round(scale(exp_aft.baseline), digits = 6) ===
    round(scale(exp_ph.baseline), digits = 6) ===
    round(exp(βₑ[1]) , digits = 6)

# test predictions

new_X = DataFrame("X1" => rand(Xs, n), "X2" => rand(Xs, n),
    "X3" => rand(Xs, n));
R"
library(survival)
weib = survreg(Surv($T_Weib) ~ ., data = $X, dist = 'weibull')
exp = survreg(Surv($T_Exp) ~ ., data = $X, dist = 'exponential')
R_lp_w = predict(weib, $new_X, type = 'lp')
R_lp_e = predict(exp, $new_X, type = 'lp')
R_q_w = predict(weib, $new_X, type = 'quantile', p = 0.2)
R_q_e = predict(exp, $new_X, type = 'quantile', p = 0.2)
";

pred_e = predict_Parametric(exp_aft, new_X);
pred_w = predict_Parametric(weib_aft, new_X);

# note we return lp not including intercept
@assert all(round.(rcopy(R"R_lp_e"), digits = 6) .==
round.(pred_e.lp .+ exp_aft.coefficients.betas[1], digits = 6))
@assert all(round.(rcopy(R"R_lp_w"), digits = 6) .==
round.(pred_w.lp .+ weib_aft.coefficients.betas[1], digits = 6))

@assert all(round.(cdf.(predict_Parametric(weib_aft, new_X).distr, rcopy(R"R_q_w")), digits = 6) .== 0.2)
@assert all(round.(cdf.(predict_Parametric(exp_aft, new_X).distr, rcopy(R"R_q_e")), digits = 6) .== 0.2)
