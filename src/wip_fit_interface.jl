using Survival, DataFrames, Distributions
n = 100
Δ = fill(1, n);
X = Matrix(DataFrame("X" => fill(1, n)));
T_Weib = rand(Weibull(2, 5), n);
T_Exp = rand(Exponential(5), n);
Y = Surv(T_Exp, Δ, "right");
weib_ph = ph(X, Y, Exponential);
data = DataFrame("X" => fill(1, n), "Y" => Y)
ph(@formula(Y ~ X), data, Exponential) ## FIXME
Table
