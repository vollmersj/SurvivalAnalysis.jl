Surv(1.0, false, "left")
Surv(1.0, true, "left")
Surv(1.0, false, "right")
Surv(1.0, false, "interval")
Surv(1.0, 2.0)

using Distributions

Surv.(rand(Uniform(), 10), rand(Binomial(1, 0.5), 10) .== 1, "right")
Surv.(rand(Uniform(), 10), rand(Binomial(1, 0.5), 10) .== 1, "left")
Surv.(rand(Uniform(), 10), rand(Uniform(), 10))

survs = Surv.(round.(rand(Uniform(1, 5), 10)), rand(Binomial(1, 0.5), 10) .== 1, "right")
ot = OutcomeTimes(survs)
ot.times
uniqueTimes(ot)
