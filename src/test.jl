include("src/Surv.jl")

using Random: seed!
seed!(1)
et = Surv.(round.(rand(Uniform(1, 10), 10)), rand(Binomial(), 10) .== 1, "right")

# km = kaplan(et)
# confint.(Ref(km), km.times)

na = nelson(et)
confint.(Ref(na), na.times)

plot(na)
