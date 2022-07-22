abstract type FittedParametric <: SurvivalModel end

coef(obj::FittedParametric) = obj.coefficients
scale(obj::FittedParametric) = scale(obj.distribution)
shape(obj::FittedParametric) = shape(obj.distribution)
