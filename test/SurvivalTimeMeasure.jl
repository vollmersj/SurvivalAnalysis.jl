n = 20;
δ = rand(Bernoulli(), n);
t = randn(n);
Y = Surv(t, δ, :r);
t̂ = randn(n);

@testset "MSE" begin
    l = MSE(Y, t̂)
    @test l isa MSE
    @test l.mean == mean(l.losses)
    @test l.losses == (t[δ] - t̂[δ]).^2
    @test l.se == std(l.losses) / sqrt(sum(δ))
end

@testset "MAE" begin
    l = MAE(Y, t̂)
    @test l isa MAE
    @test l.mean == mean(l.losses)
    @test l.losses == abs.(t[δ] - t̂[δ])
    @test l.se == std(l.losses) / sqrt(sum(δ))
end

@testset "RMSE" begin
    l = RMSE(Y, t̂)
    @test l isa RMSE
    @test l.mean == sqrt(mean(l.losses))
    @test l.losses == (t[δ] - t̂[δ]).^2
    @test l.se == (std(l.losses) / sqrt(sum(δ))) / (2 * sqrt(mean(l.losses)))
end
