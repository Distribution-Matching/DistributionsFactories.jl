using Test
using DistributionsFactories
using Distributions

include("test_mean_var.jl")


@testset "Mean Var tests" begin
    @test test_mean_var_FDist()
end

