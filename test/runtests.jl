using Test
using DistributionsFactories
using Distributions

include("test_mean_var.jl")


@testset "Mean Var tests" begin
    @test test_mean_var_Beta()
    @test test_mean_var_Chisq()
    @test test_mean_var_Erlang()
    @test test_mean_var_Exponential()
    @test test_mean_var_FDist()
end

