using Test
using DistributionsFactories
using Distributions

include("test_mean_var.jl")
include("test_numerical_aux_solvers.jl")


@testset "Mean Var tests" begin
    @test test_mean_var_Beta()
    @test test_mean_var_Chisq()
    @test_broken test_mean_var_Erlang()
    @test_broken test_mean_var_Exponential()
    @test test_mean_var_FDist()
end

@testset "solve_beta_ratio" begin
    @test test_solve_beta_ratio_known_values()
    @test test_solve_beta_ratio_roundtrip()
    @test test_solve_beta_ratio_edge_cases()
    @test test_solve_beta_ratio_errors()
end

