using Test
using DistributionsFactories
using Distributions

include("test_mean_var.jl")
include("test_mean_variants.jl")
include("test_quantile.jl")
include("test_numerical_aux_solvers.jl")
include("test_available_distributions.jl")
include("test_dist_on_support.jl")
include("test_dist_macro.jl")


@testset "Mean Var tests" begin
    @test test_mean_var_Beta()
    @test test_mean_var_Normal()
    @test test_mean_var_Uniform()
    @test test_mean_var_Logistic()
    @test test_mean_var_Laplace()
    @test test_mean_var_LogNormal()
    @test test_mean_var_Chisq()
    @test test_mean_var_Erlang()
    @test test_mean_var_Exponential()
    @test test_mean_var_Gamma()
    @test test_mean_var_Gumbel()
    @test test_mean_var_InverseGamma()
    @test test_mean_var_Poisson()
    @test test_mean_var_NegativeBinomial()
    @test test_mean_var_Binomial()
    @test_broken test_mean_var_Weibull_simple_bracketing()
    @test_broken test_mean_var_Weibull_oscar_garcia()
    @test_broken test_mean_var_Weibull_methods_agree()
    @test test_mean_var_TDist_instance()
    @test test_mean_var_TDist_instance_returns_affine()
    @test test_mean_var_TDist_instance_low_dof_errors()
    @test test_mean_var_FDist()
    @test test_mean_var_Geometric()
    @test test_mean_var_Pareto()
    @test test_mean_var_SymTriangularDist()
    @test test_mean_var_DiscreteUniform()
    @test test_mean_var_Cauchy_throws()
    @test test_mean_var_feasibility_rejections()
    @test test_mean_var_Geometric_dist_from_mean()
    @test test_dist_from_mean()
    @test test_dist_from_var_and_std()
    @test test_dist_from_mean_second_moment()
end

@testset "Mean variant wrappers" begin
    @test test_mean_std_roundtrip()
    @test test_mean_cv_roundtrip()
    @test test_mean_scv_roundtrip()
    @test test_variants_consistency()
    @test test_exists_variants()
end

@testset "Quantile-based construction" begin
    @test test_quantile_exponential()
    @test test_median_exponential()
    @test test_q1_q3_exponential()
    @test test_quantiles_location_scale()
    @test test_quantiles_gamma()
    @test test_quantiles_beta()
    @test test_mean_quantile_gamma()
    @test test_mean_quantile_beta()
    @test test_median_iqr_normal()
    @test test_mean_median_gamma()
end

@testset "available_distributions" begin
    @test test_available_by_interval()
    @test test_available_by_realinterval()
    @test test_available_by_range()
    @test test_available_with_mean_var()
    @test test_available_with_mean_std()
    @test test_available_with_mean_cv()
    @test test_available_kwargs_only()
    @test test_available_kwargs_requires_dispersion()
end

@testset "dist_on_support" begin
    @test test_support_beta_scaled()
    @test test_support_uniform_scaled()
    @test test_support_gamma_shifted()
    @test test_support_gamma_flipped()
    @test test_support_binomial_shifted()
    @test test_support_discrete_uniform_shifted()
    @test test_support_normal_truncated()
    @test test_support_gamma_truncated()
    @test test_support_errors()
end

@testset "@dist macro and PartialDist" begin
    @test test_dist_macro_bare_type()
    @test test_dist_macro_full_instance()
    @test test_dist_macro_partial()
    @test test_partial_dist_from_mean()
    @test test_partial_normal_fix_sigma()
    @test test_partial_dist_from_mean_var()
    @test test_partial_all_missing_delegates()
    @test test_partial_beta()
    @test test_instance_tdist()
    @test test_type_via_macro()
end

@testset "solve_beta_ratio" begin
    @test test_solve_beta_ratio_known_values()
    @test test_solve_beta_ratio_roundtrip()
    @test test_solve_beta_ratio_edge_cases()
    @test test_solve_beta_ratio_errors()
end

