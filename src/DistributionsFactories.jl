"""
    DistributionsFactories

A Julia package for constructing probability distributions from moment specifications
(mean, variance) and quantile specifications. Builds on
[Distributions.jl](https://github.com/JuliaStats/Distributions.jl).

The main entry point is [`dist_from_mean_var`](@ref), which constructs a distribution
of a given type from a target mean and variance. Additional interfaces include
[`dist_from_mean_std`](@ref), [`dist_from_mean_cv`](@ref), [`dist_from_quantile`](@ref),
and [`dist_from_mean`](@ref).
"""
module DistributionsFactories

using Distributions
using SpecialFunctions
using Polynomials
using Roots
using QuadGK
using IntervalSets

# Marker types for distributions not in Distributions.jl
abstract type AbstractDistributionsFactoriesType <: Distribution{Univariate, Continuous} end

struct FoldedNormal <: AbstractDistributionsFactoriesType end

abstract type AbstractDistributionsFactoriesDiscreteType <: Distribution{Univariate, Discrete} end

struct DiscreteTriangular <: AbstractDistributionsFactoriesDiscreteType end
struct DiscreteSymmetricTriangular <: AbstractDistributionsFactoriesDiscreteType end
struct TruncatedPoisson <: AbstractDistributionsFactoriesDiscreteType end

include("numerical_aux_solvers.jl")
include("exists_dist_from_mean_var.jl")
include("dist_from_mean_var.jl")
include("dist_from_variants.jl")
include("dist_from_mean.jl")
include("dist_from_quantile.jl")
include("dist_on_support.jl")
include("partial_dist.jl")
include("dist_macro.jl")
include("available_distributions.jl")

export dist_from_mean_var_on_support,
        dist_from_mean_var,
        dist_from_mean_std,
        dist_from_mean_cv,
        dist_from_mean_scv,
        dist_from_mean_second_moment,
        dist_from_var,
        dist_from_std,
        dist_from_mean,
        dist_from_quantile,
        dist_from_quantiles,
        dist_from_median,
        dist_from_q1,
        dist_from_q3,
        dist_from_median_iqr,
        dist_from_q1_q3,
        dist_from_mean_quantile,
        dist_from_mean_median,
        exists_dist_from_mean_var,
        exists_dist_from_mean_std,
        exists_dist_from_mean_cv,
        exists_dist_from_mean_scv,
        exists_dist_from_mean_second_moment,
        available_distributions,
        PartialDist,
        @dist,
        ..,
        FoldedNormal,
        DiscreteTriangular,
        DiscreteSymmetricTriangular,
        TruncatedPoisson

end # end the module