"""
    DistributionsFactories

A Julia package for constructing probability distributions from moment specifications
(mean, variance) and quantile specifications, on standard or arbitrary supports.
Builds on [Distributions.jl](https://github.com/JuliaStats/Distributions.jl).

The primary entry point is [`make_dist`](@ref). See also [`dist_exists`](@ref),
[`available_distributions`](@ref), and the [`@dist`](@ref) macro.
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

include("solvers.jl")
include("moment_matching.jl")
include("quantile_matching.jl")
include("partial_dist.jl")
include("support.jl")
include("api.jl")

# Public API
export make_dist,
        dist_exists,
        available_distributions,
        @dist,
        PartialDist,
        fixed_params,
        free_params,
        ..,
        FoldedNormal,
        DiscreteTriangular,
        DiscreteSymmetricTriangular,
        TruncatedPoisson

end # end the module
