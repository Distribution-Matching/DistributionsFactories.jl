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
using Random

# Distribution types not provided by Distributions.jl. Real
# subtypes of `ContinuousUnivariateDistribution` / `DiscreteUnivariateDistribution`
# with the full Distributions.jl interface — kept under `extensions/` so it's
# clear which types we add and which are upstream.
include("extensions/folded_normal.jl")
include("extensions/discrete_symmetric_triangular.jl")
include("extensions/discrete_triangular.jl")

include("solvers.jl")
include("langevin.jl")
include("exists_dist.jl")
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
        DistSpec,
        fixed_params,
        free_params,
        ..,
        FoldedNormal,
        DiscreteTriangular,
        DiscreteSymmetricTriangular

end # end the module
