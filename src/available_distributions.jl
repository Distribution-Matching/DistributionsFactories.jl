# Registry of distributions by support type

const _SUPPORT_DISTRIBUTIONS = Dict{Symbol, Vector{Any}}(
    :real => [Normal, TDist, Logistic, Laplace, Gumbel, SymTriangularDist],
    :positive => [Gamma, Erlang, Exponential, LogNormal, Weibull, Frechet,
                  Chi, Chisq, Rayleigh, FDist, InverseGamma, Pareto, FoldedNormal],
    :unit => [Beta, Uniform],
    :integer_nonneg => [Poisson, NegativeBinomial, Geometric],
    :integer_bounded => [Binomial, DiscreteUniform],
)

const _VALID_SUPPORTS = keys(_SUPPORT_DISTRIBUTIONS)

"""
    available_distributions(support::Symbol)

List the distribution types available for a given support.

Supported values for `support`:
- `:real` — distributions on (-∞, ∞)
- `:positive` — distributions on [0, ∞)
- `:unit` — distributions on [0, 1]
- `:integer_nonneg` — discrete distributions on {0, 1, 2, …}
- `:integer_bounded` — discrete distributions on {0, …, n}

Returns a `Vector` of distribution types.

# Example
```julia
available_distributions(:positive)
# [Gamma, Erlang, Exponential, LogNormal, Weibull, ...]
```
"""
function available_distributions(support::Symbol)
    support ∈ _VALID_SUPPORTS || throw(ArgumentError(
        "Unknown support :$support. Valid options: $(join(_VALID_SUPPORTS, ", ", " and "))"))
    return copy(_SUPPORT_DISTRIBUTIONS[support])
end

"""
    available_distributions(support::Symbol, μ̄, σ̄²)

List the distribution types on the given `support` that can be constructed with
mean `μ̄` and variance `σ̄²`. Tries each candidate via `exists_dist_from_mean_var`
and returns only those that are feasible.

Returns a `Vector` of distribution types.

# Example
```julia
available_distributions(:positive, 5.0, 3.0)
# [Gamma, LogNormal, Weibull, Frechet, InverseGamma]
```
"""
function available_distributions(support::Symbol, μ̄::Number, σ̄²::Number)
    candidates = available_distributions(support)
    feasible = []
    for D in candidates
        try
            exists_dist_from_mean_var(D, μ̄, σ̄²)
            push!(feasible, D)
        catch
        end
    end
    return feasible
end

"""
    available_distributions(μ̄, σ̄²)

List all distribution types (across all supports) that can be constructed with
mean `μ̄` and variance `σ̄²`.

Returns a `Vector` of distribution types.

# Example
```julia
available_distributions(0.5, 0.05)
# [Normal, Logistic, Laplace, ..., Gamma, ..., Beta, Uniform]
```
"""
function available_distributions(μ̄::Number, σ̄²::Number)
    feasible = []
    for support in _VALID_SUPPORTS
        append!(feasible, available_distributions(support, μ̄, σ̄²))
    end
    return feasible
end
