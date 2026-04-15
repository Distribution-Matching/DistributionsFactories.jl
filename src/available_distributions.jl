# Registry of distributions by support type and discovery API

const _SUPPORT_REAL = [Normal, TDist, Logistic, Laplace, Gumbel, SymTriangularDist]
const _SUPPORT_POSITIVE = [Gamma, Erlang, Exponential, LogNormal, Weibull, Frechet,
                           Chi, Chisq, Rayleigh, FDist, InverseGamma, Pareto, FoldedNormal]
const _SUPPORT_UNIT = [Beta, Uniform]
const _SUPPORT_INTEGER_NONNEG = [Poisson, NegativeBinomial, Geometric]
const _SUPPORT_INTEGER_BOUNDED = [Binomial, DiscreteUniform]

# --- Support classification from endpoints ---

function _classify_support(lo, hi)
    if lo == -Inf && hi == Inf
        return :real
    elseif lo == 0 && hi == Inf
        return :positive
    elseif lo == 0 && hi == 1
        return :unit
    else
        return :bounded
    end
end

function _candidates_for_support(support::Symbol)
    support === :real     && return copy(_SUPPORT_REAL)
    support === :positive && return copy(_SUPPORT_POSITIVE)
    support === :unit     && return copy(_SUPPORT_UNIT)
    support === :integer_nonneg  && return copy(_SUPPORT_INTEGER_NONNEG)
    support === :integer_bounded && return copy(_SUPPORT_INTEGER_BOUNDED)
    return []
end

# IntervalSets.jl intervals
function _candidates_for(I::AbstractInterval)
    lo, hi = endpoints(I)
    support = _classify_support(lo, hi)
    support === :bounded && return []
    return _candidates_for_support(support)
end

# Distributions.jl RealInterval
function _candidates_for(I::Distributions.RealInterval)
    support = _classify_support(I.lb, I.ub)
    support === :bounded && return []
    return _candidates_for_support(support)
end

# Integer ranges (discrete)
function _candidates_for(r::AbstractUnitRange)
    lo, hi = first(r), last(r)
    if lo == 0 && hi == typemax(eltype(r))
        return copy(_SUPPORT_INTEGER_NONNEG)
    else
        return copy(_SUPPORT_INTEGER_BOUNDED)
    end
end

# --- Convert keyword specifications to (μ̄, σ̄²) ---

function _resolve_mean_var(; mean=nothing, var=nothing, std=nothing,
                            cv=nothing, scv=nothing, second_moment=nothing,
                            kwargs...)
    μ̄ = mean
    σ̄² = var

    if μ̄ !== nothing
        if std !== nothing
            σ̄² === nothing || throw(ArgumentError("Cannot specify both var and std"))
            σ̄² = std^2
        elseif cv !== nothing
            σ̄² === nothing || throw(ArgumentError("Cannot specify both var and cv"))
            σ̄² = (cv * μ̄)^2
        elseif scv !== nothing
            σ̄² === nothing || throw(ArgumentError("Cannot specify both var and scv"))
            σ̄² = scv * μ̄^2
        elseif second_moment !== nothing
            σ̄² === nothing || throw(ArgumentError("Cannot specify both var and second_moment"))
            σ̄² = second_moment - μ̄^2
        end
    end

    return μ̄, σ̄²
end

function _filter_feasible(candidates, μ̄, σ̄²)
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

# --- Public API ---

"""
    available_distributions(support; kwargs...)

List distribution types whose support matches the given interval or range.
Accepts IntervalSets intervals (`0..Inf`), Distributions.jl `RealInterval`
(e.g. from `support(some_dist)`), or integer ranges (`0:10`).

When keyword arguments specifying moments are provided, returns only those
distributions that are feasible for the given specification.

# Supported interval forms
- `0..Inf` or `support(Gamma(1,1))` — positive distributions
- `-Inf..Inf` or `support(Normal())` — real-line distributions
- `0..1` or `support(Beta(2,3))` — unit-interval distributions
- `0:10` — bounded discrete distributions

# Keyword arguments for moment specification
- `mean` and `var` — mean and variance
- `mean` and `std` — mean and standard deviation
- `mean` and `cv` — mean and coefficient of variation
- `mean` and `scv` — mean and squared coefficient of variation
- `mean` and `second_moment` — mean and E[X²]

# Examples
```julia
available_distributions(0..Inf)
available_distributions(0..Inf, mean=5.0, var=3.0)
available_distributions(support(Beta(2,3)), mean=0.5, std=0.2)
available_distributions(0:10, mean=5.0, var=2.0)
```
"""
function available_distributions(support; kwargs...)
    candidates = _candidates_for(support)
    μ̄, σ̄² = _resolve_mean_var(; kwargs...)
    if μ̄ !== nothing && σ̄² !== nothing
        return _filter_feasible(candidates, μ̄, σ̄²)
    end
    return candidates
end

"""
    available_distributions(; kwargs...)

List all distribution types (across all supports) that can be constructed with
the given moment specification. At least `mean` and one dispersion keyword must
be provided.

# Examples
```julia
available_distributions(mean=5.0, var=3.0)
available_distributions(mean=0.5, std=0.2)
```
"""
function available_distributions(; kwargs...)
    μ̄, σ̄² = _resolve_mean_var(; kwargs...)
    (μ̄ !== nothing && σ̄² !== nothing) || throw(ArgumentError(
        "Must provide mean and a dispersion measure (var, std, cv, scv, or second_moment)"))
    all_types = vcat(_SUPPORT_REAL, _SUPPORT_POSITIVE, _SUPPORT_UNIT,
                     _SUPPORT_INTEGER_NONNEG, _SUPPORT_INTEGER_BOUNDED)
    return _filter_feasible(all_types, μ̄, σ̄²)
end
