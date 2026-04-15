# Unified API: make_dist(D; kwargs...) and dist_exists(D; kwargs...)
#
# Typed moment/quantile specs enable multiple dispatch over different
# keyword combinations, keeping the router minimal.

# --- Spec types ---

struct MeanVarSpec
    μ̄::Float64
    σ̄²::Float64
end

struct MeanSpec
    μ̄::Float64
end

struct VarSpec
    σ̄²::Float64
end

struct QuantileSpec
    p::Float64
    q::Float64
end

struct TwoQuantileSpec
    p1::Float64
    q1::Float64
    p2::Float64
    q2::Float64
end

struct MeanQuantileSpec
    μ̄::Float64
    p::Float64
    q::Float64
end

# --- Spec parser (the only "router") ---

function _moment_spec(; mean=nothing, var=nothing, std=nothing,
                       cv=nothing, scv=nothing, second_moment=nothing,
                       median=nothing, q1=nothing, q3=nothing, iqr=nothing,
                       quantiles=nothing)
    # Resolve variance from alternative dispersion measures
    σ̄² = var
    if σ̄² === nothing && std !== nothing
        σ̄² = std^2
    elseif σ̄² === nothing && mean !== nothing && cv !== nothing
        σ̄² = (cv * mean)^2
    elseif σ̄² === nothing && mean !== nothing && scv !== nothing
        σ̄² = scv * mean^2
    elseif σ̄² === nothing && mean !== nothing && second_moment !== nothing
        σ̄² = second_moment - mean^2
    end

    # Quantile-based specs
    if quantiles !== nothing
        length(quantiles) == 2 || throw(ArgumentError("quantiles must be a vector of 2 (p, q) tuples"))
        (p1, q1_val), (p2, q2_val) = quantiles
        return TwoQuantileSpec(p1, q1_val, p2, q2_val)
    end
    if q1 !== nothing && q3 !== nothing
        return TwoQuantileSpec(0.25, q1, 0.75, q3)
    end
    if median !== nothing && iqr !== nothing
        return TwoQuantileSpec(0.25, median - iqr/2, 0.75, median + iqr/2)
    end
    if mean !== nothing && median !== nothing
        return MeanQuantileSpec(mean, 0.5, median)
    end
    if median !== nothing
        return QuantileSpec(0.5, median)
    end
    if q1 !== nothing
        return QuantileSpec(0.25, q1)
    end
    if q3 !== nothing
        return QuantileSpec(0.75, q3)
    end

    # Moment-based specs
    if mean !== nothing && σ̄² !== nothing
        return MeanVarSpec(mean, σ̄²)
    elseif mean !== nothing
        return MeanSpec(mean)
    elseif σ̄² !== nothing
        return VarSpec(σ̄²)
    end

    throw(ArgumentError("Must provide at least one moment or quantile specification"))
end

# --- make_dist: entry point ---

"""
    make_dist(D; mean, var, std, cv, scv, second_moment, median, q1, q3, iqr, quantiles, support)

Construct a distribution of type `D` (or `PartialDist` via `@dist`) from the given
moment or quantile specification. All parameters are keyword arguments.

# Moment specifications
- `mean` and `var` — mean and variance
- `mean` and `std` — mean and standard deviation
- `mean` and `cv` — mean and coefficient of variation
- `mean` and `scv` — mean and squared coefficient of variation
- `mean` and `second_moment` — mean and E[X²]
- `mean` alone — for 1-parameter distributions or `PartialDist`
- `var` alone — for 1-parameter distributions or `PartialDist`

# Quantile specifications
- `median` — single quantile (1-parameter distributions)
- `q1` and `q3` — first and third quartiles
- `median` and `iqr` — median and interquartile range
- `mean` and `median` — mean and median
- `quantiles=[(p1,q1), (p2,q2)]` — two arbitrary quantiles

# Support
- `support=a..b` — place distribution on arbitrary domain (IntervalSets or RealInterval)
- `support=a:b` — discrete shift

# Examples
```julia
make_dist(Gamma, mean=5.0, var=3.0)
make_dist(Gamma, mean=5.0, cv=0.5)
make_dist(Exponential, mean=3.0)
make_dist(Exponential, median=2.0)
make_dist(Normal, q1=10.0, q3=30.0)
make_dist(Beta, mean=0.4, median=0.35)
make_dist(Beta, mean=3.5, var=0.5, support=2..7)
make_dist(@dist(Gamma(3.0, _)), mean=5.0)
make_dist(@dist(Logistic(2.0, _)), var=22.3)
```
"""
function make_dist(D; support=nothing, kwargs...)
    spec = _moment_spec(; kwargs...)
    if support === nothing
        return _make_dist(D, spec)
    else
        return _make_dist_on_support(D, spec, support)
    end
end

# --- Dispatch on spec type: moments ---

_make_dist(D, s::MeanVarSpec) = dist_from_mean_var(D, s.μ̄, s.σ̄²)
_make_dist(D, s::MeanSpec) = dist_from_mean(D, s.μ̄)
_make_dist(D, s::VarSpec) = dist_from_var(D, s.σ̄²)

# --- Dispatch on spec type: quantiles ---

_make_dist(D, s::QuantileSpec) = dist_from_quantile(D, s.p, s.q)
_make_dist(D, s::TwoQuantileSpec) = dist_from_quantiles(D, s.p1, s.q1, s.p2, s.q2)
_make_dist(D, s::MeanQuantileSpec) = dist_from_mean_quantile(D, s.μ̄, s.p, s.q)

# --- Dispatch on spec type: PartialDist ---

_make_dist(p::PartialDist, s::MeanVarSpec) = dist_from_mean_var(p, s.μ̄, s.σ̄²)
_make_dist(p::PartialDist, s::MeanSpec) = dist_from_mean(p, s.μ̄)
_make_dist(p::PartialDist, s::VarSpec) = dist_from_var(p, s.σ̄²)

# --- With support ---

function _make_dist_on_support(D, s::MeanVarSpec, support)
    return dist_from_mean_var_on_support(D, s.μ̄, s.σ̄², support=support)
end

# --- dist_exists ---

"""
    dist_exists(D; mean, var, std, cv, scv, second_moment)

Check whether a distribution of type `D` can be constructed with the given
moment specification. Returns `true` or throws `DomainError`.

# Examples
```julia
dist_exists(Beta, mean=0.5, var=0.1)    # true
dist_exists(Exponential, mean=2.5, var=1.5)  # throws DomainError
```
"""
function dist_exists(D; kwargs...)
    spec = _moment_spec(; kwargs...)
    return _dist_exists(D, spec)
end

_dist_exists(D, s::MeanVarSpec) = exists_dist_from_mean_var(D, s.μ̄, s.σ̄²)
_dist_exists(D, s::MeanSpec) = true  # 1-param distributions always exist if mean is valid
_dist_exists(D, s::VarSpec) = true
