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

struct ModeSpec
    mode::Float64
end

struct MeanModeSpec
    μ̄::Float64
    mode::Float64
end

struct ModeVarSpec
    mode::Float64
    σ̄²::Float64
end

struct ModeQuantileSpec
    mode::Float64
    p::Float64
    q::Float64
end

struct ModeIQRSpec
    mode::Float64
    iqr::Float64
end

struct MeanVarModeSpec
    μ̄::Float64
    σ̄²::Float64
    mode::Float64
end

# --- Spec parser (the only "router") ---

function _moment_spec(; mean=nothing, var=nothing, std=nothing,
                       cv=nothing, scv=nothing, second_moment=nothing,
                       median=nothing, q1=nothing, q3=nothing, iqr=nothing,
                       quantiles=nothing, mode=nothing)
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

    # Mode-based specs
    if mode !== nothing && mean !== nothing && σ̄² !== nothing
        return MeanVarModeSpec(mean, σ̄², mode)
    end
    if mode !== nothing && mean !== nothing && σ̄² === nothing
        return MeanModeSpec(mean, mode)
    end
    if mode !== nothing && σ̄² !== nothing
        return ModeVarSpec(mode, σ̄²)
    end
    if mode !== nothing && median !== nothing
        return ModeQuantileSpec(mode, 0.5, median)
    end
    if mode !== nothing && q1 !== nothing
        return ModeQuantileSpec(mode, 0.25, q1)
    end
    if mode !== nothing && q3 !== nothing
        return ModeQuantileSpec(mode, 0.75, q3)
    end
    if mode !== nothing && iqr !== nothing
        return ModeIQRSpec(mode, iqr)
    end
    if mode !== nothing && mean === nothing && σ̄² === nothing
        return ModeSpec(mode)
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

Construct a distribution of type `D` (or `DistSpec` via `@dist`) from the given
moment or quantile specification. All parameters are keyword arguments.

# Moment specifications
- `mean` and `var` — mean and variance
- `mean` and `std` — mean and standard deviation
- `mean` and `cv` — mean and coefficient of variation
- `mean` and `scv` — mean and squared coefficient of variation
- `mean` and `second_moment` — mean and E[X²]
- `mean` alone — for 1-parameter distributions or `DistSpec`
- `var` alone — for 1-parameter distributions or `DistSpec`

# Mode-based specifications (for distributions parameterised by a mode)
- `mode` alone — for 1-parameter mode-fixed distributions (e.g. `Rayleigh`)
- `mean` and `mode` — for 2-parameter distributions (e.g. `Gamma`, `Beta`)
- `mode` and `var` — for 2-parameter distributions
- `mean`, `var`, and `mode` — for 3-parameter distributions like `TriangularDist`
  and `DiscreteTriangular`
- `mode` and `iqr` — for symmetric 2-parameter distributions
- `mode` and `median` (or `q1`, `q3`) — for 2-parameter distributions

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
make_dist(TriangularDist, mean=5.0, var=2.0, mode=4.0)
make_dist(DiscreteTriangular, mean=5.0, var=2.0, mode=5)
make_dist(FoldedNormal, mean=2.5, var=1.2)
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

# --- Dispatch on spec type: mode ---

_make_dist(D, s::ModeSpec) = dist_from_mode(D, s.mode)
_make_dist(D, s::MeanModeSpec) = dist_from_mean_mode(D, s.μ̄, s.mode)
_make_dist(D, s::ModeVarSpec) = dist_from_mode_var(D, s.mode, s.σ̄²)
_make_dist(D, s::ModeQuantileSpec) = dist_from_mode_quantile(D, s.mode, s.p, s.q)
_make_dist(D, s::ModeIQRSpec) = dist_from_mode_iqr(D, s.mode, s.iqr)
_make_dist(D, s::MeanVarModeSpec) = dist_from_mean_var_mode(D, s.μ̄, s.σ̄², s.mode)

# --- Dispatch on spec type: quantiles ---

_make_dist(D, s::QuantileSpec) = dist_from_quantile(D, s.p, s.q)
_make_dist(D, s::TwoQuantileSpec) = dist_from_quantiles(D, s.p1, s.q1, s.p2, s.q2)
_make_dist(D, s::MeanQuantileSpec) = dist_from_mean_quantile(D, s.μ̄, s.p, s.q)

# --- Dispatch on spec type: DistSpec ---

_make_dist(p::DistSpec, s::MeanVarSpec) = dist_from_mean_var(p, s.μ̄, s.σ̄²)
_make_dist(p::DistSpec, s::MeanSpec) = dist_from_mean(p, s.μ̄)
_make_dist(p::DistSpec, s::VarSpec) = dist_from_var(p, s.σ̄²)

# --- With support ---

function _make_dist_on_support(D, s::MeanVarSpec, support)
    return dist_from_mean_var_on_support(D, s.μ̄, s.σ̄², support=support)
end

# --- dist_exists ---

"""
    dist_exists(D; mean, var, std, cv, scv, second_moment, support) -> Bool

Pure predicate: can a distribution of type `D` be constructed with the given
moment specification, optionally on a given `support`? Returns `true` or
`false`; never throws for infeasible moments. (Use [`make_dist`](@ref) to
actually build the distribution — it throws `DomainError` with a reason when
this predicate would have returned `false`.)

# Examples
```julia
dist_exists(Beta, mean=0.5, var=0.1)                    # true
dist_exists(Exponential, mean=2.5, var=1.5)              # false (needs σ̄² = μ̄²)
dist_exists(Beta, mean=3.5, var=0.5, support=2..7)       # true
dist_exists(Pareto, mean=3.0, support=5..Inf)            # false (μ_std ≤ 0)
```
"""
function dist_exists(D; support=nothing, kwargs...)::Bool
    spec = _moment_spec(; kwargs...)
    if support === nothing
        return _dist_exists(D, spec)
    else
        return _dist_exists_on_support(D, spec, support)
    end
end

_dist_exists(D, s::MeanVarSpec)::Bool = exists_dist_from_mean_var(D, s.μ̄, s.σ̄²)
_dist_exists(D, s::MeanSpec)::Bool = true
_dist_exists(D, s::VarSpec)::Bool = true

# Mode- and quantile-augmented specs: there is no generic feasibility predicate,
# so try the constructor and return whether it succeeds. This matches the spirit
# of the existing non-throwing predicate (`true`/`false`, never errors).
_try_construct(f, args...) = try; f(args...); true; catch; false; end

_dist_exists(D, s::MeanVarModeSpec)::Bool = _try_construct(dist_from_mean_var_mode, D, s.μ̄, s.σ̄², s.mode)
_dist_exists(D, s::MeanModeSpec)::Bool    = _try_construct(dist_from_mean_mode, D, s.μ̄, s.mode)
_dist_exists(D, s::ModeVarSpec)::Bool     = _try_construct(dist_from_mode_var, D, s.mode, s.σ̄²)
_dist_exists(D, s::ModeQuantileSpec)::Bool = _try_construct(dist_from_mode_quantile, D, s.mode, s.p, s.q)
_dist_exists(D, s::ModeIQRSpec)::Bool     = _try_construct(dist_from_mode_iqr, D, s.mode, s.iqr)
_dist_exists(D, s::ModeSpec)::Bool        = _try_construct(dist_from_mode, D, s.mode)
_dist_exists(D, s::QuantileSpec)::Bool    = _try_construct(dist_from_quantile, D, s.p, s.q)
_dist_exists(D, s::TwoQuantileSpec)::Bool = _try_construct(dist_from_quantiles, D, s.p1, s.q1, s.p2, s.q2)
_dist_exists(D, s::MeanQuantileSpec)::Bool = _try_construct(dist_from_mean_quantile, D, s.μ̄, s.p, s.q)

function _dist_exists_on_support(D, s::MeanVarSpec, support)::Bool
    lo, hi = _support_endpoints(support)
    μ_std, σ²_std = _moments_to_standard(D, s.μ̄, s.σ̄², lo, hi)
    return exists_dist_from_mean_var(D, μ_std, σ²_std)
end

function _dist_exists_on_support(D, s::MeanSpec, support)::Bool
    lo, hi = _support_endpoints(support)
    μ_std, _ = _moments_to_standard(D, s.μ̄, 1.0, lo, hi)
    if μ_std ≤ 0 && _natural_support(D) === :positive
        return false
    end
    return true
end

function _support_endpoints(support)
    if support isa AbstractInterval
        return endpoints(support)
    elseif support isa Distributions.RealInterval
        return (support.lb, support.ub)
    elseif support isa AbstractUnitRange
        return (Float64(first(support)), Float64(last(support)))
    else
        throw(ArgumentError("Unsupported support type: $(typeof(support))"))
    end
end



# Registry of distributions by support type and discovery API

const _SUPPORT_REAL = [Normal, TDist, Logistic, Laplace, Gumbel, SymTriangularDist]
const _SUPPORT_POSITIVE = [Gamma, Erlang, Exponential, LogNormal, Weibull, Frechet,
                           Chi, Chisq, Rayleigh, FDist, InverseGamma, Pareto, FoldedNormal]
const _SUPPORT_UNIT = [Beta, Uniform]
const _SUPPORT_INTEGER_NONNEG = [Poisson, NegativeBinomial, Geometric]
const _SUPPORT_INTEGER_BOUNDED = [Binomial, DiscreteUniform, DiscreteSymmetricTriangular]

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

_filter_feasible(candidates, μ̄, σ̄²) =
    filter(D -> exists_dist_from_mean_var(D, μ̄, σ̄²), candidates)

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
