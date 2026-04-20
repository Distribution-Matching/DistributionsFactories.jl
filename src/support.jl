# Construct distributions on arbitrary supports via affine transform or truncation.
#
# Decision logic:
#   - If the requested support has the same "shape" as the natural support
#     (e.g. [a,∞) for a [0,∞) distribution), use an affine transform.
#   - If the requested support is strictly contained in the natural support
#     (e.g. [0,1] for a (-∞,∞) distribution), use truncation.

# --- Natural support classification ---

_natural_support(::Type{<:Union{Normal,TDist,Logistic,Laplace,Gumbel,SymTriangularDist}}) = :real
_natural_support(::Type{<:Union{Gamma,Erlang,Exponential,LogNormal,Weibull,Frechet,
                                Chi,Chisq,Rayleigh,FDist,InverseGamma,Pareto}}) = :positive
_natural_support(::Type{<:FoldedNormal}) = :positive
_natural_support(::Type{<:Union{Beta,Uniform}}) = :unit
_natural_support(::Type{<:Union{Binomial,DiscreteUniform,DiscreteSymmetricTriangular,DiscreteTriangular}}) = :integer_bounded
_natural_support(::Type{<:Union{Poisson,NegativeBinomial,Geometric}}) = :integer_nonneg

# DistSpec delegates to the underlying distribution type
_natural_support(::DistSpec{D}) where {D} = _natural_support(D)

function _requested_support_shape(lo, hi)
    if lo == -Inf && hi == Inf
        return :real
    elseif isfinite(lo) && hi == Inf
        return :half_right    # [a, ∞)
    elseif lo == -Inf && isfinite(hi)
        return :half_left     # (-∞, b]
    elseif isfinite(lo) && isfinite(hi)
        return :bounded       # [a, b]
    end
    throw(ArgumentError("Invalid support endpoints: ($lo, $hi)"))
end

# --- Public entry point ---

"""
    dist_from_mean_var_on_support(D, μ̄, σ̄²; support)

Construct a distribution of type `D` with mean `μ̄` and variance `σ̄²` on the given
`support`. Accepts IntervalSets intervals (`support=2..7`), Distributions.jl
`RealInterval` (e.g. `support=support(some_dist)`), or integer ranges (`support=10:15`).

The function automatically determines whether to use an **affine transform** or
**truncation** based on the relationship between `D`'s natural support and the
requested support:

- **Affine**: when the requested support has the same shape as the natural one
  (e.g. Beta on `[2,7]`, Gamma on `[3,∞)`). Returns a `LocationScale` wrapper.
- **Truncation**: when the requested support is strictly contained in the natural
  one (e.g. Normal on `[0,1]`, Gamma on `[0,10]`). Returns a `Truncated` wrapper.

# Examples
```julia
dist_from_mean_var(Beta, 3.5, 0.5, support=2..7)
dist_from_mean_var(Gamma, 8.0, 3.0, support=3..Inf)
dist_from_mean_var(Gamma, 5.0, 3.0, support=-Inf..8)
dist_from_mean_var(Normal, 0.5, 0.04, support=0..1)
dist_from_mean_var(Binomial, 12.0, 1.2, support=10:15)
```
"""
function dist_from_mean_var_on_support(D, μ̄, σ̄²; support)
    if support isa AbstractInterval
        lo, hi = endpoints(support)
        return _dist_on_support(D, μ̄, σ̄², lo, hi)
    elseif support isa Distributions.RealInterval
        return _dist_on_support(D, μ̄, σ̄², support.lb, support.ub)
    elseif support isa AbstractUnitRange
        return _dist_on_support_discrete(D, μ̄, σ̄², support)
    else
        throw(ArgumentError("Unsupported `support` type: $(typeof(support)). Use an interval (a..b), RealInterval, or integer range (a:b)."))
    end
end

# --- Continuous support logic ---

function _dist_on_support(D, μ̄, σ̄², lo, hi)
    natural = _natural_support(D)
    requested = _requested_support_shape(lo, hi)

    # --- Real-line distributions ---
    if natural === :real
        if requested === :real
            return dist_from_mean_var(D, μ̄, σ̄²)
        else
            return _truncated_from_mean_var(D, μ̄, σ̄², lo, hi)
        end

    # --- Positive distributions [0,∞) ---
    elseif natural === :positive
        if requested === :half_right
            return _affine_shift(D, μ̄, σ̄², lo)
        elseif requested === :half_left
            return _affine_flip(D, μ̄, σ̄², hi)
        elseif requested === :bounded
            if lo >= 0
                return _truncated_from_mean_var(D, μ̄, σ̄², lo, hi)
            else
                throw(ArgumentError("Cannot place a $D (support [0,∞)) on [$lo, $hi] with lo < 0"))
            end
        elseif requested === :real
            throw(ArgumentError("Cannot place a $D (support [0,∞)) on (-∞,∞)"))
        end

    # --- Unit-interval distributions [0,1] ---
    elseif natural === :unit
        if requested === :bounded
            return _affine_scale(D, μ̄, σ̄², lo, hi)
        else
            throw(ArgumentError("Cannot place a $D (support [0,1]) on an unbounded interval"))
        end
    end

    throw(ArgumentError("Unsupported combination: $D on ($lo, $hi)"))
end

# --- Discrete support logic ---

function _dist_on_support_discrete(D, μ̄, σ̄², r::AbstractUnitRange)
    a = first(r)
    natural = _natural_support(D)

    if natural === :integer_bounded
        d = dist_from_mean_var(D, μ̄ - a, σ̄²)
        return a + d
    elseif natural === :integer_nonneg
        throw(ArgumentError("$D on a bounded range is not yet supported"))
    end

    throw(ArgumentError("$D does not have discrete support"))
end

# --- Affine transforms ---

function _build_on_standard(D, μ_std, σ²_std)
    # For DistSpec with 1 free param, use mean only (overdetermined with mean+var)
    if D isa DistSpec && length(free_params(D)) == 1
        return dist_from_mean(D, μ_std)
    else
        return dist_from_mean_var(D, μ_std, σ²_std)
    end
end

function _affine_shift(D, μ̄, σ̄², a)
    d = _build_on_standard(D, μ̄ - a, σ̄²)
    return Float64(a) + d
end

function _affine_flip(D, μ̄, σ̄², b)
    d = _build_on_standard(D, b - μ̄, σ̄²)
    return Float64(b) + (-1.0) * d
end

function _affine_scale(D, μ̄, σ̄², a, b)
    w = b - a
    μ_std = (μ̄ - a) / w
    σ²_std = σ̄² / w^2
    d = _build_on_standard(D, μ_std, σ²_std)
    return Float64(a) + Float64(w) * d
end

# --- Moment transform to standard support (for feasibility checking) ---

function _moments_to_standard(D, μ̄, σ̄², lo, hi)
    natural = _natural_support(D)
    requested = _requested_support_shape(lo, hi)

    if natural === :real
        requested === :real && return (μ̄, σ̄²)
        # Truncation feasibility: just check μ̄ is in bounds
        return (μ̄, σ̄²)
    elseif natural === :positive
        if requested === :half_right
            return (μ̄ - lo, σ̄²)
        elseif requested === :half_left
            return (hi - μ̄, σ̄²)
        elseif requested === :bounded && lo >= 0
            return (μ̄, σ̄²)
        end
    elseif natural === :unit
        if requested === :bounded
            w = hi - lo
            return ((μ̄ - lo) / w, σ̄² / w^2)
        end
    end

    throw(ArgumentError("Cannot check feasibility for $D on ($lo, $hi)"))
end

# --- Truncation ---

function _truncated_from_mean_var(D, μ̄, σ̄², lo, hi)
    return _solve_truncated_mean_var(D, Float64(lo), Float64(hi), Float64(μ̄), Float64(σ̄²))
end
