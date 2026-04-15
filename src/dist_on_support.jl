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
_natural_support(::Type{FoldedNormal}) = :positive
_natural_support(::Type{<:Union{Beta,Uniform}}) = :unit

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

# --- Entry point ---

"""
    dist_from_mean_var(D, μ̄, σ̄², support)

Construct a distribution of type `D` with mean `μ̄` and variance `σ̄²` on the given
`support` interval. Accepts IntervalSets intervals (`2..7`, `3..Inf`) or
Distributions.jl `RealInterval`.

The function automatically determines whether to use an **affine transform** or
**truncation** based on the relationship between `D`'s natural support and the
requested support:

- **Affine**: when the requested support has the same shape as the natural one
  (e.g. Beta on `[2,7]`, Gamma on `[3,∞)`). Returns a `LocationScale` wrapper.
- **Truncation**: when the requested support is strictly contained in the natural
  one (e.g. Normal on `[0,1]`, Gamma on `[0,10]`). Returns a `Truncated` wrapper.

# Examples
```julia
dist_from_mean_var(Beta, 3.5, 0.5, 2..7)        # affine: Beta scaled to [2,7]
dist_from_mean_var(Gamma, 5.0, 3.0, 2..Inf)     # affine: Gamma shifted to [2,∞)
dist_from_mean_var(Gamma, 5.0, 3.0, -Inf..8)    # affine: flipped Gamma on (-∞,8]
dist_from_mean_var(Normal, 0.5, 0.04, 0..1)     # truncation: Normal on [0,1]
dist_from_mean_var(Gamma, 3.0, 1.0, 0..10)      # truncation: Gamma on [0,10]
```
"""
function dist_from_mean_var(D::Type{<:Distribution}, μ̄::Number, σ̄²::Number, support::AbstractInterval)
    lo, hi = endpoints(support)
    return _dist_on_support(D, μ̄, σ̄², lo, hi)
end

function dist_from_mean_var(D::Type{<:Distribution}, μ̄::Number, σ̄²::Number, support::Distributions.RealInterval)
    return _dist_on_support(D, μ̄, σ̄², support.lb, support.ub)
end

function _dist_on_support(D, μ̄, σ̄², lo, hi)
    natural = _natural_support(D)
    requested = _requested_support_shape(lo, hi)

    # --- Real-line distributions ---
    if natural === :real
        if requested === :real
            return dist_from_mean_var(D, μ̄, σ̄²)
        else
            # Any restriction of (-∞,∞) is truncation
            return _truncated_from_mean_var(D, μ̄, σ̄², lo, hi)
        end

    # --- Positive distributions [0,∞) ---
    elseif natural === :positive
        if requested === :half_right
            # [a, ∞): shift by a
            return _affine_shift(D, μ̄, σ̄², lo)
        elseif requested === :half_left
            # (-∞, b]: flip and shift
            return _affine_flip(D, μ̄, σ̄², hi)
        elseif requested === :bounded
            # [a, b]: truncation (possibly after shift)
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
            # [a, b]: affine transform from [0,1] to [a,b]
            return _affine_scale(D, μ̄, σ̄², lo, hi)
        else
            throw(ArgumentError("Cannot place a $D (support [0,1]) on an unbounded interval"))
        end
    end

    throw(ArgumentError("Unsupported combination: $D on ($lo, $hi)"))
end

# --- Affine transforms ---

function _affine_shift(D, μ̄, σ̄², a)
    # [0,∞) → [a,∞): Y = a + X, so X has mean μ̄-a, var σ̄²
    d = dist_from_mean_var(D, μ̄ - a, σ̄²)
    return LocationScale(Float64(a), 1.0, d)
end

function _affine_flip(D, μ̄, σ̄², b)
    # [0,∞) → (-∞,b]: Y = b - X, so X has mean b-μ̄, var σ̄²
    d = dist_from_mean_var(D, b - μ̄, σ̄²)
    return LocationScale(Float64(b), -1.0, d; check_args=false)
end

function _affine_scale(D, μ̄, σ̄², a, b)
    # [0,1] → [a,b]: Y = a + (b-a)X, so X has mean (μ̄-a)/(b-a), var σ̄²/(b-a)²
    w = b - a
    μ_std = (μ̄ - a) / w
    σ²_std = σ̄² / w^2
    d = dist_from_mean_var(D, μ_std, σ²_std)
    return LocationScale(Float64(a), Float64(w), d)
end

# --- Truncation ---

function _truncated_from_mean_var(D, μ̄, σ̄², lo, hi)
    return _solve_truncated_mean_var(D, Float64(lo), Float64(hi), Float64(μ̄), Float64(σ̄²))
end
