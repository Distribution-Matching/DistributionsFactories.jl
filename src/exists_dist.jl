# Feasibility predicates for moment-based distribution construction.
#
# Architecture:
#   - `_why_not_dist_from_mean_var(D, μ̄, σ̄²)` is the workhorse. It returns
#     `nothing` when `D` admits a distribution with mean `μ̄` and variance `σ̄²`,
#     otherwise a human-readable reason string. One method per supported dist.
#   - `exists_dist_from_mean_var(D, μ̄, σ̄²) :: Bool` is the public pure
#     predicate: `isnothing(_why_not_...)`.
#   - `_require_dist_from_mean_var(D, μ̄, σ̄²)` is the internal helper used by
#     constructors in moment_matching.jl: if infeasible, it throws
#     `DomainError(reason)`.
#
# Design rule: predicates never throw on well-formed numeric input. Only the
# constructors throw, and only with the reason string `_why_not_...` supplies.

"""
    exists_dist_from_mean_var(D, μ̄, σ̄²) -> Bool

Predicate: can a distribution of type `D` be constructed with mean `μ̄` and
variance `σ̄²`? Returns `false` for unsupported families and infeasible
`(μ̄, σ̄²)` pairs. Never throws on well-formed numeric input.

Use [`dist_from_mean_var`](@ref) to actually build the distribution; that
function throws `DomainError` with a reason when `exists_dist_from_mean_var`
returns `false`.
"""
exists_dist_from_mean_var(D, μ̄::Number, σ̄²::Number) =
    isnothing(_why_not_dist_from_mean_var(D, μ̄, σ̄²))

# Default: no method means not supported. Returning a reason (rather than
# silently `false` with no context) matters for the constructor path, which
# surfaces this message in DomainError.
function _why_not_dist_from_mean_var(D, μ̄::Number, σ̄²::Number)
    return "$D: distribution not supported"
end

# Used by constructors. If infeasible, raise DomainError carrying the reason.
function _require_dist_from_mean_var(D, μ̄::Number, σ̄²::Number)
    reason = _why_not_dist_from_mean_var(D, μ̄, σ̄²)
    isnothing(reason) || throw(DomainError((μ̄, σ̄²), reason))
    return nothing
end

# Shared across every continuous dist: σ̄² must be strictly positive.
function _why_not_positive_var(name, σ̄²)
    σ̄² > 0 || return "$name: the condition σ̄² > 0 is not satisfied"
    return nothing
end


# --- per-distribution feasibility ---

function _why_not_dist_from_mean_var(::Type{Beta}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("Beta", σ̄²); isnothing(r) || return r
    μ̄ > 0           || return "Beta: the condition μ̄ > 0 is not satisfied"
    μ̄ < 1           || return "Beta: the condition μ̄ < 1 is not satisfied"
    σ̄² < μ̄*(1-μ̄)   || return "Beta: the condition σ̄² < μ̄(1-μ̄) is not satisfied"
    return nothing
end

function _why_not_dist_from_mean_var(::Type{Uniform}, μ̄::Number, σ̄²::Number)
    return _why_not_positive_var("Uniform", σ̄²)
end

function _why_not_dist_from_mean_var(::Type{Normal}, μ̄::Number, σ̄²::Number)
    return _why_not_positive_var("Normal", σ̄²)
end

function _why_not_dist_from_mean_var(::Type{TDist}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("TDist", σ̄²); isnothing(r) || return r
    isapprox(μ̄, 0; atol=1e-12) || return "TDist: the condition μ̄ = 0 is not satisfied"
    σ̄² > 1                     || return "TDist: the condition σ̄² > 1 is not satisfied"
    return nothing
end

function _why_not_dist_from_mean_var(d::TDist, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("TDist", σ̄²); isnothing(r) || return r
    ν = dof(d)
    ν > 2 || return "TDist instance: the condition ν > 2 is not satisfied (ν=$ν). Variance is undefined for ν ≤ 2."
    return nothing
end

# Cauchy has no defined mean or variance: no (μ̄, σ̄²) input can produce one.
function _why_not_dist_from_mean_var(::Type{Cauchy}, μ̄::Number, σ̄²::Number)
    return "Cauchy: distribution has no finite mean or variance; construct via quantiles instead"
end

function _why_not_dist_from_mean_var(::Type{Logistic}, μ̄::Number, σ̄²::Number)
    return _why_not_positive_var("Logistic", σ̄²)
end

function _why_not_dist_from_mean_var(::Type{Laplace}, μ̄::Number, σ̄²::Number)
    return _why_not_positive_var("Laplace", σ̄²)
end

function _why_not_dist_from_mean_var(::Type{LogNormal}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("LogNormal", σ̄²); isnothing(r) || return r
    μ̄ > 0 || return "LogNormal: the condition μ̄ > 0 is not satisfied"
    return nothing
end

function _why_not_dist_from_mean_var(::Type{Chisq}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("Chisq", σ̄²); isnothing(r) || return r
    μ̄ > 0         || return "Chisq: the condition μ̄ > 0 is not satisfied"
    isinteger(μ̄)  || return "Chisq: the condition μ̄ ∈ ℕ is not satisfied"
    σ̄² == 2μ̄     || return "Chisq: the condition σ̄² = 2μ̄ is not satisfied"
    return nothing
end

function _why_not_dist_from_mean_var(::Type{Exponential}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("Exponential", σ̄²); isnothing(r) || return r
    μ̄ > 0 || return "Exponential: the condition μ̄ > 0 is not satisfied"
    isapprox(σ̄², μ̄^2; rtol=1e-10) || return "Exponential: the condition σ̄² = μ̄² is not satisfied"
    return nothing
end

function _why_not_dist_from_mean_var(::Type{Gamma}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("Gamma", σ̄²); isnothing(r) || return r
    μ̄ > 0 || return "Gamma: the condition μ̄ > 0 is not satisfied"
    return nothing
end

function _why_not_dist_from_mean_var(::Type{Erlang}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("Erlang", σ̄²); isnothing(r) || return r
    μ̄ > 0 || return "Erlang: the condition μ̄ > 0 is not satisfied"
    return nothing
end

function _why_not_dist_from_mean_var(::Type{Frechet}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("Frechet", σ̄²); isnothing(r) || return r
    μ̄ > 0 || return "Frechet: the condition μ̄ > 0 is not satisfied"
    return nothing
end

function _why_not_dist_from_mean_var(::Type{Weibull}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("Weibull", σ̄²); isnothing(r) || return r
    μ̄ > 0 || return "Weibull: the condition μ̄ > 0 is not satisfied"
    return nothing
end

function _why_not_dist_from_mean_var(::Type{Gumbel}, μ̄::Number, σ̄²::Number)
    return _why_not_positive_var("Gumbel", σ̄²)
end

function _why_not_dist_from_mean_var(::Type{Chi}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("Chi", σ̄²); isnothing(r) || return r
    isapprox(μ̄, √(2)*gamma((μ̄^2+σ̄²+1)/2)/gamma((μ̄^2+σ̄²)/2); rtol=1e-10, atol=1e-12) ||
        return "Chi: the condition μ̄ = √(2)Γ((μ̄²+σ̄²+1)/2)/Γ((μ̄²+σ̄²)/2) is not satisfied"
    return nothing
end

function _why_not_dist_from_mean_var(::Type{Rayleigh}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("Rayleigh", σ̄²); isnothing(r) || return r
    isapprox(√(σ̄²)/μ̄, √((4-π)/π); rtol=1e-10) ||
        return "Rayleigh: the condition CV = √((4-π)/π) is not satisfied"
    return nothing
end

function _why_not_dist_from_mean_var(::Type{FDist}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("FDist", σ̄²); isnothing(r) || return r
    μ̄ > 1 || return "FDist: the condition μ̄ > 1 is not satisfied"
    μ̄ < 2 || return "FDist: the condition μ̄ < 2 is not satisfied"
    σ̄² > μ̄^2*(μ̄-1)/(2-μ̄) || return "FDist: the condition σ̄² > μ̄²(μ̄-1)/(2-μ̄) is not satisfied"
    return nothing
end

function _why_not_dist_from_mean_var(::Type{InverseGamma}, μ̄::Number, σ̄²::Number)
    return _why_not_positive_var("InverseGamma", σ̄²)
end

function _why_not_dist_from_mean_var(::Type{Binomial}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("Binomial", σ̄²); isnothing(r) || return r
    σ̄² < μ̄ || return "Binomial: the condition μ̄ > σ̄² is not satisfied"
    μ̄ > 0   || return "Binomial: the condition μ̄ > 0 is not satisfied"
    n_raw = μ̄^2 / (μ̄ - σ̄²)
    isapprox(n_raw, round(n_raw); rtol=1e-8, atol=1e-8) ||
        return "Binomial: the condition μ̄²/(μ̄-σ̄²) ∈ ℕ is not satisfied (got n ≈ $n_raw)"
    return nothing
end

function _why_not_dist_from_mean_var(::Type{Poisson}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("Poisson", σ̄²); isnothing(r) || return r
    μ̄ == σ̄² || return "Poisson: the condition μ̄ = σ̄² is not satisfied"
    return nothing
end

function _why_not_dist_from_mean_var(::Type{NegativeBinomial}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("NegativeBinomial", σ̄²); isnothing(r) || return r
    μ̄ < σ̄² || return "NegativeBinomial: the condition μ̄ < σ̄² is not satisfied"
    return nothing
end

function _why_not_dist_from_mean_var(::Type{Pareto}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("Pareto", σ̄²); isnothing(r) || return r
    μ̄ > 0 || return "Pareto: the condition μ̄ > 0 is not satisfied"
    CV² = σ̄² / μ̄^2
    α = 1 + √(1 + 1/CV²)
    α > 2 || return "Pareto: the condition α > 2 is not satisfied (variance is infinite for α ≤ 2)"
    return nothing
end

function _why_not_dist_from_mean_var(::Type{FoldedNormal}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("FoldedNormal", σ̄²); isnothing(r) || return r
    μ̄ > 0 || return "FoldedNormal: the condition μ̄ > 0 is not satisfied"
    return nothing
end

function _why_not_dist_from_mean_var(::Type{Geometric}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("Geometric", σ̄²); isnothing(r) || return r
    μ̄ > 0 || return "Geometric: the condition μ̄ > 0 is not satisfied"
    isapprox(σ̄², μ̄ * (1 + μ̄); rtol=1e-10) ||
        return "Geometric: the condition σ̄² = μ̄(1+μ̄) is not satisfied"
    return nothing
end

# --- Truncated families sharing the Langevin envelope ---

function _why_not_truncexp_envelope(name::AbstractString, lo::Real, hi::Real,
                                    μ̄::Real, σ̄²::Real)
    (lo < μ̄ < hi) || return "Truncated $name: μ̄ must be in ($lo, $hi)"
    σ²_max = _truncexp_max_var(lo, hi, μ̄)
    # Non-strict with a relative tolerance: Laplace achieves the envelope
    # exactly when its mode falls outside [lo, hi] (the density reduces to a
    # truncated exponential there), and Normal/Logistic approach it
    # arbitrarily closely. The user-facing moments will be computed by
    # quadrature + truncation-normalization + inv-Langevin Newton; accumulated
    # relative error sits around 1e-10 in the worst cases. A 1e-8 tolerance
    # keeps clearly-infeasible points rejected without spuriously flagging
    # moments that came from a real truncated Normal/Laplace/Logistic.
    if σ̄² > σ²_max * (1 + 1e-8)
        return "Truncated $name: σ̄² = $σ̄² exceeds the Langevin feasibility " *
               "boundary σ²_max ≈ $(σ²_max) at μ̄ = $μ̄ on [$lo, $hi]. " *
               "The Normal/Laplace/Logistic families share this truncated-" *
               "exponential upper envelope; use a heavier-tailed family " *
               "(e.g. Student-t) to exceed it."
    end
    return nothing
end

# Half-truncated location-scale families with exponential-decay tails
# (Normal/Laplace/Logistic) on [lo, ∞) or (-∞, hi] reach an Exponential as the
# parent location is pushed to ∓∞. The right tail of an Exponential has
# CV = 1, so the feasibility region is σ̄² < (μ̄ - lo)² (resp. (hi - μ̄)²).
# See `exploration/feasible_half_truncated.jl` for the empirical confirmation
# and `report.tex` for the derivation.
function _why_not_half_trunc_exp_envelope(name::AbstractString,
                                          lo::Real, hi::Real,
                                          μ̄::Real, σ̄²::Real)
    if isfinite(lo) && !isfinite(hi)            # half-truncated below
        μ̄ > lo || return "Truncated $name: μ̄ must be > lo=$lo"
        gap = μ̄ - lo
        σ̄² < gap^2 || return "Truncated $name on [$lo, ∞): σ̄² must be < (μ̄ - lo)² = $(gap^2) (exponential-tail bound)"
    elseif !isfinite(lo) && isfinite(hi)        # half-truncated above
        μ̄ < hi || return "Truncated $name: μ̄ must be < hi=$hi"
        gap = hi - μ̄
        σ̄² < gap^2 || return "Truncated $name on (-∞, $hi]: σ̄² must be < (hi - μ̄)² = $(gap^2) (exponential-tail bound)"
    end
    return nothing
end

function _why_not_dist_from_mean_var(d::Truncated{<:Normal}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("Normal", σ̄²); isnothing(r) || return r
    lo, hi = extrema(d)
    return isfinite(lo) && isfinite(hi) ?
        _why_not_truncexp_envelope("Normal", lo, hi, μ̄, σ̄²) :
        _why_not_half_trunc_exp_envelope("Normal", lo, hi, μ̄, σ̄²)
end

function _why_not_dist_from_mean_var(d::Truncated{<:Laplace}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("Laplace", σ̄²); isnothing(r) || return r
    lo, hi = extrema(d)
    if isfinite(lo) && isfinite(hi)
        return _why_not_truncexp_envelope("Laplace", lo, hi, μ̄, σ̄²)
    end
    # Half-truncated Laplace: the boundary σ̄² = gap² is *attained* (not just
    # approached) by parents with μp ≤ lo (resp. ≥ hi), where the truncated
    # distribution is exactly an Exponential. Allow equality with a small
    # numerical tolerance.
    if isfinite(lo) && !isfinite(hi)
        μ̄ > lo || return "Truncated Laplace: μ̄ must be > lo=$lo"
        gap = μ̄ - lo
        σ̄² ≤ gap^2 * (1 + 1e-10) || return "Truncated Laplace on [$lo, ∞): σ̄² must be ≤ (μ̄ - lo)² = $(gap^2) (exponential bound, attained)"
    elseif !isfinite(lo) && isfinite(hi)
        μ̄ < hi || return "Truncated Laplace: μ̄ must be < hi=$hi"
        gap = hi - μ̄
        σ̄² ≤ gap^2 * (1 + 1e-10) || return "Truncated Laplace on (-∞, $hi]: σ̄² must be ≤ (hi - μ̄)² = $(gap^2) (exponential bound, attained)"
    end
    return nothing
end

function _why_not_dist_from_mean_var(d::Truncated{<:Logistic}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("Logistic", σ̄²); isnothing(r) || return r
    lo, hi = extrema(d)
    return isfinite(lo) && isfinite(hi) ?
        _why_not_truncexp_envelope("Logistic", lo, hi, μ̄, σ̄²) :
        _why_not_half_trunc_exp_envelope("Logistic", lo, hi, μ̄, σ̄²)
end

# Half-truncated Student-t on [lo, ∞) or (-∞, hi]. The right tail of t_ν is
# Pareto-like with tail index ν; pushing the parent location to ∓∞ leaves a
# Pareto whose CV is √(ν/(ν-2)) for ν > 2. So the feasibility region is
# σ̄² < ν/(ν-2) · (μ̄ - lo)² (resp. (hi - μ̄)²); infeasible for ν ≤ 2.
# See `exploration/feasible_half_truncated_t_ceiling.{jl,pdf}`.
function _why_not_dist_from_mean_var(d::Truncated{<:TDist}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("TDist", σ̄²); isnothing(r) || return r
    ν = dof(d.untruncated)
    ν > 2 || return "Truncated TDist: requires ν > 2 for variance to exist (got ν=$ν)"
    lo, hi = extrema(d)
    if isfinite(lo) && isfinite(hi)
        # Two-sided: defer to the doubly-truncated polynomial-tail dome
        # (not yet implemented as a feasibility predicate). Fall back to the
        # universal Bhatia–Davis ceiling, which is loose but correct.
        σ̄² < (hi - μ̄) * (μ̄ - lo) || return "Truncated TDist on [$lo, $hi]: σ̄² must be < (hi-μ̄)(μ̄-lo) (Bhatia–Davis universal bound; tight Pareto-tail dome not yet implemented)"
        return nothing
    end
    cap = ν / (ν - 2)
    if isfinite(lo) && !isfinite(hi)
        μ̄ > lo || return "Truncated TDist: μ̄ must be > lo=$lo"
        gap = μ̄ - lo
        σ̄² < cap * gap^2 || return "Truncated TDist on [$lo, ∞) with ν=$ν: σ̄² must be < ν/(ν-2)·(μ̄-lo)² = $(cap*gap^2) (Pareto-tail bound)"
    elseif !isfinite(lo) && isfinite(hi)
        μ̄ < hi || return "Truncated TDist: μ̄ must be < hi=$hi"
        gap = hi - μ̄
        σ̄² < cap * gap^2 || return "Truncated TDist on (-∞, $hi] with ν=$ν: σ̄² must be < ν/(ν-2)·(hi-μ̄)² = $(cap*gap^2) (Pareto-tail bound)"
    end
    return nothing
end

# Truncated Poisson has only λ as a free parameter. Mean+var is overdetermined;
# only the mean constraint is enforced here (variance becomes a derived check
# inside the constructor). Feasibility: μ̄ in the open truncation interval.
function _why_not_dist_from_mean_var(d::Truncated{<:Poisson}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("Poisson", σ̄²); isnothing(r) || return r
    lo, hi = d.lower, d.upper
    (lo < μ̄ < hi) || return "Truncated Poisson: μ̄ must be in ($lo, $hi)"
    return nothing
end

# --- Distributions whose feasibility rules aren't encoded yet ---

function _why_not_dist_from_mean_var(::Type{TriangularDist}, μ̄::Number, σ̄²::Number)
    return "TriangularDist: feasibility rule not yet implemented"
end

function _why_not_dist_from_mean_var(::Type{SymTriangularDist}, μ̄::Number, σ̄²::Number)
    return _why_not_positive_var("SymTriangularDist", σ̄²)
end

# DiscreteSymmetricTriangular: var = n(n+2)/6 with n ∈ ℕ₀. Solving:
#   n = -1 + √(1 + 6σ̄²). Feasible iff μ̄ ∈ ℤ and that n is a non-negative integer.
function _why_not_dist_from_mean_var(::Type{DiscreteSymmetricTriangular}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("DiscreteSymmetricTriangular", σ̄²); isnothing(r) || return r
    isapprox(μ̄, round(μ̄); atol=1e-8) || return "DiscreteSymmetricTriangular: μ̄ must be an integer (got $μ̄)"
    σ̄² ≥ 0 || return "DiscreteSymmetricTriangular: σ̄² must be ≥ 0"
    n_raw = -1 + √(1 + 6σ̄²)
    (isapprox(n_raw, round(n_raw); atol=1e-8) && round(n_raw) ≥ 0) ||
        return "DiscreteSymmetricTriangular: half-width n = -1 + √(1+6σ̄²) must be a non-negative integer (got n ≈ $n_raw)"
    return nothing
end

# DiscreteTriangular has 3 integer parameters and 2 moment constraints; mean+var
# alone is underdetermined. Surface that to the user (the mean+var+mode factory
# is the right entry point).
function _why_not_dist_from_mean_var(::Type{DiscreteTriangular}, μ̄::Number, σ̄²::Number)
    return "DiscreteTriangular: mean+var alone is underdetermined (3 integer params); supply `mode` as well"
end

function _why_not_dist_from_mean_var(::Type{DiscreteUniform}, μ̄::Number, σ̄²::Number)
    r = _why_not_positive_var("DiscreteUniform", σ̄²); isnothing(r) || return r
    n_raw = -1 + √(1 + 12 * σ̄²)
    (isapprox(n_raw, round(n_raw); atol=1e-8) && round(n_raw) ≥ 0) ||
        return "DiscreteUniform: n = b - a must be a non-negative integer (got n ≈ $n_raw)"
    n = round(Int, n_raw)
    a_raw = μ̄ - n / 2
    isapprox(a_raw, round(a_raw); atol=1e-8) ||
        return "DiscreteUniform: lower bound a must be an integer (got a ≈ $a_raw)"
    return nothing
end
