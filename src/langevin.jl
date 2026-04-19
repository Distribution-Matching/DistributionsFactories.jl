# Langevin-function feasibility boundary for truncated Normal / Laplace / Logistic.
#
# Background (see exploration/report.tex in the companion exploration repo):
# for any continuous, unimodal, location–scale family with exponential tails
# (Normal, Laplace, Logistic) truncated to a bounded interval [a, b], the set of
# achievable (mean, variance) pairs forms a dome whose upper envelope coincides
# with the truncated-exponential family
#
#     p_λ(x) = e^{-λ x} / Z(λ),   x ∈ [a, b],   λ ∈ ℝ.
#
# Centering by c = (a+b)/2 and writing w = (b-a)/2, its moments satisfy
#
#     mean(λ)    = c − w · L(λ w),
#     var(λ)     = w² · L'(λ w),
#
# where L(x) = coth(x) − 1/x is the Langevin function. Inverting, the maximum
# variance achievable at mean μ̄ is
#
#     σ²_max(μ̄) = w² · L'( L^{-1}( (c − μ̄) / w ) ).
#
# A truncated-Normal / Laplace / Logistic with mean μ̄ and variance σ̄² exists
# iff a < μ̄ < b and σ̄² < σ²_max(μ̄).

"""
    langevin(x)

Langevin function `L(x) = coth(x) − 1/x`.

For `|x| < 1e-4` we use the Maclaurin expansion
`L(x) = x/3 − x³/45 + 2x⁵/945 − x⁷/4725 + O(x⁹)`
(derived from the well-known series of `coth(x) = 1/x + x/3 − x³/45 + …`).
The expansion avoids the `coth(x) − 1/x` cancellation at small `x`, where
both terms individually blow up like `1/x` but their difference is `O(x)`.
The cutoff `1e-4` keeps the next neglected term `O(x⁹/93555)` below `1e-36`,
well inside double precision.
"""
function langevin(x::Real)
    if abs(x) < 1e-4
        return x/3 - x^3/45 + 2x^5/945 - x^7/4725
    else
        return coth(x) - 1/x
    end
end

"""
    langevin_deriv(x)

Derivative `L'(x) = 1/x² − 1/sinh²(x)`.

For `|x| < 1e-3` we use
`L'(x) = 1/3 − x²/15 + 2x⁴/189 − x⁶/675 + O(x⁸)`
(term-by-term derivative of the series above). As with `L`, the closed form
suffers from catastrophic cancellation near zero — `1/x²` and `1/sinh²(x)`
individually behave like `1/x²` — while the true value tends to `1/3`. The
cutoff is slightly larger than in `langevin` because the series here is
dominated by constant and `O(x²)` terms, so the cancellation bites sooner.
"""
function langevin_deriv(x::Real)
    if abs(x) < 1e-3
        return 1/3 - x^2/15 + 2x^4/189 - x^6/675
    else
        s = sinh(x)
        return 1/x^2 - 1/(s*s)
    end
end

"""
    inv_langevin(y)

Inverse Langevin function `L^{-1}(y)` on `y ∈ (-1, 1)`.

The initial guess is Cohen's (1991) [3/2] Padé approximant
`y(3 − y²)/(1 − y²)`, which is exact to `O(y³)` near the origin (matching
`L^{-1}(y) = 3y + (9/5) y³ + …`) and captures the `±1` pole structure at the
domain edges, giving fewer-than-five-iterations Newton convergence across the
whole interval. Refinement is plain Newton on `L(z) − y = 0` using
`L'` above; we cap at 50 iterations as a safety net but exit as soon as
`|L(z) − y| < 1e-14`.
"""
function inv_langevin(y::Real)
    abs(y) < 1 || throw(DomainError(y, "inv_langevin: |y| must be < 1"))
    y == 0 && return zero(float(y))
    z = y * (3 - y^2) / (1 - y^2)
    for _ in 1:50
        f = langevin(z) - y
        abs(f) < 1e-14 && return z
        z -= f / langevin_deriv(z)
    end
    return z
end

"""
    _truncexp_max_var(a, b, μ̄)

Upper-envelope variance at mean `μ̄` for the truncated-exponential family on
`[a, b]`. This is the shared feasibility boundary for truncated Normal,
Laplace, and Logistic: a distribution in any of these families with mean `μ̄`
and variance `σ̄²` exists iff `σ̄² < _truncexp_max_var(a, b, μ̄)`.
"""
function _truncexp_max_var(a::Real, b::Real, μ̄::Real)
    a < b || throw(ArgumentError("Require a < b, got a=$a, b=$b"))
    a < μ̄ < b || throw(DomainError(μ̄, "μ̄ must lie strictly in ($a, $b)"))
    c = (a + b) / 2
    w = (b - a) / 2
    z = inv_langevin((c - μ̄) / w)
    return w^2 * langevin_deriv(z)
end
