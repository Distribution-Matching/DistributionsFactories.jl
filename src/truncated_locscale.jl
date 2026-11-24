# Fast moment matching for the truncated Normal and the truncated
# location-scale Student-t.
#
#   truncated_locscale_moments(family, μ, σ, a, b) -> (mean, var)
#   truncated_locscale_params(family, a, b, m, v)  -> (μ, σ)
#
# `family` is the type `Normal` or an instance `TDist(ν)` (ν > 0, possibly
# non-integer). Bounds may each be finite or infinite.
#
# Moments. Write α = (a - μ)/σ, β = (b - μ)/σ and let Y be the standard law
# (N(0,1) or t_ν) truncated to [α, β], with mass Z = F(β) - F(α).
#
#   Normal:  E[Y]  = (φ(α) - φ(β)) / Z,
#            E[Y²] = 1 + (α φ(α) - β φ(β)) / Z.
#
#   t_ν:     with g(t) = (ν + t²) f_ν(t) one has g'(t) = -(ν - 1) t f_ν(t), so
#            E[Y]  = (g(α) - g(β)) / ((ν - 1) Z),
#            and integrating t · (t f_ν) by parts with the same identity,
#            E[Y²] = (ν - [t g(t)]_α^β / Z) / (ν - 2).
#            This avoids the t_{ν-2} CDF and has no cancellation as ν → ∞,
#            where it tends to the Normal formulas. On a bounded interval both
#            identities hold for every ν ∉ {1, 2}; with an infinite end they
#            need ν > 1 (mean) and ν > 2 (second moment).
#
# Stability. Moments are returned relative to the mode r = clamp(0, α, β) of
# the truncated standard law, as δ = E[Y] - r, so that the mean in user
# coordinates is clamp(μ, a, b) + σ δ, without cancellation even when |μ| is
# huge. The closed forms lose accuracy when the spread of Y is small relative
# to |r|: far in a tail, or on an interval that is narrow relative to σ (for
# instance a truncated normal close to the boundary curve of its feasible
# region, where σ ≫ b - a). There the moments are computed instead by
# 64-point Gauss–Legendre quadrature of exp(ℓ(r + s) - ℓ(r)) on the window
# where the log-density is within 60 of its maximum. The integrand is then
# smooth on the window, the quadrature is accurate to about 1e-15, and the
# log-density difference is evaluated without cancellation.
#
# Everything is generic over `Real`, so the moments can be differentiated
# with ForwardDiff. The CDFs are wrapped so that `Dual` arguments get their
# derivative from the density.

# ---------------------------------------------------------------- exceptions

"""
    InfeasibleMomentsError(family, lower, upper, mean, var, proven, msg)

Thrown by [`truncated_locscale_params`](@ref) and by `make_dist` /
`dist_from_mean_var` for truncated Normal, Laplace, Logistic and Student-t
templates when no member of the family, truncated to `[lower, upper]`, has
the requested `mean` and `var`.

`proven` is `true` when the target violates a known feasibility condition,
and `false` when the condition is not known (two-sided Student-t) and the
solver failed to converge; in the latter case the target is presumably,
but not provably, infeasible.

Catch it to fall back to another construction:
```julia
try
    μ, σ = truncated_locscale_params(TDist(4), 0.0, 1.0, m, v)
catch e
    e isa InfeasibleMomentsError || rethrow()
    # fall back
end
```
"""
struct InfeasibleMomentsError <: Exception
    family::String
    lower::Float64
    upper::Float64
    mean::Float64
    var::Float64
    proven::Bool
    msg::String
end

function Base.showerror(io::IO, e::InfeasibleMomentsError)
    print(io, "InfeasibleMomentsError: ", e.family, " truncated to [",
          e.lower, ", ", e.upper, "] with mean ", e.mean, " and variance ",
          e.var, ": ", e.msg)
end

_infeasible(fam, a, b, m, v, msg; proven::Bool = true) =
    throw(InfeasibleMomentsError(fam, a, b, m, v, proven, msg))

_famname(::Type{Normal}) = "Normal"
_famname(d::TDist) = "TDist($(dof(d)))"

# ---------------------------------------------------------------- quadrature

# Gauss–Legendre nodes and weights on [-1, 1] by Newton's method on P_n.
function _gauss_legendre(n::Int)
    x = zeros(n); w = zeros(n)
    for i in 1:n
        z = cos(π * (i - 0.25) / (n + 0.5))
        pp = 0.0
        for _ in 1:100
            p1, p2 = 1.0, 0.0
            for j in 1:n
                p3 = p2; p2 = p1
                p1 = ((2j - 1) * z * p2 - (j - 1) * p3) / j
            end
            pp = n * (z * p1 - p2) / (z^2 - 1)
            z1 = z
            z = z1 - p1 / pp
            if abs(z - z1) < 1e-15
                p1, p2 = 1.0, 0.0
                for j in 1:n
                    p3 = p2; p2 = p1
                    p1 = ((2j - 1) * z * p2 - (j - 1) * p3) / j
                end
                pp = n * (z * p1 - p2) / (z^2 - 1)
                break
            end
        end
        x[i] = z
        w[i] = 2 / ((1 - z^2) * pp^2)
    end
    return Tuple(x), Tuple(w)
end

const _GL_NODES, _GL_WEIGHTS = _gauss_legendre(64)
const _GL_DROP = 60.0   # window: log-density within 60 of its maximum

# Mean offset and variance of the density ∝ exp(Δℓ(s)) on [slo, shi], by
# Gauss–Legendre quadrature (two passes, so the variance is centred).
function _gl_moments(Δℓ, slo, shi)
    T = promote_type(typeof(slo), typeof(shi))
    h = (shi - slo) / 2
    c = (shi + slo) / 2
    Z = zero(T); S1 = zero(T)
    @inbounds for i in eachindex(_GL_NODES)
        s = c + h * _GL_NODES[i]
        e = _GL_WEIGHTS[i] * exp(Δℓ(s))
        Z += e; S1 += e * s
    end
    δ = S1 / Z
    S2 = zero(T)
    @inbounds for i in eachindex(_GL_NODES)
        s = c + h * _GL_NODES[i]
        S2 += _GL_WEIGHTS[i] * exp(Δℓ(s)) * (s - δ)^2
    end
    return δ, S2 / Z
end

# Integral of f over [lo, hi] (lo ≤ 0 ≤ hi, both finite) by composite
# Gauss–Legendre on geometric panels with edges 0, ±1, ±2, ±4, …, so that
# slowly decaying (heavy-tailed) integrands on very wide intervals are
# resolved. Generic over Real (no adaptive step selection on Duals).
function _gl_panels(f, lo, hi)
    acc = zero(f(zero(promote_type(typeof(lo), typeof(hi)))))
    for (sgn, L) in ((1, hi), (-1, -lo))
        x0 = 0.0
        while x0 < L
            x1f = max(1.0, 2x0)
            last = x1f ≥ L
            x1 = last ? L : x1f
            h = (x1 - x0) / 2; c = (x1 + x0) / 2
            @inbounds for i in eachindex(_GL_NODES)
                acc += h * _GL_WEIGHTS[i] * f(sgn * (c + h * _GL_NODES[i]))
            end
            last && break
            x0 = x1f
        end
    end
    return acc
end

function _gl_panels_moments(Δℓ, lo, hi)
    f(s) = exp(Δℓ(s))
    Z = _gl_panels(f, lo, hi)
    δ = _gl_panels(s -> s * f(s), lo, hi) / Z
    return δ, _gl_panels(s -> (s - δ)^2 * f(s), lo, hi) / Z
end

# ---------------------------------------------------------------- standard laws

# Densities and (complementary) CDFs. The CDFs are evaluated in Float64 and,
# for ForwardDiff.Dual arguments, differentiated through the density, since
# `beta_inc` (used by the t CDF) is not generic.
_npdf(x::Real) = exp(-x^2 / 2) / sqrt(2π)
_ncdf(x::Real) = cdf(Normal(), x)
_nccdf(x::Real) = ccdf(Normal(), x)

_t_logc(ν) = loggamma((ν + 1) / 2) - loggamma(ν / 2) - log(ν * π) / 2
_tpdf(ν, x::Real) = exp(_t_logc(ν) - (ν + 1) / 2 * log1p(x^2 / ν))
# t tail P(T > x) for x ≥ 0 via the regularised incomplete beta function,
# with the complementary argument passed explicitly. Distributions' `ccdf`
# for TDist loses accuracy far in the tail (relative error ~1e-8 at x = 3e4,
# ~1e-5 at x = 1e6 for ν = 3); this form stays at ~1e-15.
function _ttail(ν, x::Real)
    isinf(x) && return zero(float(x))
    x2 = x^2
    return beta_inc(ν / 2, oftype(float(x), 0.5), ν / (ν + x2), x2 / (ν + x2))[1] / 2
end
_tcdf(ν, x::Real) = x ≤ 0 ? _ttail(ν, -x) : 1 - _ttail(ν, x)
_tccdf(ν, x::Real) = x ≥ 0 ? _ttail(ν, x) : 1 - _ttail(ν, -x)

function _ncdf(x::ForwardDiff.Dual{T}) where {T}
    v = ForwardDiff.value(x)
    return ForwardDiff.Dual{T}(_ncdf(v), _npdf(v) * ForwardDiff.partials(x))
end
function _nccdf(x::ForwardDiff.Dual{T}) where {T}
    v = ForwardDiff.value(x)
    return ForwardDiff.Dual{T}(_nccdf(v), -_npdf(v) * ForwardDiff.partials(x))
end
function _tcdf(ν, x::ForwardDiff.Dual{T}) where {T}
    v = ForwardDiff.value(x)
    return ForwardDiff.Dual{T}(_tcdf(ν, v), _tpdf(ν, v) * ForwardDiff.partials(x))
end
function _tccdf(ν, x::ForwardDiff.Dual{T}) where {T}
    v = ForwardDiff.value(x)
    return ForwardDiff.Dual{T}(_tccdf(ν, v), -_tpdf(ν, v) * ForwardDiff.partials(x))
end

# Mass of [α, β], using the tail on the side where the interval lies.
_nmass(α, β) = α > 0 ? _nccdf(α) - _nccdf(β) : _ncdf(β) - _ncdf(α)
_tmass(ν, α, β) = α > 0 ? _tccdf(ν, α) - _tccdf(ν, β) : _tcdf(ν, β) - _tcdf(ν, α)

# Internal family keys: `_NormalKey()` for the Normal and ν for t_ν. (A
# singleton instead of the type `Normal`, so that Julia specialises on it.)
struct _NormalKey end

# Closed forms: (E[Y], E[Y²]) of the standard law truncated to [α, β].
function _closed_moments(::_NormalKey, α, β)
    Z = _nmass(α, β)
    φa = isfinite(α) ? _npdf(α) : zero(Z)
    φb = isfinite(β) ? _npdf(β) : zero(Z)
    aφa = isfinite(α) ? α * φa : zero(Z)
    bφb = isfinite(β) ? β * φb : zero(Z)
    return (φa - φb) / Z, 1 + (aφa - bφb) / Z
end

function _closed_moments(ν::Real, α, β)
    Z = _tmass(ν, α, β)
    ga = isfinite(α) ? (ν + α^2) * _tpdf(ν, α) : zero(Z)
    gb = isfinite(β) ? (ν + β^2) * _tpdf(ν, β) : zero(Z)
    tga = isfinite(α) ? α * ga : zero(Z)
    tgb = isfinite(β) ? β * gb : zero(Z)
    return (ga - gb) / ((ν - 1) * Z), (ν - (tgb - tga) / Z) / (ν - 2)
end

# Log-density difference Δℓ(s) = ℓ(r + s) - ℓ(r), without cancellation.
_Δℓ(::_NormalKey, r) = s -> -s * (r + s / 2)
_Δℓ(ν::Real, r) = s -> -(ν + 1) / 2 * log1p(s * (2r + s) / (ν + r^2))

# Half-width of the window around 0 outside of which the log-density is more
# than _GL_DROP below its value at r.
_window(::_NormalKey, r) = sqrt(r^2 + 2_GL_DROP)
_window(ν::Real, r) = sqrt(ν * ((1 + r^2 / ν) * exp(2_GL_DROP / (ν + 1)) - 1))

# ν at a pole of the t closed forms (only relevant on a bounded interval:
# with an infinite end, ν > 2 and the pole at ν = 2 is a genuine divergence
# of the variance, which the closed form captures).
_near_pole(::_NormalKey) = false
_near_pole(ν::Real) = abs(ν - 1) < 1e-4 || abs(ν - 2) < 1e-4

# For large ν the incomplete-beta t CDF loses accuracy (about 1e-9 at
# ν = 10⁴), while the law is close to normal with a short high-density
# window, so quadrature is used throughout.
_prefer_quadrature(::_NormalKey) = false
_prefer_quadrature(ν::Real) = ν ≥ 200

# _std_moments(fam, α, β, slo, shi) -> (δ, v)
#
# For the standard law of `fam` (`_NormalKey()` or ν for t_ν) truncated to
# [α, β], return δ = E[Y] - r and v = Var(Y), where r = clamp(0, α, β) is the
# mode. slo = α - r and shi = β - r are passed separately, computed in user
# coordinates by the caller, because β - α loses accuracy when |α| is large.
# Closed forms where they are well conditioned, Gauss–Legendre quadrature on
# the high-density window otherwise.
function _std_moments(fam, α, β, slo, shi)
    r = clamp(zero(α), α, β)
    Δℓ = _Δℓ(fam, r)
    R = _window(fam, r)
    lo, hi = max(slo, -R - r), min(shi, R - r)
    bounded = isfinite(slo) && isfinite(shi)
    # A bounded interval that is narrow on the scale of the law (width ≤ 1),
    # or narrow relative to its distance from 0: the closed forms cancel and
    # the integrand is smooth on the window.
    narrow = bounded &&
             (shi - slo ≤ 1 || (α * β > 0 && shi - slo ≤ 4 * min(abs(α), abs(β))))
    if narrow || _prefer_quadrature(fam)
        return _gl_moments(Δℓ, lo, hi)
    elseif bounded && _near_pole(fam)
        # Heavy tail (ν ≈ 1 or 2): the variance is dominated by the far end
        # of the interval, so integrate over all of it, not the window.
        return _gl_panels_moments(Δℓ, slo, shi)
    end
    m, E2 = _closed_moments(fam, α, β)
    v = E2 - m^2
    # Cancellation in E2 - m² costs about eps·max(1, r²)/v in relative terms.
    v > 1e-3 * max(one(v), r^2) && return m - r, v
    return _gl_moments(Δℓ, lo, hi)
end

# ---------------------------------------------------------------- public moments

_famkey(::Type{Normal}) = _NormalKey()
_famkey(d::TDist) = dof(d)

"""
    truncated_locscale_moments(family, μ, σ, a, b) -> (mean, var)

Mean and variance of the location-scale law `μ + σ X`, with `X` standard
Normal (`family = Normal`) or Student-t (`family = TDist(ν)`), truncated to
`[a, b]`. Either bound may be infinite.

Closed forms (see the source for the derivation), with a Gauss–Legendre
fallback in regimes where they cancel; the results are accurate to about
1e-13 relative in all regimes, and generic over `Real`, so they can be
differentiated with ForwardDiff. For the t, a finite mean needs `ν > 1` and
a finite variance `ν > 2` when a bound is infinite; on a bounded interval
every `ν > 0` is allowed.
"""
truncated_locscale_moments(::Type{Normal}, μ::Real, σ::Real, a::Real, b::Real) =
    _tls_moments(_NormalKey(), μ, σ, a, b)
truncated_locscale_moments(d::TDist, μ::Real, σ::Real, a::Real, b::Real) =
    _tls_moments(dof(d), μ, σ, a, b)

function _tls_moments(key, μ::Real, σ::Real, a::Real, b::Real)
    a < b || throw(ArgumentError("need a < b, got a = $a, b = $b"))
    σ > 0 || throw(DomainError(σ, "σ must be positive"))
    if key isa Real && !(isfinite(a) && isfinite(b))
        key > 2 || throw(DomainError(key,
            "TDist(ν) truncated to an unbounded interval needs ν > 2 for a finite variance"))
    end
    T = promote_type(typeof(μ), typeof(σ), Float64)
    α = isfinite(a) ? (a - μ) / σ : T(a)
    β = isfinite(b) ? (b - μ) / σ : T(b)
    mode = clamp(μ, a, b)
    slo = isfinite(a) ? (a - mode) / σ : T(a)
    shi = isfinite(b) ? (b - mode) / σ : T(b)
    δ, v = _std_moments(key, α, β, slo, shi)
    return mode + σ * δ, σ^2 * v
end

# ---------------------------------------------------------------- solvers

# Damped Newton on the 2×2 system F(x) = 0, with the Jacobian from ForwardDiff
# dual numbers. `valid(x1, x2)` rejects steps leaving the parameter domain.
# Returns (x1, x2, converged).
struct _NewtonTag end

function _newton2(resid, x1::Float64, x2::Float64, valid;
                  tol::Float64, maxiter::Int)
    for _ in 1:maxiter
        d1 = ForwardDiff.Dual{_NewtonTag}(x1, 1.0, 0.0)
        d2 = ForwardDiff.Dual{_NewtonTag}(x2, 0.0, 1.0)
        F1, F2 = resid(d1, d2)
        f1, f2 = ForwardDiff.value(F1), ForwardDiff.value(F2)
        isfinite(f1) && isfinite(f2) || return x1, x2, false
        nF = max(abs(f1), abs(f2))
        nF < tol && return x1, x2, true
        J11, J12 = ForwardDiff.partials(F1)
        J21, J22 = ForwardDiff.partials(F2)
        det = J11 * J22 - J12 * J21
        (isfinite(det) && det != 0) || return x1, x2, false
        s1 = -(J22 * f1 - J12 * f2) / det
        s2 = -(J11 * f2 - J21 * f1) / det
        t = 1.0
        accepted = false
        for _ in 1:60
            y1, y2 = x1 + t * s1, x2 + t * s2
            if valid(y1, y2)
                g1, g2 = resid(y1, y2)
                if isfinite(g1) && isfinite(g2) &&
                   max(abs(g1), abs(g2)) < (1 - 1e-4 * t) * nF
                    accepted = true
                    break
                end
            end
            t /= 2
        end
        # No decrease along the Newton direction: accept if already at the
        # attainable accuracy (the residual floor near the boundary curve).
        accepted || return x1, x2, nF < 1e3 * tol
        x1 += t * s1
        x2 += t * s2
    end
    return x1, x2, false
end

# Normal on a canonical domain, in natural parameters θ = (μ/σ², -1/(2σ²)).
# Near the boundary curve μ ≈ -λσ² diverges but θ tends to (-λ, 0⁻), so
# Newton stays well conditioned. Start: μ = m, σ² = v.
function _canonical_params(::Type{Normal}, a, b, m, v; tol, maxiter)
    sd = sqrt(v)
    function resid(t1, t2)
        σ2 = -1 / (2t2)
        mm, vv = _tls_moments(_NormalKey(), t1 * σ2, sqrt(σ2), a, b)
        return (mm - m) / sd, (vv - v) / v
    end
    θ1, θ2, ok = _newton2(resid, m / v, -1 / (2v), (y1, y2) -> y2 < 0;
                          tol = tol, maxiter = maxiter)
    s2 = -1 / (2θ2)   # (not σ2: the closure above would then share it)
    return θ1 * s2, sqrt(s2), ok
end

# Student-t on a canonical domain, in (μ, log σ). Start: μ = m,
# σ² = (ν - 2)/ν · v (σ² = v when ν ≤ 2).
function _canonical_params(d::TDist, a, b, m, v; tol, maxiter)
    ν = dof(d)
    sd = sqrt(v)
    function resid(x1, x2)
        mm, vv = _tls_moments(ν, x1, exp(x2), a, b)
        return (mm - m) / sd, (vv - v) / v
    end
    σ0 = ν > 2 ? sqrt((ν - 2) / ν * v) : sd
    x1, x2, ok = _newton2(resid, m, log(σ0), (y1, y2) -> 0 < exp(y2) < Inf;
                          tol = tol, maxiter = maxiter)
    return x1, exp(x2), ok
end

# Known feasibility conditions. Returns `nothing` or a reason string.
function _feasibility(::Type{Normal}, a, b, m, v)
    a < m < b || return "the mean must lie strictly inside the interval"
    v > 0 || return "the variance must be positive"
    if isfinite(a) && isfinite(b)
        vmax = _truncexp_max_var(a, b, m)
        v < vmax || return "the variance must be below the truncated-exponential " *
                           "bound v_I(mean) = $vmax (the boundary is not attained)"
    elseif isfinite(a)
        v < (m - a)^2 || return "the standard deviation must be below mean - lower = $(m - a)"
    elseif isfinite(b)
        v < (b - m)^2 || return "the standard deviation must be below upper - mean = $(b - m)"
    end
    return nothing
end

function _feasibility(d::TDist, a, b, m, v)
    ν = dof(d)
    a < m < b || return "the mean must lie strictly inside the interval"
    v > 0 || return "the variance must be positive"
    if isfinite(a) && isfinite(b)
        # The exact two-sided region is open (GitHub issue #1); Bhatia–Davis
        # is necessary for every law on [a, b].
        v < (b - m) * (m - a) ||
            return "the variance must be below the Bhatia–Davis bound (upper - mean)(mean - lower)"
        return nothing
    end
    ν > 2 || return "an unbounded interval needs ν > 2 for a finite variance"
    cap = ν / (ν - 2)
    if isfinite(a)
        v < cap * (m - a)^2 ||
            return "the variance must be below ν/(ν - 2)·(mean - lower)² = $(cap * (m - a)^2)"
    elseif isfinite(b)
        v < cap * (b - m)^2 ||
            return "the variance must be below ν/(ν - 2)·(upper - mean)² = $(cap * (b - m)^2)"
    end
    return nothing
end

"""
    truncated_locscale_params(family, a, b, m, v; tol=1e-11, maxiter=100) -> (μ, σ)

Location `μ` and scale `σ` such that `μ + σ X`, with `X` standard Normal
(`family = Normal`) or Student-t (`family = TDist(ν)`, `ν > 0`), truncated to
`[a, b]`, has mean `m` and variance `v`. Either bound may be infinite. Returns
plain numbers, for use in inner loops; `make_dist` builds the distribution.

The problem is shifted and scaled to a canonical interval (`[-1/2, 1/2]`,
`[0, ∞)`, or reflected), and solved by damped Newton iteration with an exact
Jacobian (ForwardDiff through the closed-form moments of
[`truncated_locscale_moments`](@ref)). The Normal is solved in its natural
parameters, which keeps Newton well conditioned up to the boundary of the
feasible region. `tol` is on the residuals relative to `√v` (mean) and `v`
(variance); if the line search stalls (very close to the boundary of the
feasible region, where the moments themselves carry roundoff) a residual
below `1000 tol` is accepted.

Throws [`InfeasibleMomentsError`](@ref) if the target violates a known
feasibility condition, or, for the two-sided Student-t (whose feasible region
is not yet characterised), if Newton does not converge.
"""
truncated_locscale_params(::Type{Normal}, a::Real, b::Real, m::Real, v::Real;
                          tol::Real = 1e-11, maxiter::Int = 100) =
    _tls_params(Normal, a, b, m, v, Float64(tol), maxiter)
truncated_locscale_params(d::TDist, a::Real, b::Real, m::Real, v::Real;
                          tol::Real = 1e-11, maxiter::Int = 100) =
    _tls_params(d, a, b, m, v, Float64(tol), maxiter)

function _tls_params(fam::F, a::Real, b::Real, m::Real, v::Real,
                     tol::Float64, maxiter::Int) where {F}
    a < b || throw(ArgumentError("need a < b, got a = $a, b = $b"))
    a, b, m, v = Float64(a), Float64(b), Float64(m), Float64(v)
    reason = _feasibility(fam, a, b, m, v)
    reason === nothing || _infeasible(_famname(fam), a, b, m, v, reason)

    if !isfinite(a) && !isfinite(b)
        fam === Normal && return m, sqrt(v)
        ν = dof(fam)
        return m, sqrt((ν - 2) / ν * v)
    end

    # Canonical problem, solved, then mapped back.
    if isfinite(a) && isfinite(b)
        c, w = (a + b) / 2, b - a
        μc, σc, ok = _canonical_params(fam, -0.5, 0.5, (m - c) / w, v / w^2;
                                       tol = tol, maxiter = maxiter)
        μ, σ = c + w * μc, w * σc
    elseif isfinite(a)
        sd = sqrt(v)
        μc, σc, ok = _canonical_params(fam, 0.0, Inf, (m - a) / sd, 1.0;
                                       tol = tol, maxiter = maxiter)
        μ, σ = a + sd * μc, sd * σc
    else
        # (-∞, b]: reflect to [-b, ∞) (both families are symmetric).
        sd = sqrt(v)
        μc, σc, ok = _canonical_params(fam, 0.0, Inf, (b - m) / sd, 1.0;
                                       tol = tol, maxiter = maxiter)
        μ, σ = b - sd * μc, sd * σc
    end

    ok || _infeasible(_famname(fam), a, b, m, v,
        "Newton did not converge within $maxiter iterations" *
        (fam isa TDist && isfinite(a) && isfinite(b) ?
            "; the two-sided t feasible region is not yet characterised, " *
            "so the target is presumably outside it" :
            "; the target is close to the boundary of the feasible region"),
        proven = false)
    return μ, σ
end
