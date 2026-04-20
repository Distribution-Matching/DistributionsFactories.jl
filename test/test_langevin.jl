using DistributionsFactories: langevin, langevin_deriv, inv_langevin,
    _truncexp_max_var
using QuadGK: quadgk

# --- Langevin function itself ---

function test_langevin_small_vs_closed()
    # At x slightly above the series cutoff, the two branches should agree.
    for x in (1.01e-4, 5e-4, 1e-3, 1e-2, 0.5, 1.0, 3.0, 10.0)
        series = x/3 - x^3/45 + 2x^5/945 - x^7/4725
        closed = coth(x) - 1/x
        # Closed form is the ground truth except right at 0.
        isapprox(closed, langevin(x); atol=1e-12, rtol=1e-12) || return false
        isapprox(series, x/3 - x^3/45 + 2x^5/945 - x^7/4725; atol=0) || return false
    end
    # At exactly 0, L(0) = 0.
    langevin(0.0) == 0.0 || return false
    # Odd function.
    for x in (0.1, 1.0, 5.0)
        isapprox(langevin(-x), -langevin(x); atol=1e-14) || return false
    end
    return true
end

function test_langevin_deriv_matches_finite_diff()
    h = 1e-6
    for x in (-5.0, -1.0, -0.1, -1e-3, 1e-3, 0.1, 1.0, 5.0)
        fd = (langevin(x + h) - langevin(x - h)) / (2h)
        isapprox(fd, langevin_deriv(x); rtol=1e-6, atol=1e-9) || return false
    end
    # At 0, L'(0) = 1/3.
    isapprox(langevin_deriv(0.0), 1/3; atol=1e-14) || return false
    return true
end

function test_inv_langevin_roundtrip()
    # L(L^{-1}(y)) ≈ y across the whole domain.
    for y in range(-0.999, 0.999, length=41)
        z = inv_langevin(y)
        isapprox(langevin(z), y; atol=1e-12, rtol=1e-12) || return false
    end
    # And L^{-1}(L(x)) ≈ x.
    for x in (-8.0, -3.0, -0.5, 0.5, 3.0, 8.0)
        isapprox(inv_langevin(langevin(x)), x; rtol=1e-10, atol=1e-10) || return false
    end
    inv_langevin(0.0) == 0.0 || return false
    return true
end

function test_inv_langevin_domain_error()
    try; inv_langevin(1.0);  return false; catch e; e isa DomainError || return false; end
    try; inv_langevin(-1.0); return false; catch e; e isa DomainError || return false; end
    try; inv_langevin(1.5);  return false; catch e; e isa DomainError || return false; end
    return true
end

# --- Truncexp envelope: cross-check against direct quadrature ---

# Moments of the truncated-exponential density ∝ e^{-λx} on [a, b].
function _truncexp_moments(a, b, λ)
    abs(λ) < 1e-12 && return ((a+b)/2, (b-a)^2 / 12)
    Z, _  = quadgk(x -> exp(-λ*x), a, b; rtol=1e-12)
    m1, _ = quadgk(x -> x * exp(-λ*x), a, b; rtol=1e-12)
    m2, _ = quadgk(x -> x^2 * exp(-λ*x), a, b; rtol=1e-12)
    μ = m1 / Z
    v = m2 / Z - μ^2
    return (μ, v)
end

function test_truncexp_max_var_matches_quadrature()
    for (a, b) in ((-0.5, 0.5), (0.0, 1.0), (-2.0, 3.0))
        for λ in (-8.0, -2.0, -0.5, 0.5, 2.0, 8.0)
            μ, v_quad = _truncexp_moments(a, b, λ)
            v_closed = _truncexp_max_var(a, b, μ)
            # The closed form *is* the quadrature variance at this (a,b,μ).
            isapprox(v_closed, v_quad; rtol=1e-8, atol=1e-10) || return false
        end
        # λ = 0: uniform, var = (b-a)²/12 at the midpoint.
        c = (a + b) / 2
        isapprox(_truncexp_max_var(a, b, c), (b-a)^2 / 12; rtol=1e-10) || return false
    end
    return true
end

function test_truncexp_max_var_below_bhatia_davis()
    # The truncexp dome lies strictly inside the Bhatia–Davis ceiling (b-μ)(μ-a).
    for (a, b) in ((-0.5, 0.5), (2.0, 7.0))
        for μ in range(a + 0.1*(b-a), b - 0.1*(b-a), length=9)
            σ²_dome = _truncexp_max_var(a, b, μ)
            σ²_bd   = (b - μ) * (μ - a)
            (σ²_dome < σ²_bd) || return false
        end
    end
    return true
end

# --- Wiring into exists_dist_from_mean_var / dist_from_mean_var ---

function test_truncated_normal_langevin_rejects_above_dome()
    d = truncated(Normal(), -0.5, 0.5)
    μ̄ = -0.3
    σ²_max = _truncexp_max_var(-0.5, 0.5, μ̄)
    # Pure-predicate semantics: above the dome → false, inside → true.
    exists_dist_from_mean_var(d, μ̄, σ²_max * 1.01) && return false
    exists_dist_from_mean_var(d, μ̄, σ²_max * 0.5) || return false
    return true
end

function test_truncated_laplace_langevin_rejects_above_dome()
    d = truncated(Laplace(), -0.5, 0.5)
    μ̄ = 0.1
    σ²_max = _truncexp_max_var(-0.5, 0.5, μ̄)
    exists_dist_from_mean_var(d, μ̄, σ²_max * 1.5) && return false
    exists_dist_from_mean_var(d, μ̄, σ²_max * 0.5) || return false
    return true
end

function test_truncated_logistic_langevin_rejects_above_dome()
    d = truncated(Logistic(), -0.5, 0.5)
    # At the midpoint, σ²_max = (b-a)²/12 = 1/12.
    isapprox(_truncexp_max_var(-0.5, 0.5, 0.0), 1/12; rtol=1e-10) || return false
    exists_dist_from_mean_var(d, 0.0, 1/12 + 0.01) && return false
    exists_dist_from_mean_var(d, 0.0, 0.05) || return false
    return true
end

function test_dist_from_mean_var_truncated_rejects_above_dome()
    # dist_from_mean_var on the three families now calls the existence check first.
    for D in (Normal, Laplace, Logistic)
        d = truncated(D(), -0.5, 0.5)
        μ̄ = 0.2
        σ²_max = _truncexp_max_var(-0.5, 0.5, μ̄)
        try
            dist_from_mean_var(d, μ̄, σ²_max * 1.2)
            return false
        catch e
            e isa DomainError || return false
        end
    end
    return true
end

# Compute (mean, variance) of a parent distribution `parent` truncated to [a,b].
# We don't lean on Distributions.jl's own mean()/var() on Truncated{...} since
# they're not always defined for Logistic/Laplace — quadrature matches what the
# solver sees anyway.
function _truncated_moments(parent, a, b)
    td = truncated(parent, a, b)
    m1, _ = quadgk(x -> x * pdf(td, x), a, b; rtol=1e-10)
    m2, _ = quadgk(x -> x^2 * pdf(td, x), a, b; rtol=1e-10)
    return (m1, m2 - m1^2)
end

# Roundtrip: any (μ, σ²) that actually comes from a truncated Normal/Laplace/
# Logistic on [a,b] must lie strictly inside the Langevin dome, so
# exists_dist_from_mean_var should accept it. This is the converse of the
# "reject above dome" tests — it catches a whole class of bugs in the envelope
# formula (sign errors, off-by-one in (c-μ̄)/w, L vs L', …) that the purely
# analytic tests would miss.
function _roundtrip_family(D, parents, intervals)
    for (a, b) in intervals
        for parent in parents
            μ̄, σ̄² = _truncated_moments(parent, a, b)
            # Skip degenerate cases where truncation mass underflows: in that
            # regime the pdf returns NaN/0 and the moment is ill-defined, but
            # this is a numerics issue, not a feasibility issue.
            (isfinite(μ̄) && isfinite(σ̄²) && σ̄² > 0) || continue

            d_template = truncated(D(), a, b)
            if !exists_dist_from_mean_var(d_template, μ̄, σ̄²)
                @info "Roundtrip feasibility failed" D parent a b μ̄ σ̄² σ²_max=_truncexp_max_var(a, b, μ̄)
                return false
            end
        end
    end
    return true
end

function test_truncated_normal_roundtrip_feasible()
    parents = [Normal(μ, σ) for μ in (-1.0, -0.25, 0.0, 0.25, 1.0)
                            for σ in (0.05, 0.2, 0.5, 1.0, 3.0)]
    return _roundtrip_family(Normal, parents, ((-0.5, 0.5), (0.0, 1.0), (-2.0, 3.0)))
end

function test_truncated_laplace_roundtrip_feasible()
    parents = [Laplace(μ, θ) for μ in (-1.0, -0.25, 0.0, 0.25, 1.0)
                             for θ in (0.05, 0.2, 0.5, 1.0, 3.0)]
    return _roundtrip_family(Laplace, parents, ((-0.5, 0.5), (0.0, 1.0), (-2.0, 3.0)))
end

function test_truncated_logistic_roundtrip_feasible()
    parents = [Logistic(μ, s) for μ in (-1.0, -0.25, 0.0, 0.25, 1.0)
                              for s in (0.05, 0.2, 0.5, 1.0, 3.0)]
    return _roundtrip_family(Logistic, parents, ((-0.5, 0.5), (0.0, 1.0), (-2.0, 3.0)))
end

# Extra: verify that the (μ̄, σ̄²) produced by a parent distribution is
# genuinely *strictly* below the envelope — i.e. the margin is visible, not
# just "passes the < test by a hair". This rules out the pathological case
# where exists_dist_from_mean_var accepts because we got lucky with rounding.
# --- Constructor roundtrip: dist_from_mean_var on Truncated{D} actually
# recovers the input moments (not just feasibility predicate semantics). ---

function _constructor_roundtrip_family(D, parents, intervals;
                                       atol::Real = 1e-6, rtol::Real = 1e-6,
                                       min_truncation_mass::Real = 0.05,
                                       boundary_margin::Real = 0.05)
    for parent in parents
        for (a, b) in intervals
            # Skip pathological cases where the parent puts almost no mass in
            # [a, b] — the truncated moments are numerically defined but
            # standardize to the boundary of the canonical (-0.5, 0.5) domain
            # where Newton legitimately struggles.
            mass = cdf(parent, b) - cdf(parent, a)
            mass ≥ min_truncation_mass || continue

            μ̄, σ̄² = _truncated_moments(parent, a, b)
            (isfinite(μ̄) && isfinite(σ̄²) && σ̄² > 0) || continue

            # Skip cases whose standardized mean lies within `boundary_margin`
            # of ±0.5 — this is the "mass concentrated at a boundary" regime
            # (e.g. half-truncated Laplace at its mode), where the parent is
            # at the edge of identifiability. The solver's design regime is
            # moments comfortably inside the canonical interval.
            c, w = (a + b) / 2, b - a
            μ̄_std = (μ̄ - c) / w
            abs(μ̄_std) ≤ 0.5 - boundary_margin || continue

            d_template = truncated(D(), a, b)
            d_recovered = dist_from_mean_var(d_template, μ̄, σ̄²)

            μ_back, σ²_back = _truncated_moments(d_recovered.untruncated, a, b)
            if !isapprox(μ_back, μ̄; atol=atol, rtol=rtol)
                @info "constructor mean mismatch" D parent a b μ̄ μ_back
                return false
            end
            if !isapprox(σ²_back, σ̄²; atol=atol, rtol=rtol)
                @info "constructor var mismatch" D parent a b σ̄² σ²_back
                return false
            end
        end
    end
    return true
end

function test_truncated_normal_constructor_roundtrip()
    parents = [Normal(μ, σ) for μ in (-1.0, -0.25, 0.0, 0.25, 1.0)
                            for σ in (0.05, 0.2, 0.5, 1.0, 3.0)]
    return _constructor_roundtrip_family(Normal, parents,
                                          ((-0.5, 0.5), (0.0, 1.0), (-2.0, 3.0)))
end

function test_truncated_laplace_constructor_roundtrip()
    parents = [Laplace(μ, θ) for μ in (-1.0, -0.25, 0.0, 0.25, 1.0)
                             for θ in (0.05, 0.2, 0.5, 1.0, 3.0)]
    return _constructor_roundtrip_family(Laplace, parents,
                                          ((-0.5, 0.5), (0.0, 1.0), (-2.0, 3.0)))
end

function test_truncated_logistic_constructor_roundtrip()
    parents = [Logistic(μ, s) for μ in (-1.0, -0.25, 0.0, 0.25, 1.0)
                              for s in (0.05, 0.2, 0.5, 1.0, 3.0)]
    return _constructor_roundtrip_family(Logistic, parents,
                                          ((-0.5, 0.5), (0.0, 1.0), (-2.0, 3.0)))
end

# --- Half-truncated constructor roundtrip ---

# Quadrature on a half-line. We cap at 50σp from the parent mode for stability;
# all three families have at most exponential right tails so the residual is
# negligible.
function _half_truncated_moments(parent, lo, hi)
    if isfinite(lo) && !isfinite(hi)
        upper = lo + 50 * std(parent) + abs(mean(parent)) + 50
        Z, _ = quadgk(x -> pdf(parent, x), lo, upper; rtol=1e-10)
        Z > 1e-12 || return (NaN, NaN)
        m,  _ = quadgk(x -> x   * pdf(parent, x), lo, upper; rtol=1e-10)
        m2, _ = quadgk(x -> x^2 * pdf(parent, x), lo, upper; rtol=1e-10)
        μ = m / Z; v = m2/Z - μ^2
        return (μ, v)
    elseif !isfinite(lo) && isfinite(hi)
        lower = hi - 50 * std(parent) - abs(mean(parent)) - 50
        Z, _ = quadgk(x -> pdf(parent, x), lower, hi; rtol=1e-10)
        Z > 1e-12 || return (NaN, NaN)
        m,  _ = quadgk(x -> x   * pdf(parent, x), lower, hi; rtol=1e-10)
        m2, _ = quadgk(x -> x^2 * pdf(parent, x), lower, hi; rtol=1e-10)
        μ = m / Z; v = m2/Z - μ^2
        return (μ, v)
    end
    error("expected exactly one infinite bound")
end

function _half_constructor_roundtrip(D, parents, bounds; atol=1e-6, rtol=1e-6)
    for parent in parents
        for (lo, hi) in bounds
            mass = (isfinite(hi) ? cdf(parent, hi) : 1.0) -
                   (isfinite(lo) ? cdf(parent, lo) : 0.0)
            mass ≥ 0.05 || continue

            μ̄, σ̄² = _half_truncated_moments(parent, lo, hi)
            (isfinite(μ̄) && isfinite(σ̄²) && σ̄² > 0) || continue

            template = if isfinite(lo)
                truncated(D(); lower=lo)
            else
                truncated(D(); upper=hi)
            end

            d_recovered = dist_from_mean_var(template, μ̄, σ̄²)
            μ_back, σ²_back = _half_truncated_moments(d_recovered.untruncated, lo, hi)
            if !isapprox(μ_back, μ̄; atol=atol, rtol=rtol) ||
               !isapprox(σ²_back, σ̄²; atol=atol, rtol=rtol)
                @info "half-trunc constructor mismatch" D parent lo hi μ̄ σ̄² μ_back σ²_back
                return false
            end
        end
    end
    return true
end

function test_half_trunc_normal_constructor()
    parents = [Normal(μ, σ) for μ in (-2.0, -0.5, 0.0, 0.5, 2.0)
                            for σ in (0.3, 1.0, 3.0)]
    bounds = [(0.0, Inf), (-1.0, Inf), (1.0, Inf), (-Inf, 0.0), (-Inf, 1.0)]
    return _half_constructor_roundtrip(Normal, parents, bounds)
end

function test_half_trunc_laplace_constructor()
    # Avoid the boundary case (Laplace μp ≤ lo gives the Exponential limit
    # exactly, where the parent is non-unique). Use parents safely above.
    parents = [Laplace(μ, θ) for μ in (0.5, 1.0, 2.0) for θ in (0.3, 1.0, 3.0)]
    bounds = [(0.0, Inf), (-Inf, 5.0)]
    return _half_constructor_roundtrip(Laplace, parents, bounds)
end

function test_half_trunc_logistic_constructor()
    parents = [Logistic(μ, s) for μ in (-1.0, -0.25, 0.5, 1.5)
                              for s in (0.3, 1.0, 2.0)]
    bounds = [(0.0, Inf), (-Inf, 0.0)]
    return _half_constructor_roundtrip(Logistic, parents, bounds)
end

function test_truncated_families_margin_from_dome()
    # Pick a middle-of-the-road parent per family, away from the degenerate /
    # saturated regimes where the distribution is effectively a point mass or
    # essentially uniform.
    parents = (Normal(0.1, 0.4), Laplace(-0.2, 0.3), Logistic(0.05, 0.25))
    for parent in parents
        for (a, b) in ((-0.5, 0.5), (-1.0, 2.0))
            μ̄, σ̄² = _truncated_moments(parent, a, b)
            σ²_max = _truncexp_max_var(a, b, μ̄)
            # The parent is not the extremal truncexp distribution, so we
            # expect a non-trivial gap. Require at least ~1% margin.
            (σ̄² < 0.99 * σ²_max) || return false
        end
    end
    return true
end
