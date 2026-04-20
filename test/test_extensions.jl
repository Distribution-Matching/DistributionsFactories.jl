# Tests for distributions defined in src/extensions/ and the factory paths
# that produce them.

using DistributionsFactories: _truncated_poisson_moments
using QuadGK
using Random

# --- FoldedNormal -----------------------------------------------------------

function test_folded_normal_interface()
    d = FoldedNormal(2.0, 1.0)
    isapprox(mean(d), 2.0169814; atol=1e-6) || return false
    isapprox(var(d), 0.9317860; atol=1e-6) || return false
    insupport(d, 1.0) || return false
    !insupport(d, -0.5) || return false
    pdf(d, -1.0) == 0 || return false
    cdf(d, 0.0) == 0 || return false
    isapprox(quantile(d, cdf(d, 1.5)), 1.5; atol=1e-6) || return false
    return true
end

function test_folded_normal_half_normal()
    # μ = 0 reduces to half-normal: mean = σ √(2/π), var = σ²(1 - 2/π)
    σ = 1.5
    d = FoldedNormal(0.0, σ)
    isapprox(mean(d), σ * √(2/π); atol=1e-10) || return false
    isapprox(var(d), σ^2 * (1 - 2/π); atol=1e-10) || return false
    return true
end

function test_folded_normal_pdf_integrates_to_one()
    d = FoldedNormal(1.5, 0.7)
    I, _ = quadgk(x -> pdf(d, x), 0.0, 30.0)
    return isapprox(I, 1.0; atol=1e-6)
end

function test_folded_normal_factory_roundtrip()
    d = make_dist(FoldedNormal, mean=2.5, var=1.2)
    isa(d, FoldedNormal) || return false
    isapprox(mean(d), 2.5; atol=1e-6) || return false
    isapprox(var(d), 1.2; atol=1e-6) || return false
    return true
end

function test_folded_normal_rand_in_support()
    rng = MersenneTwister(42)
    d = FoldedNormal(0.5, 1.0)
    samples = [rand(rng, d) for _ in 1:1000]
    return all(s -> s ≥ 0, samples)
end

# --- DiscreteSymmetricTriangular --------------------------------------------

function test_dst_closed_form_moments()
    for n in 0:6, μ in (-3, 0, 7)
        d = DiscreteSymmetricTriangular(μ, n)
        # Compute mean/var by direct sum from PMF.
        pts = (μ - n):(μ + n)
        ps = [pdf(d, k) for k in pts]
        isapprox(sum(ps), 1.0; atol=1e-12) || return false
        m_emp = sum(k * pdf(d, k) for k in pts)
        v_emp = sum((k - m_emp)^2 * pdf(d, k) for k in pts)
        isapprox(mean(d), m_emp; atol=1e-12) || return false
        isapprox(var(d), v_emp; atol=1e-12) || return false
    end
    return true
end

function test_dst_factory_feasible()
    # n=3 → var = 3·5/6 = 2.5
    d = make_dist(DiscreteSymmetricTriangular, mean=4, var=2.5)
    return d == DiscreteSymmetricTriangular(4, 3)
end

function test_dst_factory_rejects_nonint_n()
    try
        make_dist(DiscreteSymmetricTriangular, mean=0, var=1.0)
        return false
    catch e
        return e isa DomainError
    end
end

function test_dst_quantile_inverse_cdf()
    d = DiscreteSymmetricTriangular(0, 4)
    for k in -4:4
        cd = cdf(d, k)
        quantile(d, cd) == k || return false
    end
    return true
end

# --- DiscreteTriangular -----------------------------------------------------

function test_disc_triangular_pmf_normalises()
    d = DiscreteTriangular(2, 8, 5)
    return isapprox(sum(pdf(d, k) for k in 2:8), 1.0; atol=1e-12)
end

function test_disc_triangular_mode_is_max_pmf()
    d = DiscreteTriangular(0, 10, 7)
    return argmax([pdf(d, k) for k in 0:10]) - 1 == 7  # 0-indexed → value 7
end

function test_disc_triangular_factory_close_to_target()
    # Tolerate approximation since (a, b) get rounded to integers.
    d = make_dist(DiscreteTriangular, mean=5.0, var=2.0, mode=5)
    isa(d, DiscreteTriangular) || return false
    isapprox(mean(d), 5.0; atol=0.5) || return false
    isapprox(var(d), 2.0; atol=1.0) || return false
    return true
end

# --- TriangularDist mean+var+mode -------------------------------------------

function test_triangular_dist_factory_exact()
    d = make_dist(TriangularDist, mean=5.0, var=2.0, mode=5.0)
    isa(d, TriangularDist) || return false
    isapprox(mean(d), 5.0; atol=1e-10) || return false
    isapprox(var(d), 2.0; atol=1e-10) || return false
    isapprox(d.c, 5.0; atol=1e-10) || return false
    return true
end

function test_triangular_dist_factory_asymmetric()
    # Mode off-centre.
    d = make_dist(TriangularDist, mean=4.0, var=1.5, mode=3.0)
    return isa(d, TriangularDist) &&
           isapprox(mean(d), 4.0; atol=1e-10) &&
           isapprox(var(d), 1.5; atol=1e-10) &&
           d.a ≤ 3.0 ≤ d.b
end

function test_triangular_dist_factory_infeasible_mode()
    # Mode outside the resulting (a, b) range → DomainError.
    try
        make_dist(TriangularDist, mean=0.0, var=10.0, mode=20.0)
        return false
    catch e
        return e isa DomainError
    end
end

# --- Truncated{Poisson} factory ---------------------------------------------

function test_truncated_poisson_mean_only()
    template = truncated(Poisson(0.0); lower=0, upper=10)
    d = DistributionsFactories.dist_from_mean(template, 3.0)
    m, _ = _truncated_poisson_moments(d.untruncated.λ, 0.0, 10.0)
    return isapprox(m, 3.0; atol=1e-6)
end

function test_truncated_poisson_mean_var_consistent()
    template = truncated(Poisson(0.0); lower=0, upper=10)
    m, v = _truncated_poisson_moments(2.5, 0.0, 10.0)  # consistent target
    d = make_dist(template, mean=m, var=v)
    return isapprox(d.untruncated.λ, 2.5; rtol=1e-4)
end

function test_truncated_poisson_inconsistent_var_throws()
    template = truncated(Poisson(0.0); lower=0, upper=10)
    try
        make_dist(template, mean=3.0, var=10.0)  # var not achievable
        return false
    catch e
        return e isa ArgumentError
    end
end

function test_truncated_poisson_mean_outside_bounds_rejected()
    template = truncated(Poisson(0.0); lower=2, upper=8)
    return !exists_dist_from_mean_var(template, 1.0, 1.0) &&
           !exists_dist_from_mean_var(template, 9.0, 1.0)
end

# --- Half-truncated feasibility (exponential-tail families) -----------------

function test_half_trunc_normal_lower_feasible()
    d = truncated(Normal(0,1); lower=0)         # support [0, ∞)
    # σ̄² < (μ̄ - 0)² = μ̄² should be accepted
    exists_dist_from_mean_var(d, 2.0, 1.0)  || return false   # CV = 0.5 ✓
    exists_dist_from_mean_var(d, 5.0, 4.0)  || return false   # CV = 0.4 ✓
    return true
end

function test_half_trunc_normal_lower_above_bound_rejected()
    d = truncated(Normal(0,1); lower=0)
    # σ̄² ≥ μ̄² should be rejected (CV ≥ 1)
    exists_dist_from_mean_var(d, 2.0, 4.0)  && return false   # CV = 1
    exists_dist_from_mean_var(d, 2.0, 5.0)  && return false   # CV > 1
    return true
end

function test_half_trunc_normal_lower_with_offset()
    d = truncated(Normal(0,1); lower=3)         # support [3, ∞)
    # Bound becomes σ̄² < (μ̄ - 3)²
    exists_dist_from_mean_var(d, 5.0, 3.0) || return false    # 3 < 4 ✓
    exists_dist_from_mean_var(d, 5.0, 5.0) && return false    # 5 ≥ 4 ✗
    exists_dist_from_mean_var(d, 2.5, 1.0) && return false    # μ̄ < lo
    return true
end

function test_half_trunc_normal_upper_feasible()
    d = truncated(Normal(0,1); upper=0)         # support (-∞, 0]
    # σ̄² < (0 - μ̄)² = μ̄²
    exists_dist_from_mean_var(d, -2.0, 1.0) || return false   # gap=2, var<4 ✓
    exists_dist_from_mean_var(d, -2.0, 4.0) && return false   # var = gap² ✗
    exists_dist_from_mean_var(d,  1.0, 1.0) && return false   # μ̄ > hi
    return true
end

function test_half_trunc_laplace_logistic()
    for D in (Laplace, Logistic)
        d_lo = truncated(D(); lower=0)
        d_hi = truncated(D(); upper=0)
        exists_dist_from_mean_var(d_lo,  2.0, 1.0) || return false
        exists_dist_from_mean_var(d_hi, -2.0, 1.0) || return false
        # Above the bound: definitely infeasible for both
        exists_dist_from_mean_var(d_lo,  2.0, 5.0) && return false
        exists_dist_from_mean_var(d_hi, -2.0, 5.0) && return false
    end
    # Boundary semantics differ between the two families:
    # - Laplace attains σ̄² = gap² exactly (parent μp ≤ lo gives Exponential)
    # - Logistic only approaches it asymptotically
    exists_dist_from_mean_var(truncated(Laplace();  lower=0), 2.0, 4.0)  || return false
    exists_dist_from_mean_var(truncated(Logistic(); lower=0), 2.0, 4.0)  && return false
    return true
end

function test_half_trunc_locscale_two_sided_still_uses_dome()
    # Make sure adding the half-truncated branch didn't break the existing
    # two-sided Langevin envelope.
    d = truncated(Normal(0,1); lower=-0.5, upper=0.5)
    return exists_dist_from_mean_var(d, 0.0, 0.05) &&        # inside dome
           !exists_dist_from_mean_var(d, 0.0, 0.5)           # above dome
end

# --- Truncated Student-t feasibility ---------------------------------------

function test_trunc_tdist_low_dof_rejected()
    d = truncated(TDist(2.0); lower=0)          # variance undefined for ν ≤ 2
    return !exists_dist_from_mean_var(d, 1.0, 1.0)
end

function test_trunc_tdist_pareto_bound()
    # ν = 4: ceiling is √(ν/(ν-2)) = √2 ≈ 1.414
    d = truncated(TDist(4.0); lower=0)
    # σ̄² < ν/(ν-2)·μ̄² = 2·μ̄²
    exists_dist_from_mean_var(d, 2.0, 7.5) || return false   # 7.5 < 8 ✓
    exists_dist_from_mean_var(d, 2.0, 8.5) && return false   # 8.5 > 8 ✗
    return true
end

function test_trunc_tdist_recovers_normal_limit()
    # Large ν → exponential-tail bound (CV → 1)
    d = truncated(TDist(1000.0); lower=0)
    exists_dist_from_mean_var(d, 2.0, 3.9) || return false    # well below μ̄²
    exists_dist_from_mean_var(d, 2.0, 4.1) && return false    # just above
    return true
end

function test_trunc_tdist_upper_side()
    d = truncated(TDist(5.0); upper=0)
    # ν = 5: cap = 5/3
    exists_dist_from_mean_var(d, -3.0, 14.0) || return false  # < 5/3·9 = 15
    exists_dist_from_mean_var(d, -3.0, 16.0) && return false  # > 15
    return true
end

# --- Weibull / Frechet constructor roundtrip ---

function test_weibull_constructor_roundtrip()
    # Sweep parents and verify dist_from_mean_var recovers them.
    for k in (0.7, 1.0, 1.5, 2.5, 4.0), λ in (0.5, 1.0, 3.0)
        parent = Weibull(k, λ)
        μ̄, σ̄² = mean(parent), var(parent)
        d = dist_from_mean_var(Weibull, μ̄, σ̄²)
        isapprox(mean(d), μ̄; rtol=1e-4) || return false
        isapprox(var(d), σ̄²; rtol=1e-4) || return false
    end
    return true
end

function test_frechet_constructor_roundtrip()
    # Frechet: shape α > 2 needed for variance.
    for α in (2.5, 3.0, 5.0, 10.0), s in (0.5, 1.0, 3.0)
        parent = Frechet(α, s)
        μ̄, σ̄² = mean(parent), var(parent)
        d = dist_from_mean_var(Frechet, μ̄, σ̄²)
        isapprox(mean(d), μ̄; rtol=1e-4) || return false
        isapprox(var(d), σ̄²; rtol=1e-4) || return false
    end
    return true
end

# --- InverseGamma μ̄ > 0 enforcement ---

function test_inverse_gamma_negative_mean_rejected()
    return !exists_dist_from_mean_var(InverseGamma, -1.0, 1.0)
end

# --- Truncated{<:Poisson} non-integer bound rejection ---

function test_truncated_poisson_noninteger_bounds_rejected()
    template = truncated(Poisson(0.0); lower=0.5, upper=10)
    try
        make_dist(template, mean=3.0, var=2.0)
        return false
    catch e
        return e isa ArgumentError
    end
end

# --- Truncated{<:TDist} factory (half-truncated only) -----------------------

function _trunc_tdist_half_moments(parent, lo, hi)
    if isfinite(lo)
        Z, _ = quadgk(x -> pdf(parent, x), lo, Inf; rtol=1e-10)
        m, _ = quadgk(x -> x*pdf(parent, x), lo, Inf; rtol=1e-10)
        m2,_ = quadgk(x -> x^2*pdf(parent, x), lo, Inf; rtol=1e-10)
    else
        Z, _ = quadgk(x -> pdf(parent, x), -Inf, hi; rtol=1e-10)
        m, _ = quadgk(x -> x*pdf(parent, x), -Inf, hi; rtol=1e-10)
        m2,_ = quadgk(x -> x^2*pdf(parent, x), -Inf, hi; rtol=1e-10)
    end
    return m/Z, m2/Z - (m/Z)^2
end

function test_trunc_tdist_factory_half_below_roundtrip()
    for ν in (2.5, 3.0, 5.0, 10.0), μp in (-1.0, 0.0, 0.5, 2.0)
        parent = μp + 1.0 * TDist(ν)
        μ̄, σ̄² = _trunc_tdist_half_moments(parent, 0.0, Inf)
        template = truncated(TDist(ν); lower=0)
        d = dist_from_mean_var(template, μ̄, σ̄²)
        μ_back, σ²_back = _trunc_tdist_half_moments(d.untruncated, 0.0, Inf)
        isapprox(μ_back, μ̄; rtol=1e-4)  || return false
        isapprox(σ²_back, σ̄²; rtol=1e-4) || return false
    end
    return true
end

function test_trunc_tdist_factory_half_above_roundtrip()
    for ν in (3.0, 5.0), μp in (-2.0, -0.5, 1.0)
        parent = μp + 1.0 * TDist(ν)
        μ̄, σ̄² = _trunc_tdist_half_moments(parent, -Inf, 0.0)
        template = truncated(TDist(ν); upper=0)
        d = dist_from_mean_var(template, μ̄, σ̄²)
        μ_back, σ²_back = _trunc_tdist_half_moments(d.untruncated, -Inf, 0.0)
        isapprox(μ_back, μ̄; rtol=1e-4)  || return false
        isapprox(σ²_back, σ̄²; rtol=1e-4) || return false
    end
    return true
end

function test_trunc_tdist_factory_two_sided_not_implemented()
    template = truncated(TDist(5); lower=-1, upper=1)
    try
        dist_from_mean_var(template, 0.0, 0.5)
        return false
    catch e
        return e isa ArgumentError
    end
end
