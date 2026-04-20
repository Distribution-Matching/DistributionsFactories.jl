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
