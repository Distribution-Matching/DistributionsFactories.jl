# Tests for truncated_locscale_moments / truncated_locscale_params
# (closed-form truncated Normal and Student-t moments, fast solver,
# InfeasibleMomentsError).

using ForwardDiff
using QuadGK

# Reference moments by adaptive quadrature in s = x - mode, with the density
# normalised at the mode, so that narrow and far-tail cases are resolved.
function _ref_tls_moments(fam, μ, σ, a, b)
    parent = fam === Normal ? Normal(μ, σ) : μ + σ * fam
    r = clamp(μ, a, b)
    lr = logpdf(parent, r)
    f(s) = exp(logpdf(parent, r + s) - lr)
    lo, hi = a - r, b - r
    if fam === Normal                      # clip where the density is < e^-80
        R = σ * sqrt(2 * 80 + ((r - μ) / σ)^2)
        lo, hi = max(lo, μ - R - r), min(hi, μ + R - r)
    end
    pts = lo < 0 < hi ? (lo, 0.0, hi) : (lo, hi)
    Z = quadgk(f, pts...; rtol=1e-12)[1]
    δ = quadgk(s -> s * f(s), pts...; rtol=1e-12)[1] / Z
    v = quadgk(s -> (s - δ)^2 * f(s), pts...; rtol=1e-12)[1] / Z
    return r + δ, v
end

const _TLS_BOUNDS = ((-1.0, 2.0), (0.0, Inf), (-Inf, 0.5), (-Inf, Inf), (5.0, 5.2))

function test_tls_moments_vs_quadgk()
    for fam in (Normal, TDist(0.7), TDist(1.5), TDist(3.0), TDist(4.5), TDist(7.5),
                TDist(30.0), TDist(1e4)),
        μ in (-3.0, -0.5, 0.0, 0.4, 2.0, 6.0), σ in (0.05, 0.3, 1.0, 4.0, 50.0),
        (a, b) in _TLS_BOUNDS
        # Unbounded intervals need ν > 2. The quadrature reference is slow on
        # heavy tails, so for the t they use ν ≥ 3 and a coarser (μ, σ) grid.
        if fam isa TDist && !(isfinite(a) && isfinite(b))
            dof(fam) ≥ 3 || continue
            (μ in (-0.5, 2.0) && σ in (0.3, 4.0)) || continue
        end
        m, v = truncated_locscale_moments(fam, μ, σ, a, b)
        mr, vr = _ref_tls_moments(fam, μ, σ, a, b)
        sd = sqrt(vr)
        if !(abs(m - mr) ≤ 1e-9 * sd && abs(v - vr) ≤ 1e-9 * vr)
            @info "moments mismatch" fam μ σ a b m mr v vr
            return false
        end
    end
    return true
end

# Regimes where the plain closed forms cancel: a truncated Normal near the
# boundary curve on [0, 1] (σ ≫ 1, μ = -λσ²), far tails, narrow intervals.
function test_tls_moments_hard_regimes()
    cases = [(Normal, -2 * 300.0^2, 300.0, 0.0, 1.0),
             (Normal, -3 * 1e3^2, 1e3, 0.0, 1.0),
             (Normal, 5.0, 0.01, 0.0, 1.0),
             (Normal, -40.0, 1.0, 0.0, Inf),
             (Normal, 0.3, 1e4, 0.0, 1.0),
             (TDist(30.0), -2 * 300.0^2, 300.0, 0.0, 1.0),
             (TDist(4.5), 400.0, 1.0, 0.0, 100.0),
             (TDist(3.0), -1e4, 1.0, 0.0, 1e-4),
             (TDist(1e3), -300.0, 1.0, 0.0, Inf)]
    for (fam, μ, σ, a, b) in cases
        m, v = truncated_locscale_moments(fam, μ, σ, a, b)
        mr, vr = _ref_tls_moments(fam, μ, σ, a, b)
        if !(abs(m - mr) ≤ 1e-9 * sqrt(vr) && abs(v - vr) ≤ 1e-9 * vr)
            @info "hard-regime mismatch" fam μ σ a b m mr v vr
            return false
        end
    end
    return true
end

# Far in the tail the half-line t tends to a Pareto law from below: Var/gap²
# increases to ν/(ν-2) as ≈ ν/(ν-2) - c/α² (c = 2.4 for ν = 3). This needs the
# t tail probability to full relative accuracy at large arguments.
function test_tls_t_far_tail_below_pareto_bound()
    for ν in (3.0, 10.0), α in (1e2, 1e3, 3.5e4, 1e6)
        m, v = truncated_locscale_moments(TDist(ν), 0.0, 1.0, α, Inf)
        excess = v / (m - α)^2 - ν / (ν - 2)
        excess < 0 || return false
    end
    m, v = truncated_locscale_moments(TDist(3.0), 0.0, 1.0, 1e3, Inf)
    return isapprox(v / (m - 1e3)^2 - 3, -2.39999444572848e-6; rtol=1e-6)
end

function test_tls_forwarddiff_vs_finite_differences()
    for fam in (Normal, TDist(3.5), TDist(12.0)),
        (μ, σ, a, b) in ((0.3, 0.7, -0.5, 0.5), (1.0, 2.0, 0.0, Inf),
                         (-2.0, 0.5, -Inf, 0.0), (-2 * 50.0^2, 50.0, 0.0, 1.0))
        f(x) = collect(truncated_locscale_moments(fam, x[1], x[2], a, b))
        J = ForwardDiff.jacobian(f, [μ, σ])
        for j in 1:2
            # Richardson-extrapolated central differences.
            h = 1e-3 * max(1.0, abs((μ, σ)[j]))
            cd(h) = (xp = [μ, σ]; xm = [μ, σ]; xp[j] += h; xm[j] -= h; (f(xp) - f(xm)) / (2h))
            fd = (4cd(h / 2) - cd(h)) / 3
            # Absolute floor: roundoff in the moments themselves, divided by
            # the parameter scale (derivatives can be ~1e-14 of the moments).
            floor = 1e-12 .* abs.(f([μ, σ])) ./ max(1.0, abs((μ, σ)[j]))
            all(abs.(J[:, j] .- fd) .≤ 1e-6 .* abs.(fd) .+ floor) || begin
                @info "derivative mismatch" fam μ σ a b j J[:, j] fd
                return false
            end
        end
    end
    return true
end

function test_tls_params_roundtrip()
    for fam in (Normal, TDist(3.0), TDist(4.5), TDist(30.0)),
        μ in (-2.0, -0.3, 0.4, 1.5), σ in (0.2, 1.0, 3.0),
        (a, b) in ((-1.0, 2.0), (0.0, Inf), (-Inf, 0.5))
        m, v = truncated_locscale_moments(fam, μ, σ, a, b)
        μh, σh = truncated_locscale_params(fam, a, b, m, v)
        if !(isapprox(μh, μ; atol=1e-7 * σ) && isapprox(σh, σ; rtol=1e-7))
            @info "roundtrip mismatch" fam μ σ a b μh σh
            return false
        end
    end
    return true
end

# Targets close to the boundary of the feasible region (Normal).
function test_tls_params_near_boundary()
    for mfrac in (0.02, 0.2, 0.5, 0.8, 0.98), q in (0.99, 0.9999, 0.999999)
        a, b = 0.0, 1.0
        m = a + mfrac * (b - a)
        v = q * DistributionsFactories._truncexp_max_var(a, b, m)
        μ, σ = truncated_locscale_params(Normal, a, b, m, v)
        mm, vv = truncated_locscale_moments(Normal, μ, σ, a, b)
        (abs(mm - m) ≤ 1e-9 * sqrt(v) && abs(vv - v) ≤ 1e-9 * v) || return false
    end
    for z in 1 .+ (1e-6, 1e-3, 0.1, 10.0)         # half-line, σ̄ = 1, mean z
        μ, σ = truncated_locscale_params(Normal, 0.0, Inf, z, 1.0)
        mm, vv = truncated_locscale_moments(Normal, μ, σ, 0.0, Inf)
        (abs(mm - z) ≤ 1e-9 && abs(vv - 1) ≤ 1e-9) || return false
    end
    return true
end

function test_tls_large_nu_matches_normal()
    for (μ, σ, a, b) in ((0.3, 0.7, -0.5, 0.5), (1.0, 2.0, 0.0, Inf), (-1.0, 1.0, -Inf, 0.0))
        mN, vN = truncated_locscale_moments(Normal, μ, σ, a, b)
        mT, vT = truncated_locscale_moments(TDist(1e8), μ, σ, a, b)
        (abs(mT - mN) ≤ 1e-6 * sqrt(vN) && abs(vT - vN) ≤ 1e-6 * vN) || return false
        μN, σN = truncated_locscale_params(Normal, a, b, mN, vN)
        μT, σT = truncated_locscale_params(TDist(1e8), a, b, mN, vN)
        (abs(μT - μN) ≤ 1e-5 * σN && abs(σT - σN) ≤ 1e-5 * σN) || return false
    end
    return true
end

_throws_infeasible(f) = try
    f(); false
catch e
    e isa InfeasibleMomentsError
end

function test_tls_infeasible_targets()
    vI = DistributionsFactories._truncexp_max_var(0.0, 1.0, 0.3)
    checks = [
        () -> truncated_locscale_params(Normal, 0.0, 1.0, 0.3, 1.01 * vI),   # above Γ
        () -> truncated_locscale_params(Normal, 0.0, 1.0, 1.2, 0.01),        # mean outside
        () -> truncated_locscale_params(Normal, 0.0, Inf, 1.0, 1.0),         # σ̄ = μ̄ - a
        () -> truncated_locscale_params(Normal, -Inf, 0.0, -1.0, 1.5),
        () -> truncated_locscale_params(TDist(5.0), 0.0, Inf, 1.0, 5 / 3 + 0.01),
        () -> truncated_locscale_params(TDist(2.0), 0.0, Inf, 1.0, 0.5),     # ν ≤ 2
        () -> truncated_locscale_params(TDist(4.0), 0.0, 1.0, 0.5, 0.26),    # Bhatia–Davis
        () -> truncated_locscale_params(TDist(4.0), 0.0, 1.0, 0.5, 0.2),     # no convergence
        () -> make_dist(truncated(Normal(); lower=0, upper=1), mean=0.3, var=1.01 * vI),
        () -> make_dist(truncated(TDist(4.0); lower=0, upper=1), mean=0.5, var=0.26),
    ]
    all(_throws_infeasible, checks) || return false
    # The error is informative and flags provable infeasibility.
    e = try truncated_locscale_params(Normal, 0.0, 1.0, 0.3, 1.01 * vI) catch e; e end
    e.proven || return false
    occursin("InfeasibleMomentsError", sprint(showerror, e)) || return false
    return true
end

function test_tls_make_dist_two_sided_t()
    d = make_dist(truncated(TDist(4.5); lower=-1, upper=2), mean=0.3, var=0.4)
    m, v = truncated_locscale_moments(TDist(4.5), d.untruncated.μ, d.untruncated.σ, -1, 2)
    return isapprox(m, 0.3; atol=1e-10) && isapprox(v, 0.4; rtol=1e-10)
end
