"""
    solve_beta_ratio(C)

Find γ > 0 such that `C = γ / B(1/γ, 1/γ)  =  γ · Γ(2/γ) / Γ(1/γ)²`.

The function f(γ) is strictly decreasing on (0,∞) with range (1/2, ∞),
so for any C > 1/2 there is a unique solution.
"""
function solve_beta_ratio(C::Real)
    C > 0.5 || throw(DomainError(C, "C must be > 1/2"))
    isfinite(C) || throw(DomainError(C, "C must be finite"))

    # Special case: f(1) = 1
    isapprox(C, 1.0; atol=1e-12) && return 1.0

    # Root-finding on the log form for numerical stability:
    #   g(γ) = log(γ) + logΓ(2/γ) - 2·logΓ(1/γ) - log(C) = 0
    logC = log(C)
    g(γ) = log(γ) + loggamma(2 / γ) - 2 * loggamma(1 / γ) - logC

    if C > 1  # γ < 1: bracket in [ε, 1]
        lo = 0.5
        while g(lo) < 0
            lo /= 2
        end
        hi = 1.0
    else      # C < 1: γ > 1: bracket in [1, M]
        lo = 1.0
        hi = 2.0
        while g(hi) > 0
            hi *= 2
        end
    end

    return find_zero(g, (lo, hi))
end

"""
    _solve_truncated_unit(D, μ̄_std, σ̄²_std; maxiter=200, tol=1e-10)

Solve the canonical 2D problem on the unit interval `[-0.5, 0.5]`:
find `(μp_std, σp_std)` such that `truncated(D(μp_std, σp_std), -0.5, 0.5)`
has mean `μ̄_std` and variance `σ̄²_std`.

Throws on non-convergence rather than silently returning a partial solution.
"""
function _solve_truncated_unit(::Type{D}, μ̄_std::Real, σ̄²_std::Real;
                               maxiter::Int = 200, tol::Real = 1e-10
                               ) where {D<:ContinuousUnivariateDistribution}
    -0.5 < μ̄_std < 0.5 || throw(DomainError(μ̄_std,
        "standardized mean must be in (-0.5, 0.5)"))
    σ̄²_std > 0 || throw(DomainError(σ̄²_std,
        "standardized variance must be > 0"))

    # Initial guess: parent μ near the standardized target mean,
    # parent σ near √(target var). Both are O(1) on the unit interval.
    x = [Float64(μ̄_std), log(√Float64(σ̄²_std))]
    h = 1e-7
    a, b = -0.5, 0.5

    function residual(μp, logσ)
        σp = exp(logσ)
        d = truncated(D(μp, σp), a, b)
        m, _  = quadgk(x -> x * pdf(d, x), a, b; rtol=1e-10)
        m2, _ = quadgk(x -> x^2 * pdf(d, x), a, b; rtol=1e-10)
        return [m - μ̄_std, m2 - m^2 - σ̄²_std]
    end

    converged = false
    for _ in 1:maxiter
        F = residual(x[1], x[2])
        if maximum(abs.(F)) < tol
            converged = true
            break
        end

        J = zeros(2, 2)
        for j in 1:2
            xp = copy(x); xp[j] += h
            Fp = residual(xp[1], xp[2])
            J[:, j] = (Fp - F) / h
        end
        dx = J \ (-F)

        step = 1.0
        for _ in 1:20
            x_new = x + step * dx
            try
                F_new = residual(x_new[1], x_new[2])
                if maximum(abs.(F_new)) < maximum(abs.(F))
                    break
                end
            catch
            end
            step *= 0.5
        end
        x .+= step * dx
    end

    converged || throw(ErrorException(
        "_solve_truncated_unit: did not converge for D=$D after $maxiter " *
        "iterations (target μ̄=$μ̄_std, σ̄²=$σ̄²_std)"))

    return x[1], exp(x[2])
end

# Location-scale families: Y = (X - c)/w preserves family membership, i.e.
# truncated(D(μp, σp), a, b) ≡ truncated(D((μp-c)/w, σp/w), -0.5, 0.5) under the
# affine map. Add new families here as needed.
const _TruncLocScale = Union{Normal, Laplace, Logistic}

"""
    _solve_truncated_half_below_unit(D, target_z; maxiter=200, tol=1e-10)

Canonical half-truncated solver. Find `(μp, σp)` such that
`truncated(D(μp, σp), 0, Inf)` has mean `target_z` and variance 1.
Used as the inner kernel by the lower- and upper-half-truncated factories
after standardization. `target_z` must be > 1 (the exponential-tail bound).
"""
function _solve_truncated_half_below_unit(::Type{D}, target_z::Real;
                                          maxiter::Int = 200, tol::Real = 1e-10
                                          ) where {D<:_TruncLocScale}
    target_z > 1 || throw(DomainError(target_z,
        "standardized mean z = (μ̄-lo)/σ̄ must be > 1 (exponential-tail bound)"))

    # Initial guess: parent at the truncation point with std 1 (≈ half-D).
    x = [0.0, 0.0]   # (μp, log σp)
    h = 1e-7

    function residual(μp, logσ)
        σp = exp(logσ)
        d = truncated(D(μp, σp), 0.0, Inf)
        # Use a finite upper cap matched to (μp, σp) for quadrature stability.
        # 50σ above max(0, μp) captures all but ~10⁻³⁰⁰ of an exponential tail.
        upper = max(0.0, μp) + 50 * σp
        m, _  = quadgk(x -> x   * pdf(d, x), 0.0, upper; rtol=1e-10)
        m2, _ = quadgk(x -> x^2 * pdf(d, x), 0.0, upper; rtol=1e-10)
        return [m - target_z, m2 - m^2 - 1.0]
    end

    converged = false
    for _ in 1:maxiter
        F = residual(x[1], x[2])
        if maximum(abs.(F)) < tol
            converged = true
            break
        end

        J = zeros(2, 2)
        for j in 1:2
            xp = copy(x); xp[j] += h
            Fp = residual(xp[1], xp[2])
            J[:, j] = (Fp - F) / h
        end
        dx = J \ (-F)

        step = 1.0
        for _ in 1:20
            x_new = x + step * dx
            try
                F_new = residual(x_new[1], x_new[2])
                if maximum(abs.(F_new)) < maximum(abs.(F))
                    break
                end
            catch
            end
            step *= 0.5
        end
        x .+= step * dx
    end

    converged || throw(ErrorException(
        "_solve_truncated_half_below_unit: did not converge for D=$D after " *
        "$maxiter iterations (target z=$target_z)"))

    return x[1], exp(x[2])
end

"""
    _solve_truncated_half_below(D, lo, target_μ, target_var)

Solve the half-truncated problem on `[lo, ∞)`: find `(μp, σp)` such that
`truncated(D(μp, σp), lo, Inf)` has mean `target_μ` and variance `target_var`.
Standardizes to the canonical interval `[0, ∞)` with `(μ̄ - lo)/σ̄, 1`,
calls [`_solve_truncated_half_below_unit`](@ref), then un-standardizes.
"""
function _solve_truncated_half_below(::Type{D}, lo::Real,
                                     target_μ::Real, target_var::Real;
                                     maxiter::Int = 200, tol::Real = 1e-10
                                     ) where {D<:_TruncLocScale}
    σ̄ = √target_var
    z = (target_μ - lo) / σ̄
    # Boundary attainment for Laplace: σ̄ = (μ̄ - lo) gives exact Exponential.
    # Pick the canonical boundary parent Laplace(lo, σ̄) directly.
    if D <: Laplace && z ≤ 1 + 1e-8
        return truncated(Laplace(lo, σ̄), lo, Inf)
    end
    μp_std, σp_std = _solve_truncated_half_below_unit(D, z;
                                                      maxiter=maxiter, tol=tol)
    μp = lo + σ̄ * μp_std
    σp = σ̄ * σp_std
    return truncated(D(μp, σp), lo, Inf)
end

"""
    _solve_truncated_half_above(D, hi, target_μ, target_var)

Solve `(-∞, hi]` by reflection: the right-flipped problem is half-below at
`-hi` with mean `-target_μ`. Returns the reflected truncated distribution.
"""
function _solve_truncated_half_above(::Type{D}, hi::Real,
                                     target_μ::Real, target_var::Real;
                                     maxiter::Int = 200, tol::Real = 1e-10
                                     ) where {D<:_TruncLocScale}
    # Reflect: if Y = -X then truncated(D(μp, σp), -∞, hi) on X corresponds to
    # truncated(D(-μp, σp), -hi, ∞) on Y. (Normal/Laplace/Logistic are all
    # symmetric, so reflecting the parent location works.)
    d_reflected = _solve_truncated_half_below(D, -hi, -target_μ, target_var;
                                              maxiter=maxiter, tol=tol)
    μp_ref, σp_ref = params(d_reflected.untruncated)
    return truncated(D(-μp_ref, σp_ref), -Inf, hi)
end

"""
    _solve_truncated_mean_var(D, lo, hi, target_μ, target_var; maxiter=200, tol=1e-10)

Standardize-and-solve for location-scale families (Normal, Laplace, Logistic).
Maps `[lo, hi] → [-0.5, 0.5]` via `Y = (X - c)/w` where `c = (lo+hi)/2`,
`w = hi - lo`, calls [`_solve_truncated_unit`](@ref) on the canonical
interval, then un-standardizes the recovered parent parameters to the
original `[lo, hi]`. Numerically much more stable than solving in user
coordinates because all moments and parameters live in O(1) range.
"""
function _solve_truncated_mean_var(::Type{D}, lo::Real, hi::Real,
                                   target_μ::Real, target_var::Real;
                                   maxiter::Int = 200, tol::Real = 1e-10
                                   ) where {D<:_TruncLocScale}
    c = (lo + hi) / 2
    w = hi - lo
    μp_std, σp_std = _solve_truncated_unit(D,
                                           (target_μ - c) / w,
                                           target_var / w^2;
                                           maxiter=maxiter, tol=tol)
    μp = c + w * μp_std
    σp = w * σp_std
    return truncated(D(μp, σp), lo, hi)
end

"""
    _solve_truncated_mean_var(D, lo, hi, target_μ, target_var; …)

Generic fallback for non-location-scale families (e.g. `Gamma`). Solves the
2D Newton problem directly in user coordinates because there is no closed-form
affine map that preserves the family on the canonical interval.
"""
function _solve_truncated_mean_var(::Type{D}, lo::Real, hi::Real,
                                   target_μ::Real, target_var::Real;
                                   x0::Vector{Float64} = [target_μ, √target_var],
                                   maxiter::Int = 200, tol::Real = 1e-10
                                   ) where {D<:ContinuousUnivariateDistribution}
    x_work = [Float64(x0[1]), log(Float64(x0[2]))]
    h = 1e-7

    function residual(μp, logσ)
        σp = exp(logσ)
        d = truncated(D(μp, σp), lo, hi)
        m, _  = quadgk(x -> x * pdf(d, x), lo, hi)
        m2, _ = quadgk(x -> x^2 * pdf(d, x), lo, hi)
        return [m - target_μ, m2 - m^2 - target_var]
    end

    for _ in 1:maxiter
        F = residual(x_work[1], x_work[2])
        if maximum(abs.(F)) < tol
            return truncated(D(x_work[1], exp(x_work[2])), lo, hi)
        end

        J = zeros(2, 2)
        for j in 1:2
            xp = copy(x_work); xp[j] += h
            Fp = residual(xp[1], xp[2])
            J[:, j] = (Fp - F) / h
        end
        dx = J \ (-F)

        step = 1.0
        for _ in 1:20
            x_new = x_work + step * dx
            try
                F_new = residual(x_new[1], x_new[2])
                if maximum(abs.(F_new)) < maximum(abs.(F))
                    break
                end
            catch
            end
            step *= 0.5
        end
        x_work .+= step * dx
    end

    # Generic path: keep the historical "return best" behaviour for now.
    return truncated(D(x_work[1], exp(x_work[2])), lo, hi)
end

"""
    _folded_normal_moments(μ_parent, σ_parent)

Compute the mean and variance of the folded normal |X| where X ~ Normal(μ_parent, σ_parent).
"""
function _folded_normal_moments(μp::Real, σp::Real)
    μ_fn = σp * √(2/π) * exp(-μp^2 / (2σp^2)) + μp * erf(μp / (σp * √2))
    var_fn = μp^2 + σp^2 - μ_fn^2
    return μ_fn, var_fn
end

"""
    _truncated_poisson_moments(λ, lo, hi) -> (mean, var)

Direct summation. `Distributions.jl` does not currently define `mean`/`var` for
`Truncated{<:Poisson}`; we sum the (renormalised) PMF over the integer support
in `[ceil(lo), floor(hi)]`.
"""
function _truncated_poisson_moments(λ::Real, lo::Real, hi::Real)
    p = Poisson(λ)
    lo_i = ceil(Int, lo)
    hi_i = floor(Int, hi)
    Z = cdf(p, hi_i) - (lo_i > 0 ? cdf(p, lo_i - 1) : 0.0)
    Z > 0 || return (NaN, NaN)
    m = 0.0
    m2 = 0.0
    for k in lo_i:hi_i
        pk = pdf(p, k) / Z
        m  += k * pk
        m2 += k^2 * pk
    end
    return (m, m2 - m^2)
end

"""
    _solve_truncated_poisson_mean(lo, hi, target_μ)

Find Poisson rate `λ > 0` such that `truncated(Poisson(λ); lower=lo, upper=hi)`
has mean `target_μ`. The truncated mean is monotone increasing in λ on the open
interval `(lo, hi)`, so a 1D root-find suffices.
"""
function _solve_truncated_poisson_mean(lo::Real, hi::Real, target_μ::Real)
    (lo < target_μ < hi) || throw(DomainError(target_μ, "target_μ must lie in ($lo, $hi)"))

    truncated_mean(λ) = _truncated_poisson_moments(λ, lo, hi)[1]

    lo_λ = max(target_μ / 4, 1e-3)
    hi_λ = max(target_μ * 4, lo_λ + 1.0)
    while truncated_mean(lo_λ) > target_μ && lo_λ > 1e-8
        lo_λ /= 2
    end
    while truncated_mean(hi_λ) < target_μ
        hi_λ *= 2
    end

    λ = find_zero(λ -> truncated_mean(λ) - target_μ, (lo_λ, hi_λ))
    return truncated(Poisson(λ); lower=lo, upper=hi)
end

"""
    _solve_folded_normal(target_μ, target_var; maxiter=200, tol=1e-10)

Find `(μ_parent, σ_parent)` such that `|Normal(μ_parent, σ_parent)|` has the given mean and variance.
Uses 2D Newton iteration.
"""
function _solve_folded_normal(target_μ::Real, target_var::Real;
                              maxiter::Int = 200, tol::Real = 1e-10)
    # Initial guess: if target_μ >> √target_var, the fold barely matters, so μp ≈ target_μ
    x = [target_μ, log(√target_var)]  # (μ_parent, log(σ_parent))
    h = 1e-7

    function residual(μp, logσ)
        σp = exp(logσ)
        m, v = _folded_normal_moments(μp, σp)
        return [m - target_μ, v - target_var]
    end

    for _ in 1:maxiter
        F = residual(x[1], x[2])

        if maximum(abs.(F)) < tol
            return x[1], exp(x[2])
        end

        J = zeros(2, 2)
        for j in 1:2
            xp = copy(x); xp[j] += h
            Fp = residual(xp[1], xp[2])
            J[:, j] = (Fp - F) / h
        end

        dx = J \ (-F)

        step = 1.0
        for _ in 1:20
            x_new = x + step * dx
            try
                F_new = residual(x_new[1], x_new[2])
                if maximum(abs.(F_new)) < maximum(abs.(F))
                    break
                end
            catch
            end
            step *= 0.5
        end
        x .+= step * dx
    end

    return x[1], exp(x[2])
end
