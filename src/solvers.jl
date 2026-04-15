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
    _solve_truncated_mean_var(D, lo, hi, target_μ, target_var; x0, maxiter=200, tol=1e-10)

Find location–scale parameters `(μ_parent, σ_parent)` such that
`truncated(D(μ_parent, σ_parent), lo, hi)` has the given mean and variance.

Uses 2D Newton iteration with numerical Jacobian.
"""
function _solve_truncated_mean_var(::Type{D}, lo::Real, hi::Real,
                                   target_μ::Real, target_var::Real;
                                   x0::Vector{Float64} = [target_μ, √target_var],
                                   maxiter::Int = 200, tol::Real = 1e-10) where {D<:ContinuousUnivariateDistribution}
    x = copy(x0)
    h = 1e-7

    function residual(μp, logσ)
        σp = exp(logσ)
        d = truncated(D(μp, σp), lo, hi)
        # Use quadrature for moments since not all truncated distributions
        # define mean()/var() in Distributions.jl
        m, _ = quadgk(x -> x * pdf(d, x), lo, hi)
        m2, _ = quadgk(x -> x^2 * pdf(d, x), lo, hi)
        return [m - target_μ, m2 - m^2 - target_var]
    end

    # Work in (μ_parent, log(σ_parent)) space to keep σ > 0
    x_work = [x[1], log(x[2])]

    for _ in 1:maxiter
        F = residual(x_work[1], x_work[2])

        if maximum(abs.(F)) < tol
            return truncated(D(x_work[1], exp(x_work[2])), lo, hi)
        end

        # Numerical Jacobian
        J = zeros(2, 2)
        for j in 1:2
            xp = copy(x_work); xp[j] += h
            Fp = residual(xp[1], xp[2])
            J[:, j] = (Fp - F) / h
        end

        dx = J \ (-F)

        # Damped step
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

    # Return best result even if not fully converged
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
