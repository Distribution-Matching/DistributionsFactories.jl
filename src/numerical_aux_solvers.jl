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
