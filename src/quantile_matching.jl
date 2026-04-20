# Quantile-based construction: core methods

# --- Helper for location-scale families ---

function _dist_from_quantiles_location_scale(::Type{D}, p1::Number, q1::Number, p2::Number, q2::Number) where {D<:ContinuousUnivariateDistribution}
    (0 < p1 < 1 && 0 < p2 < 1) || throw(DomainError("p values must be in (0,1)"))
    p1 ≠ p2 || throw(DomainError("p1 and p2 must be distinct"))
    z1 = quantile(D(), p1)
    z2 = quantile(D(), p2)
    s = (q2 - q1) / (z2 - z1)
    s > 0 || throw(DomainError("Quantile specification implies non-positive scale"))
    μ = q1 - s * z1
    return D(μ, s)
end

# --- Single quantile (1-parameter distributions) ---

"""
    dist_from_quantile(D, p, q)

Construct a 1-parameter distribution `D` such that its `p`-th quantile equals `q`.

Supported distributions: `Exponential`.

See also: [`dist_from_quantiles`](@ref), [`make_dist`](@ref)
"""
function dist_from_quantile end

function dist_from_quantile(::Type{Exponential}, p::Number, q::Number)
    (0 < p < 1) || throw(DomainError(p, "p must be in (0,1)"))
    q > 0 || throw(DomainError(q, "Exponential: quantile must be > 0"))
    θ = -q / log(1 - p)
    return Exponential(θ)
end

function dist_from_quantile(::Type{Rayleigh}, p::Number, q::Number)
    (0 < p < 1) || throw(DomainError(p, "p must be in (0,1)"))
    q > 0 || throw(DomainError(q, "Rayleigh: quantile must be > 0"))
    # Quantile of Rayleigh(σ): q = σ √(-2 log(1-p))
    σ = q / √(-2 * log(1 - p))
    return Rayleigh(σ)
end

function dist_from_quantile(::Type{Geometric}, p::Number, q::Number)
    (0 < p < 1) || throw(DomainError(p, "p must be in (0,1)"))
    q ≥ 0 && isinteger(q) || throw(DomainError(q, "Geometric: quantile must be a non-negative integer"))
    # Geometric(prob) on {0,1,…} has CDF P(X ≤ k) = 1 - (1-prob)^(k+1).
    # The exact p-th quantile is the smallest k with 1 - (1-prob)^(k+1) ≥ p,
    # i.e. (1-prob)^(k+1) ≤ 1-p. For continuous matching we set equality at
    # k+1, then `prob = 1 - (1-p)^(1/(q+1))`.
    prob = 1 - (1 - p)^(1 / (q + 1))
    return Geometric(prob)
end

# --- Two quantiles (2-parameter distributions) ---

"""
    dist_from_quantiles(D, p1, q1, p2, q2)

Construct a 2-parameter distribution `D` such that its `p1`-th quantile equals `q1`
and its `p2`-th quantile equals `q2`.

Direct formula for location-scale families (`Normal`, `Laplace`, `Logistic`, `Cauchy`, `Gumbel`).
Numerical (root-finding) for shape-scale families (`Gamma`, `Beta`).

See also: [`make_dist`](@ref)
"""
function dist_from_quantiles end

dist_from_quantiles(::Type{Normal}, p1::Number, q1::Number, p2::Number, q2::Number) =
    _dist_from_quantiles_location_scale(Normal, p1, q1, p2, q2)

dist_from_quantiles(::Type{Laplace}, p1::Number, q1::Number, p2::Number, q2::Number) =
    _dist_from_quantiles_location_scale(Laplace, p1, q1, p2, q2)

dist_from_quantiles(::Type{Logistic}, p1::Number, q1::Number, p2::Number, q2::Number) =
    _dist_from_quantiles_location_scale(Logistic, p1, q1, p2, q2)

dist_from_quantiles(::Type{Cauchy}, p1::Number, q1::Number, p2::Number, q2::Number) =
    _dist_from_quantiles_location_scale(Cauchy, p1, q1, p2, q2)

dist_from_quantiles(::Type{Gumbel}, p1::Number, q1::Number, p2::Number, q2::Number) =
    _dist_from_quantiles_location_scale(Gumbel, p1, q1, p2, q2)

# --- Two quantiles (2-parameter non-location-scale, numerical) ---

function dist_from_quantiles(::Type{Gamma}, p1::Number, q1::Number, p2::Number, q2::Number)
    (0 < p1 < 1 && 0 < p2 < 1) || throw(DomainError("p values must be in (0,1)"))
    p1 ≠ p2 || throw(DomainError("p1 and p2 must be distinct"))
    q1 > 0 && q2 > 0 || throw(DomainError("Gamma: quantiles must be > 0"))

    # Gamma is a scale family: for any α, set θ = q1 / quantile(Gamma(α,1), p1).
    # Work in log-space to keep α > 0.
    r = q2 / q1
    f(logα) = begin
        α = exp(logα)
        quantile(Gamma(α, 1.0), p2) / quantile(Gamma(α, 1.0), p1) - r
    end

    logα_sol = find_zero(f, 0.0)
    α_sol = exp(logα_sol)
    θ = q1 / quantile(Gamma(α_sol, 1.0), p1)
    return Gamma(α_sol, θ)
end

function dist_from_quantiles(::Type{LogNormal}, p1::Number, q1::Number, p2::Number, q2::Number)
    (0 < p1 < 1 && 0 < p2 < 1) || throw(DomainError("p values must be in (0,1)"))
    p1 ≠ p2 || throw(DomainError("p1 and p2 must be distinct"))
    q1 > 0 && q2 > 0 || throw(DomainError("LogNormal: quantiles must be > 0"))
    # log(X) ~ Normal(μ_log, σ_log); reduce to the location-scale Normal solver
    # on the log-quantiles.
    d_log = _dist_from_quantiles_location_scale(Normal, p1, log(q1), p2, log(q2))
    return LogNormal(d_log.μ, d_log.σ)
end

function dist_from_quantiles(::Type{Weibull}, p1::Number, q1::Number, p2::Number, q2::Number)
    (0 < p1 < 1 && 0 < p2 < 1) || throw(DomainError("p values must be in (0,1)"))
    p1 ≠ p2 || throw(DomainError("p1 and p2 must be distinct"))
    q1 > 0 && q2 > 0 || throw(DomainError("Weibull: quantiles must be > 0"))
    # Weibull(k, λ) has quantile q = λ·(-log(1-p))^(1/k). Eliminate λ:
    # q2/q1 = (log(1-p2)/log(1-p1))^(1/k). Solve for k in log-space.
    r = q2 / q1
    L = log(1 - p2) / log(1 - p1)
    L > 0 && L != 1 || throw(DomainError("Weibull: degenerate p values"))
    k = log(L) / log(r)
    k > 0 || throw(DomainError("Weibull: quantile spec implies non-positive shape"))
    λ = q1 / (-log(1 - p1))^(1 / k)
    return Weibull(k, λ)
end

function dist_from_quantiles(::Type{Pareto}, p1::Number, q1::Number, p2::Number, q2::Number)
    (0 < p1 < 1 && 0 < p2 < 1) || throw(DomainError("p values must be in (0,1)"))
    p1 ≠ p2 || throw(DomainError("p1 and p2 must be distinct"))
    q1 > 0 && q2 > 0 || throw(DomainError("Pareto: quantiles must be > 0"))
    # Pareto(α, θ) has quantile q = θ·(1-p)^(-1/α). Same elimination as Weibull.
    r = q2 / q1
    L = (1 - p2) / (1 - p1)   # < 1 if p2 > p1
    L != 1 || throw(DomainError("Pareto: degenerate p values"))
    α = -log(L) / log(r)
    α > 0 || throw(DomainError("Pareto: quantile spec implies non-positive shape"))
    θ = q1 * (1 - p1)^(1 / α)
    return Pareto(α, θ)
end

function dist_from_quantiles(::Type{Beta}, p1::Number, q1::Number, p2::Number, q2::Number)
    (0 < p1 < 1 && 0 < p2 < 1) || throw(DomainError("p values must be in (0,1)"))
    p1 ≠ p2 || throw(DomainError("p1 and p2 must be distinct"))
    (0 < q1 < 1 && 0 < q2 < 1) || throw(DomainError("Beta: quantiles must be in (0,1)"))

    # 2D Newton in log-space: solve for x = [log(α), log(β)]
    # Initial guess from normal approximation
    z1 = quantile(Normal(), p1)
    z2 = quantile(Normal(), p2)
    σ_est = max((q2 - q1) / (z2 - z1), 1e-4)
    μ_est = clamp(q1 - σ_est * z1, 0.01, 0.99)
    v_est = clamp(σ_est^2, 1e-6, μ_est * (1 - μ_est) * 0.99)
    S0 = μ_est * (1 - μ_est) / v_est - 1
    x = [log(max(μ_est * S0, 0.1)), log(max((1 - μ_est) * S0, 0.1))]

    h = 1e-7
    for _ in 1:200
        α, β = exp(x[1]), exp(x[2])
        d = Beta(α, β)
        F = [quantile(d, p1) - q1, quantile(d, p2) - q2]

        if maximum(abs.(F)) < 1e-12
            return Beta(α, β)
        end

        # Numerical Jacobian
        J = zeros(2, 2)
        for j in 1:2
            xp = copy(x); xp[j] += h
            αp, βp = exp(xp[1]), exp(xp[2])
            dp = Beta(αp, βp)
            Fp = [quantile(dp, p1) - q1, quantile(dp, p2) - q2]
            J[:, j] = (Fp - F) / h
        end

        dx = J \ (-F)
        # Damped step to stay in valid region
        step = 1.0
        for _ in 1:20
            x_new = x + step * dx
            try
                α_new, β_new = exp(x_new[1]), exp(x_new[2])
                F_new = [quantile(Beta(α_new, β_new), p1) - q1,
                         quantile(Beta(α_new, β_new), p2) - q2]
                if maximum(abs.(F_new)) < maximum(abs.(F))
                    break
                end
            catch
            end
            step *= 0.5
        end
        x .+= step * dx
    end

    α, β = exp(x[1]), exp(x[2])
    return Beta(α, β)
end

# --- Hybrid: mean + quantile ---

"""
    dist_from_mean_quantile(D, μ̄, p, q)

Numerical (root-finding). Construct a 2-parameter distribution `D` with mean `μ̄`
and `p`-th quantile equal to `q`.

Supported distributions: `Gamma`, `Beta`.

See also: [`make_dist`](@ref)
"""
function dist_from_mean_quantile end

function dist_from_mean_quantile(::Type{Gamma}, μ̄::Number, p::Number, q::Number)
    μ̄ > 0 || throw(DomainError(μ̄, "Gamma: μ̄ must be > 0"))
    (0 < p < 1) || throw(DomainError(p, "p must be in (0,1)"))
    q > 0 || throw(DomainError(q, "Gamma: quantile must be > 0"))

    # mean = α*θ, so θ = μ̄/α. Solve for α.
    α_sol = find_zero(
        α -> quantile(Gamma(α, μ̄ / α), p) - q,
        1.0
    )
    return Gamma(α_sol, μ̄ / α_sol)
end

function dist_from_mean_quantile(::Type{LogNormal}, μ̄::Number, p::Number, q::Number)
    μ̄ > 0 || throw(DomainError(μ̄, "LogNormal: μ̄ must be > 0"))
    (0 < p < 1) || throw(DomainError(p, "p must be in (0,1)"))
    q > 0 || throw(DomainError(q, "LogNormal: quantile must be > 0"))
    # mean = exp(μ + σ²/2), so μ = log(μ̄) - σ²/2. quantile constraint:
    # exp(μ + σ·z_p) = q where z_p = quantile(Normal(), p). Substitute.
    # Solve for σ via 1D root-find.
    z_p = quantile(Normal(), p)
    σ = find_zero(σ -> log(μ̄) - σ^2/2 + σ*z_p - log(q), 0.5)
    return LogNormal(log(μ̄) - σ^2/2, σ)
end

function dist_from_mean_quantile(::Type{Beta}, μ̄::Number, p::Number, q::Number)
    (0 < μ̄ < 1) || throw(DomainError(μ̄, "Beta: μ̄ must be in (0,1)"))
    (0 < p < 1) || throw(DomainError(p, "p must be in (0,1)"))
    (0 < q < 1) || throw(DomainError(q, "Beta: quantile must be in (0,1)"))

    # α = μ̄S, β = (1-μ̄)S. Solve for S (= α+β) via quantile constraint.
    logS_sol = find_zero(
        logS -> begin
            S = exp(logS)
            quantile(Beta(μ̄ * S, (1 - μ̄) * S), p) - q
        end,
        1.0
    )
    S = exp(logS_sol)
    return Beta(μ̄ * S, (1 - μ̄) * S)
end
