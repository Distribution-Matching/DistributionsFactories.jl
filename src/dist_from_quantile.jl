# Quantile-based construction: core methods

# --- Single quantile (1-parameter distributions) ---

function dist_from_quantile(::Type{Exponential}, p::Number, q::Number)
    (0 < p < 1) || throw(DomainError(p, "p must be in (0,1)"))
    q > 0 || throw(DomainError(q, "Exponential: quantile must be > 0"))
    θ = -q / log(1 - p)
    return Exponential(θ)
end

# --- Two quantiles (2-parameter location-scale families) ---

function dist_from_quantiles(::Type{Normal}, p1::Number, q1::Number, p2::Number, q2::Number)
    (0 < p1 < 1 && 0 < p2 < 1) || throw(DomainError("p values must be in (0,1)"))
    p1 ≠ p2 || throw(DomainError("p1 and p2 must be distinct"))
    z1 = quantile(Normal(), p1)
    z2 = quantile(Normal(), p2)
    σ = (q2 - q1) / (z2 - z1)
    σ > 0 || throw(DomainError("Quantile specification implies non-positive σ"))
    μ = q1 - σ * z1
    return Normal(μ, σ)
end

function dist_from_quantiles(::Type{Laplace}, p1::Number, q1::Number, p2::Number, q2::Number)
    (0 < p1 < 1 && 0 < p2 < 1) || throw(DomainError("p values must be in (0,1)"))
    p1 ≠ p2 || throw(DomainError("p1 and p2 must be distinct"))
    z1 = quantile(Laplace(), p1)
    z2 = quantile(Laplace(), p2)
    b = (q2 - q1) / (z2 - z1)
    b > 0 || throw(DomainError("Quantile specification implies non-positive scale"))
    μ = q1 - b * z1
    return Laplace(μ, b)
end

function dist_from_quantiles(::Type{Logistic}, p1::Number, q1::Number, p2::Number, q2::Number)
    (0 < p1 < 1 && 0 < p2 < 1) || throw(DomainError("p values must be in (0,1)"))
    p1 ≠ p2 || throw(DomainError("p1 and p2 must be distinct"))
    z1 = quantile(Logistic(), p1)
    z2 = quantile(Logistic(), p2)
    s = (q2 - q1) / (z2 - z1)
    s > 0 || throw(DomainError("Quantile specification implies non-positive scale"))
    μ = q1 - s * z1
    return Logistic(μ, s)
end

function dist_from_quantiles(::Type{Cauchy}, p1::Number, q1::Number, p2::Number, q2::Number)
    (0 < p1 < 1 && 0 < p2 < 1) || throw(DomainError("p values must be in (0,1)"))
    p1 ≠ p2 || throw(DomainError("p1 and p2 must be distinct"))
    z1 = quantile(Cauchy(), p1)
    z2 = quantile(Cauchy(), p2)
    σ = (q2 - q1) / (z2 - z1)
    σ > 0 || throw(DomainError("Quantile specification implies non-positive scale"))
    μ = q1 - σ * z1
    return Cauchy(μ, σ)
end

function dist_from_quantiles(::Type{Gumbel}, p1::Number, q1::Number, p2::Number, q2::Number)
    (0 < p1 < 1 && 0 < p2 < 1) || throw(DomainError("p values must be in (0,1)"))
    p1 ≠ p2 || throw(DomainError("p1 and p2 must be distinct"))
    z1 = quantile(Gumbel(), p1)
    z2 = quantile(Gumbel(), p2)
    β = (q2 - q1) / (z2 - z1)
    β > 0 || throw(DomainError("Quantile specification implies non-positive scale"))
    μ = q1 - β * z1
    return Gumbel(μ, β)
end

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

function dist_from_quantiles(::Type{Beta}, p1::Number, q1::Number, p2::Number, q2::Number)
    (0 < p1 < 1 && 0 < p2 < 1) || throw(DomainError("p values must be in (0,1)"))
    p1 ≠ p2 || throw(DomainError("p1 and p2 must be distinct"))
    (0 < q1 < 1 && 0 < q2 < 1) || throw(DomainError("Beta: quantiles must be in (0,1)"))

    # 2D Newton in log-space: solve for x = [log(α), log(β)]
    # F(x) = [quantile(Beta(exp(x1), exp(x2)), p1) - q1,
    #          quantile(Beta(exp(x1), exp(x2)), p2) - q2]

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

function dist_from_mean_quantile(::Type{Gamma}, μ::Number, p::Number, q::Number)
    μ > 0 || throw(DomainError(μ, "Gamma: μ must be > 0"))
    (0 < p < 1) || throw(DomainError(p, "p must be in (0,1)"))
    q > 0 || throw(DomainError(q, "Gamma: quantile must be > 0"))

    # mean = α*θ, so θ = μ/α. Solve for α.
    α_sol = find_zero(
        α -> quantile(Gamma(α, μ / α), p) - q,
        1.0
    )
    return Gamma(α_sol, μ / α_sol)
end

function dist_from_mean_quantile(::Type{Beta}, μ::Number, p::Number, q::Number)
    (0 < μ < 1) || throw(DomainError(μ, "Beta: μ must be in (0,1)"))
    (0 < p < 1) || throw(DomainError(p, "p must be in (0,1)"))
    (0 < q < 1) || throw(DomainError(q, "Beta: quantile must be in (0,1)"))

    # α = μS, β = (1-μ)S. Solve for S (= α+β) via quantile constraint.
    logS_sol = find_zero(
        logS -> begin
            S = exp(logS)
            quantile(Beta(μ * S, (1 - μ) * S), p) - q
        end,
        1.0
    )
    S = exp(logS_sol)
    return Beta(μ * S, (1 - μ) * S)
end
