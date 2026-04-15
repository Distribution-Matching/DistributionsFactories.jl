# --- Affine: bounded [0,1] → [a,b] ---

function test_support_beta_scaled()
    for (a, b) in [(-3.0, 2.0), (2.0, 7.0), (0.0, 10.0)]
        for α in [1.0, 2.0, 5.0]
            for β in [1.0, 3.0]
                d_std = Beta(α, β)
                μ_std, v_std = mean(d_std), var(d_std)
                w = b - a
                μ̄ = a + w * μ_std
                σ̄² = w^2 * v_std
                d = dist_from_mean_var_on_support(Beta, μ̄, σ̄²; support=a..b)
                if !isapprox(mean(d), μ̄, rtol=1e-8) || !isapprox(var(d), σ̄², rtol=1e-8)
                    @info "Beta scaled mismatch" a b α β mean(d) μ̄ var(d) σ̄²
                    return false
                end
                if !isapprox(minimum(d), a, atol=1e-10) || !isapprox(maximum(d), b, atol=1e-10)
                    @info "Beta scaled bounds mismatch" minimum(d) a maximum(d) b
                    return false
                end
            end
        end
    end
    return true
end

function test_support_uniform_scaled()
    d = dist_from_mean_var_on_support(Uniform, 0.0, 8.333333333333334; support=-5..5)
    isapprox(mean(d), 0.0, atol=1e-8) || return false
    isapprox(minimum(d), -5.0, atol=1e-10) || return false
    isapprox(maximum(d), 5.0, atol=1e-10) || return false
    return true
end

# --- Affine: positive [0,∞) → [a,∞) shift ---

function test_support_gamma_shifted()
    for a in [1.0, 3.0, -2.0]
        for α in [2.0, 5.0]
            for θ in [1.0, 2.0]
                d_std = Gamma(α, θ)
                μ̄ = a + mean(d_std)
                σ̄² = var(d_std)
                d = dist_from_mean_var_on_support(Gamma, μ̄, σ̄²; support=a..Inf)
                if !isapprox(mean(d), μ̄, rtol=1e-8) || !isapprox(var(d), σ̄², rtol=1e-8)
                    @info "Gamma shift mismatch" a α θ mean(d) μ̄
                    return false
                end
                if !isapprox(minimum(d), a, atol=1e-10)
                    @info "Gamma shift bound mismatch" minimum(d) a
                    return false
                end
            end
        end
    end
    return true
end

# --- Affine: positive [0,∞) → (-∞,b] flip ---

function test_support_gamma_flipped()
    for b in [5.0, 10.0, 0.0]
        d_std = Gamma(3.0, 1.0)
        μ̄ = b - mean(d_std)
        σ̄² = var(d_std)
        d = dist_from_mean_var_on_support(Gamma, μ̄, σ̄²; support=-Inf..b)
        if !isapprox(mean(d), μ̄, rtol=1e-8) || !isapprox(var(d), σ̄², rtol=1e-8)
            @info "Gamma flip mismatch" b mean(d) μ̄
            return false
        end
        if !isapprox(maximum(d), b, atol=1e-10)
            @info "Gamma flip bound mismatch" maximum(d) b
            return false
        end
    end
    return true
end

# --- Discrete: bounded {0,...,n} → {a,...,b} ---

function test_support_binomial_shifted()
    for a in [5, 10, -3]
        n = 8
        b = a + n
        p = 0.4
        d_std = Binomial(n, p)
        μ̄ = a + mean(d_std)
        σ̄² = var(d_std)
        d = dist_from_mean_var_on_support(Binomial, μ̄, σ̄²; support=a:b)
        if !isapprox(mean(d), μ̄, atol=1e-8) || !isapprox(var(d), σ̄², rtol=1e-8)
            @info "Binomial shift mismatch" a mean(d) μ̄
            return false
        end
        if minimum(d) != a || maximum(d) != b
            @info "Binomial shift bounds mismatch" minimum(d) a maximum(d) b
            return false
        end
    end
    return true
end

function test_support_discrete_uniform_shifted()
    d = dist_from_mean_var_on_support(DiscreteUniform, 12.5, 2.9166666666666665; support=10:15)
    isapprox(mean(d), 12.5, atol=1e-8) || return false
    minimum(d) == 10 || return false
    maximum(d) == 15 || return false
    return true
end

# --- Truncation: real-line → bounded ---

function test_support_normal_truncated()
    d = dist_from_mean_var_on_support(Normal, 0.5, 0.04; support=0..1)
    isapprox(mean(d), 0.5, rtol=1e-4) || return false
    isapprox(var(d), 0.04, rtol=1e-4) || return false
    isapprox(minimum(d), 0.0, atol=1e-10) || return false
    isapprox(maximum(d), 1.0, atol=1e-10) || return false
    return true
end

# --- Truncation: positive → bounded ---

function test_support_gamma_truncated()
    d = dist_from_mean_var_on_support(Gamma, 3.0, 1.0; support=0..10)
    # Use quadrature since mean()/var() may not be defined for Truncated{Gamma}
    m, _ = DistributionsFactories.quadgk(x -> x * pdf(d, x), 0, 10)
    m2, _ = DistributionsFactories.quadgk(x -> x^2 * pdf(d, x), 0, 10)
    v = m2 - m^2
    isapprox(m, 3.0, rtol=1e-4) || return false
    isapprox(v, 1.0, rtol=1e-4) || return false
    isapprox(minimum(d), 0.0, atol=1e-10) || return false
    isapprox(maximum(d), 10.0, atol=1e-10) || return false
    return true
end

# --- Error cases ---

function test_support_errors()
    # Unit distribution on unbounded interval
    try
        dist_from_mean_var_on_support(Beta, 0.5, 0.05; support=0..Inf)
        return false
    catch e
        e isa ArgumentError || return false
    end

    # Positive distribution on (-∞,∞)
    try
        dist_from_mean_var_on_support(Gamma, 5.0, 3.0; support=-Inf..Inf)
        return false
    catch e
        e isa ArgumentError || return false
    end

    # Positive distribution on [a,b] with a < 0
    try
        dist_from_mean_var_on_support(Gamma, 5.0, 3.0; support=-1..10)
        return false
    catch e
        e isa ArgumentError || return false
    end

    return true
end
