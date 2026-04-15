function test_mean_var_Beta()
    for α ∈ 0.5:0.5:10
        for β ∈ 0.5:0.5:10
            true_dist = Beta(α, β)
            m, v = mean(true_dist), var(true_dist)
            new_dist = dist_from_mean_var(Beta, m, v)
            if !all(isapprox.(params(true_dist), params(new_dist), atol = 1e-8))
                @info "Mismatch:", new_dist, true_dist
                return false
            end
        end
    end
    return true
end

function test_mean_var_Chisq()
    for k ∈ 1:50
        true_dist = Chisq(k)
        m, v = mean(true_dist), var(true_dist)
        new_dist = dist_from_mean_var(Chisq, m, v)
        if !all(isapprox.(params(true_dist), params(new_dist), atol = 1e-8))
            @info "Mismatch:", new_dist, true_dist
            return false
        end
    end
    return true
end

function test_mean_var_Erlang()
    for k ∈ 1:10
        for θ ∈ 0.5:0.5:5
            true_dist = Erlang(k, θ)
            m, v = mean(true_dist), var(true_dist)
            new_dist = dist_from_mean_var(Erlang, m, v)
            if !all(isapprox.(params(true_dist), params(new_dist), atol = 1e-8))
                @info "Mismatch:", new_dist, true_dist
                return false
            end
        end
    end
    return true
end

function test_mean_var_Exponential()
    for θ ∈ 0.1:0.1:10
        true_dist = Exponential(θ)
        m, v = mean(true_dist), var(true_dist)
        new_dist = dist_from_mean_var(Exponential, m, v)
        if !all(isapprox.(params(true_dist), params(new_dist), atol = 1e-8))
            @info "Mismatch:", new_dist, true_dist
            return false
        end
    end
    return true
end

function test_mean_var_Weibull_simple_bracketing()
    for k ∈ 0.5:0.5:10
        for λ ∈ 0.5:0.5:10
            true_dist = Weibull(k, λ)
            m, v = mean(true_dist), var(true_dist)
            new_dist = dist_from_mean_var(Weibull, m, v; method=:simple_bracketing)
            if !all(isapprox.(params(true_dist), params(new_dist), atol = 1e-8))
                @info "Mismatch:", new_dist, true_dist
                return false
            end
        end
    end
    return true
end

function test_mean_var_Weibull_oscar_garcia()
    for k ∈ 0.5:0.5:10
        for λ ∈ 0.5:0.5:10
            true_dist = Weibull(k, λ)
            m, v = mean(true_dist), var(true_dist)
            new_dist = dist_from_mean_var(Weibull, m, v; method=:oscar_garcia)
            if !all(isapprox.(params(true_dist), params(new_dist), atol = 1e-8))
                @info "Mismatch:", new_dist, true_dist
                return false
            end
        end
    end
    return true
end

function test_mean_var_Weibull_methods_agree()
    for k ∈ 0.5:0.5:10
        for λ ∈ 0.5:0.5:10
            true_dist = Weibull(k, λ)
            m, v = mean(true_dist), var(true_dist)
            d1 = dist_from_mean_var(Weibull, m, v; method=:simple_bracketing)
            d2 = dist_from_mean_var(Weibull, m, v; method=:oscar_garcia)
            if !all(isapprox.(params(d1), params(d2), atol = 1e-8))
                @info "Methods disagree:", d1, d2
                return false
            end
        end
    end
    return true
end

function test_mean_var_TDist_instance()
    for ν ∈ [3, 5, 10, 30, 100]
        for target_μ ∈ [-10.0, 0.0, 5.0, 42.0]
            for target_var ∈ [0.1, 1.0, 4.0, 25.0]
                d = dist_from_mean_var(TDist(ν), target_μ, target_var)
                if !isapprox(mean(d), target_μ, atol=1e-8)
                    @info "Mean mismatch:", ν, target_μ, target_var, mean(d)
                    return false
                end
                if !isapprox(var(d), target_var, rtol=1e-8)
                    @info "Var mismatch:", ν, target_μ, target_var, var(d)
                    return false
                end
            end
        end
    end
    return true
end

function test_mean_var_TDist_instance_returns_affine()
    d = dist_from_mean_var(TDist(5), 10.0, 4.0)
    return d isa LocationScale
end

function test_mean_var_TDist_instance_low_dof_errors()
    try
        dist_from_mean_var(TDist(2), 0.0, 1.0)
        return false
    catch e
        return e isa DomainError
    end
end

function test_mean_var_FDist()
    for d₁ ∈ 2:0.5:50
        for d₂ ∈ 4.5:0.5:50
            true_dist = FDist(d₁, d₂)
            m, v = mean(true_dist), var(true_dist)
            new_dist = dist_from_mean_var(FDist, m, v)
            if !all(isapprox.(params(true_dist), params(new_dist), atol = 1e-8))
                @info "Mismatch:", new_dist, true_dist
                return false
            end
        end
    end
    return true
end

function test_mean_var_Geometric()
    for p ∈ 0.1:0.1:0.9
        true_dist = Geometric(p)
        m, v = mean(true_dist), var(true_dist)
        new_dist = dist_from_mean_var(Geometric, m, v)
        if !isapprox(params(true_dist)[1], params(new_dist)[1], atol = 1e-8)
            @info "Mismatch:", new_dist, true_dist
            return false
        end
    end
    return true
end

function test_mean_var_Pareto()
    for α ∈ 3.0:0.5:10.0
        for θ ∈ 0.5:0.5:5.0
            true_dist = Pareto(α, θ)
            m, v = mean(true_dist), var(true_dist)
            new_dist = dist_from_mean_var(Pareto, m, v)
            if !isapprox(mean(new_dist), m, rtol = 1e-8) || !isapprox(var(new_dist), v, rtol = 1e-8)
                @info "Mismatch:", new_dist, true_dist
                return false
            end
        end
    end
    return true
end

function test_mean_var_SymTriangularDist()
    for μ ∈ -5.0:2.5:5.0
        for s ∈ 1.0:1.0:5.0
            true_dist = SymTriangularDist(μ, s)
            m, v = mean(true_dist), var(true_dist)
            new_dist = dist_from_mean_var(SymTriangularDist, m, v)
            if !isapprox(mean(new_dist), m, atol = 1e-8) || !isapprox(var(new_dist), v, rtol = 1e-8)
                @info "Mismatch:", new_dist, true_dist
                return false
            end
        end
    end
    return true
end

function test_mean_var_DiscreteUniform()
    for a ∈ -5:5
        for b ∈ (a+1):(a+20)
            true_dist = DiscreteUniform(a, b)
            m, v = mean(true_dist), var(true_dist)
            new_dist = dist_from_mean_var(DiscreteUniform, m, v)
            if params(true_dist) != params(new_dist)
                @info "Mismatch:", new_dist, true_dist
                return false
            end
        end
    end
    return true
end

function test_mean_var_Normal()
    for μ ∈ -5.0:2.5:5.0
        for σ ∈ 0.5:0.5:3.0
            true_dist = Normal(μ, σ)
            m, v = mean(true_dist), var(true_dist)
            new_dist = dist_from_mean_var(Normal, m, v)
            if !all(isapprox.(params(true_dist), params(new_dist), atol = 1e-8))
                @info "Mismatch:", new_dist, true_dist
                return false
            end
        end
    end
    return true
end

function test_mean_var_Uniform()
    for a ∈ -5.0:2.5:5.0
        for b ∈ (a+0.5):1.0:(a+5.0)
            true_dist = Uniform(a, b)
            m, v = mean(true_dist), var(true_dist)
            new_dist = dist_from_mean_var(Uniform, m, v)
            if !isapprox(mean(new_dist), m, atol=1e-8) || !isapprox(var(new_dist), v, rtol=1e-8)
                @info "Mismatch:", new_dist, true_dist
                return false
            end
        end
    end
    return true
end

function test_mean_var_Logistic()
    for μ ∈ -5.0:2.5:5.0
        for s ∈ 0.5:0.5:3.0
            true_dist = Logistic(μ, s)
            m, v = mean(true_dist), var(true_dist)
            new_dist = dist_from_mean_var(Logistic, m, v)
            if !all(isapprox.(params(true_dist), params(new_dist), atol = 1e-8))
                @info "Mismatch:", new_dist, true_dist
                return false
            end
        end
    end
    return true
end

function test_mean_var_Laplace()
    for μ ∈ -5.0:2.5:5.0
        for b ∈ 0.5:0.5:3.0
            true_dist = Laplace(μ, b)
            m, v = mean(true_dist), var(true_dist)
            new_dist = dist_from_mean_var(Laplace, m, v)
            if !all(isapprox.(params(true_dist), params(new_dist), atol = 1e-8))
                @info "Mismatch:", new_dist, true_dist
                return false
            end
        end
    end
    return true
end

function test_mean_var_LogNormal()
    for μ_log ∈ -1.0:0.5:2.0
        for σ_log ∈ 0.1:0.3:1.5
            true_dist = LogNormal(μ_log, σ_log)
            m, v = mean(true_dist), var(true_dist)
            new_dist = dist_from_mean_var(LogNormal, m, v)
            if !isapprox(mean(new_dist), m, rtol=1e-8) || !isapprox(var(new_dist), v, rtol=1e-8)
                @info "Mismatch:", new_dist, true_dist
                return false
            end
        end
    end
    return true
end

function test_mean_var_Gamma()
    for α ∈ 0.5:0.5:5.0
        for θ ∈ 0.5:0.5:5.0
            true_dist = Gamma(α, θ)
            m, v = mean(true_dist), var(true_dist)
            new_dist = dist_from_mean_var(Gamma, m, v)
            if !all(isapprox.(params(true_dist), params(new_dist), atol = 1e-8))
                @info "Mismatch:", new_dist, true_dist
                return false
            end
        end
    end
    return true
end

function test_mean_var_Gumbel()
    for μ ∈ -5.0:2.5:5.0
        for β ∈ 0.5:0.5:3.0
            true_dist = Gumbel(μ, β)
            m, v = mean(true_dist), var(true_dist)
            new_dist = dist_from_mean_var(Gumbel, m, v)
            if !all(isapprox.(params(true_dist), params(new_dist), atol = 1e-8))
                @info "Mismatch:", new_dist, true_dist
                return false
            end
        end
    end
    return true
end

function test_mean_var_InverseGamma()
    for α ∈ 3.0:1.0:8.0
        for β ∈ 0.5:0.5:5.0
            true_dist = InverseGamma(α, β)
            m, v = mean(true_dist), var(true_dist)
            new_dist = dist_from_mean_var(InverseGamma, m, v)
            if !all(isapprox.(params(true_dist), params(new_dist), atol = 1e-8))
                @info "Mismatch:", new_dist, true_dist
                return false
            end
        end
    end
    return true
end

function test_mean_var_Poisson()
    for λ ∈ 1.0:1.0:10.0
        true_dist = Poisson(λ)
        m, v = mean(true_dist), var(true_dist)
        new_dist = dist_from_mean_var(Poisson, m, v)
        if params(true_dist) != params(new_dist)
            @info "Mismatch:", new_dist, true_dist
            return false
        end
    end
    return true
end

function test_mean_var_NegativeBinomial()
    for r ∈ 1.0:1.0:5.0
        for p ∈ 0.1:0.1:0.9
            true_dist = NegativeBinomial(r, p)
            m, v = mean(true_dist), var(true_dist)
            new_dist = dist_from_mean_var(NegativeBinomial, m, v)
            if !all(isapprox.(params(true_dist), params(new_dist), atol = 1e-8))
                @info "Mismatch:", new_dist, true_dist
                return false
            end
        end
    end
    return true
end

function test_mean_var_Binomial()
    for n ∈ [5, 10, 20, 50]
        for p ∈ [0.2, 0.5, 0.8]
            true_dist = Binomial(n, p)
            m, v = mean(true_dist), var(true_dist)
            new_dist = dist_from_mean_var(Binomial, m, v)
            n_new, p_new = params(new_dist)
            if n_new != n || !isapprox(p_new, p, atol=1e-10)
                @info "Mismatch:", new_dist, true_dist
                return false
            end
        end
    end
    return true
end

function test_mean_var_Cauchy_throws()
    try
        dist_from_mean_var(Cauchy, 0.0, 1.0)
        return false
    catch e
        return e isa DomainError
    end
end

function test_mean_var_feasibility_rejections()
    # Beta: var >= μ(1-μ)
    try; dist_from_mean_var(Beta, 0.5, 0.3); return false; catch e; e isa DomainError || return false; end
    # Exponential: var ≠ μ²
    try; dist_from_mean_var(Exponential, 2.0, 3.0); return false; catch e; e isa DomainError || return false; end
    # Poisson: var ≠ μ
    try; dist_from_mean_var(Poisson, 3.0, 5.0); return false; catch e; e isa DomainError || return false; end
    # Geometric: var ≠ μ(1+μ)
    try; dist_from_mean_var(Geometric, 2.0, 1.0); return false; catch e; e isa DomainError || return false; end
    # Gamma: μ ≤ 0
    try; dist_from_mean_var(Gamma, -1.0, 1.0); return false; catch e; e isa DomainError || return false; end
    return true
end

function test_dist_from_mean()
    # Exponential
    d = dist_from_mean(Exponential, 3.0)
    isapprox(mean(d), 3.0, atol=1e-8) || return false
    # Poisson
    d = dist_from_mean(Poisson, 5.0)
    isapprox(mean(d), 5.0, atol=1e-8) || return false
    # Rayleigh
    d = dist_from_mean(Rayleigh, 2.0)
    isapprox(mean(d), 2.0, atol=1e-8) || return false
    # Chisq
    d = dist_from_mean(Chisq, 4)
    isapprox(mean(d), 4.0, atol=1e-8) || return false
    return true
end

function test_dist_from_var_and_std()
    # Exponential: var=4 → μ=2
    d = make_dist(Exponential, var=4.0)
    isapprox(mean(d), 2.0, atol=1e-8) || return false
    # Poisson: var=3 → μ=3
    d = make_dist(Poisson, var=3.0)
    isapprox(mean(d), 3.0, atol=1e-8) || return false
    # Geometric: var=6 → μ=2
    d = make_dist(Geometric, var=6.0)
    isapprox(mean(d), 2.0, atol=1e-8) || return false
    # std=2 → var=4 → μ=2
    d = make_dist(Exponential, std=2.0)
    isapprox(mean(d), 2.0, atol=1e-8) || return false
    return true
end

function test_dist_from_mean_second_moment()
    # Gamma: μ=3, m2=12 → var=12-9=3
    d = make_dist(Gamma, mean=3.0, second_moment=12.0)
    isapprox(mean(d), 3.0, atol=1e-8) || return false
    isapprox(var(d), 3.0, rtol=1e-8) || return false
    return true
end

function test_mean_var_Geometric_dist_from_mean()
    for p ∈ 0.1:0.1:0.9
        μ = (1 - p) / p
        new_dist = dist_from_mean(Geometric, μ)
        if !isapprox(params(new_dist)[1], p, atol = 1e-8)
            @info "Mismatch: p=$p got $(params(new_dist)[1])"
            return false
        end
    end
    return true
end