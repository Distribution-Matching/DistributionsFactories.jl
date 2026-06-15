function test_quantile_exponential()
    for θ in [0.5, 1.0, 3.0, 10.0]
        true_d = Exponential(θ)
        for p in [0.1, 0.25, 0.5, 0.75, 0.9]
            q = quantile(true_d, p)
            d = dist_from_quantile(Exponential, p, q)
            if !isapprox(params(d)[1], θ, rtol=1e-8)
                @info "Mismatch:", θ, p, params(d)
                return false
            end
        end
    end
    return true
end

function test_median_exponential()
    for θ in [0.5, 1.0, 5.0]
        true_d = Exponential(θ)
        d = make_dist(Exponential, median=median(true_d))
        if !isapprox(params(d)[1], θ, rtol=1e-8)
            @info "Mismatch:", θ, params(d)
            return false
        end
    end
    return true
end

function test_q1_q3_exponential()
    for θ in [1.0, 3.0, 7.0]
        true_d = Exponential(θ)
        d1 = make_dist(Exponential, q1=quantile(true_d, 0.25))
        d3 = make_dist(Exponential, q3=quantile(true_d, 0.75))
        if !isapprox(params(d1)[1], θ, rtol=1e-8) || !isapprox(params(d3)[1], θ, rtol=1e-8)
            @info "Mismatch:", θ, params(d1), params(d3)
            return false
        end
    end
    return true
end

function test_quantiles_location_scale()
    cases = [
        (Normal,   [(0.0, 1.0), (5.0, 2.0), (-3.0, 0.5)]),
        (Laplace,  [(0.0, 1.0), (2.0, 3.0)]),
        (Logistic, [(0.0, 1.0), (-1.0, 2.0)]),
        (Cauchy,   [(0.0, 1.0), (3.0, 2.0)]),
        (Gumbel,   [(0.0, 1.0), (5.0, 2.0)]),
    ]
    for (D, param_list) in cases
        for p in param_list
            true_d = D(p...)
            q25 = quantile(true_d, 0.25)
            q75 = quantile(true_d, 0.75)
            d = dist_from_quantiles(D, 0.25, q25, 0.75, q75)
            if !all(isapprox.(params(d), p, rtol=1e-8))
                @info "Mismatch:", D, p, params(d)
                return false
            end
        end
    end
    return true
end

function test_quantiles_gamma()
    for α in [0.5, 1.0, 3.0, 10.0]
        for θ in [0.5, 2.0, 5.0]
            true_d = Gamma(α, θ)
            q10 = quantile(true_d, 0.1)
            q90 = quantile(true_d, 0.9)
            d = dist_from_quantiles(Gamma, 0.1, q10, 0.9, q90)
            if !all(isapprox.(params(d), (α, θ), rtol=1e-6))
                @info "Mismatch:", α, θ, params(d)
                return false
            end
        end
    end
    return true
end

function test_quantiles_beta()
    for α in [0.5, 1.0, 2.0, 5.0]
        for β in [0.5, 1.0, 3.0, 8.0]
            true_d = Beta(α, β)
            q25 = quantile(true_d, 0.25)
            q75 = quantile(true_d, 0.75)
            d = dist_from_quantiles(Beta, 0.25, q25, 0.75, q75)
            if !all(isapprox.(params(d), (α, β), rtol=1e-4))
                @info "Mismatch:", α, β, params(d)
                return false
            end
        end
    end
    return true
end

function test_mean_quantile_gamma()
    for α in [1.0, 3.0, 8.0]
        for θ in [0.5, 2.0]
            true_d = Gamma(α, θ)
            d = dist_from_mean_quantile(Gamma, mean(true_d), 0.75, quantile(true_d, 0.75))
            if !all(isapprox.(params(d), (α, θ), rtol=1e-6))
                @info "Mismatch:", α, θ, params(d)
                return false
            end
        end
    end
    return true
end

function test_mean_quantile_beta()
    for α in [2.0, 3.0, 5.0]
        for β in [1.5, 3.0, 8.0]
            true_d = Beta(α, β)
            d = dist_from_mean_quantile(Beta, mean(true_d), 0.75, quantile(true_d, 0.75))
            if !all(isapprox.(params(d), (α, β), rtol=1e-4))
                @info "Mismatch:", α, β, params(d)
                return false
            end
        end
    end
    return true
end

function test_median_iqr_normal()
    for (μ, σ) in [(0.0, 1.0), (5.0, 2.0), (-3.0, 0.5)]
        true_d = Normal(μ, σ)
        med = median(true_d)
        iqr = quantile(true_d, 0.75) - quantile(true_d, 0.25)
        d = make_dist(Normal, median=med, iqr=iqr)
        if !all(isapprox.(params(d), (μ, σ), rtol=1e-8))
            @info "Mismatch:", μ, σ, params(d)
            return false
        end
    end
    return true
end

function test_mean_median_gamma()
    for α in [2.0, 4.0, 10.0]
        for θ in [1.0, 2.5]
            true_d = Gamma(α, θ)
            d = make_dist(Gamma, mean=mean(true_d), median=median(true_d))
            if !all(isapprox.(params(d), (α, θ), rtol=1e-6))
                @info "Mismatch:", α, θ, params(d)
                return false
            end
        end
    end
    return true
end

function test_quantiles_inverse_gamma()
    # Includes small-shape (α < 1) heavy-tailed cases: the scale-family solver
    # must bracket the root rather than overshoot into invalid shapes.
    for (α, θ) in [(0.2, 1.0), (0.3, 2.0), (2.0, 1.0), (3.5, 2.0), (1.5, 5.0), (5.0, 0.5)]
        d0 = InverseGamma(α, θ)
        q1 = quantile(d0, 0.1)
        q2 = quantile(d0, 0.9)
        d = dist_from_quantiles(InverseGamma, 0.1, q1, 0.9, q2)
        if !all(isapprox.(params(d), (α, θ), rtol=1e-5))
            @info "Mismatch:" α θ params(d)
            return false
        end
    end
    return true
end

function test_quantiles_uniform()
    for (a, b) in [(0.0, 1.0), (-2.0, 5.0), (3.0, 10.0), (-7.0, -1.0)]
        d0 = Uniform(a, b)
        q1 = quantile(d0, 0.2)
        q2 = quantile(d0, 0.8)
        d = dist_from_quantiles(Uniform, 0.2, q1, 0.8, q2)
        if !all(isapprox.(params(d), (a, b), rtol=1e-8, atol=1e-10))
            @info "Mismatch:" a b params(d)
            return false
        end
    end
    return true
end

function test_quantiles_frechet()
    for (α, θ) in [(2.0, 3.0), (1.5, 1.0), (4.0, 2.5), (0.8, 5.0)]
        d0 = Frechet(α, θ)
        q1 = quantile(d0, 0.25)
        q2 = quantile(d0, 0.75)
        d = dist_from_quantiles(Frechet, 0.25, q1, 0.75, q2)
        if !all(isapprox.(params(d), (α, θ), rtol=1e-6))
            @info "Mismatch:" α θ params(d)
            return false
        end
    end
    return true
end

function test_quantiles_symtriangular()
    for (μ, σ) in [(0.0, 1.0), (2.0, 3.0), (-1.0, 0.5), (5.0, 2.0)]
        d0 = SymTriangularDist(μ, σ)
        q1 = quantile(d0, 0.1)
        q2 = quantile(d0, 0.9)
        d = dist_from_quantiles(SymTriangularDist, 0.1, q1, 0.9, q2)
        if !all(isapprox.(params(d), (μ, σ), rtol=1e-8, atol=1e-10))
            @info "Mismatch:" μ σ params(d)
            return false
        end
    end
    return true
end

function test_quantile_chisq()
    for k in [1.5, 3.0, 7.2, 10.0]
        d0 = Chisq(k)
        for p in [0.2, 0.5, 0.8]
            q = quantile(d0, p)
            d = dist_from_quantile(Chisq, p, q)
            if !isapprox(params(d)[1], k, rtol=1e-5)
                @info "Mismatch:" k p params(d)
                return false
            end
        end
    end
    return true
end

function test_quantile_chi()
    for k in [1.5, 3.0, 7.2, 10.0]
        d0 = Chi(k)
        for p in [0.2, 0.5, 0.8]
            q = quantile(d0, p)
            d = dist_from_quantile(Chi, p, q)
            if !isapprox(params(d)[1], k, rtol=1e-5)
                @info "Mismatch:" k p params(d)
                return false
            end
        end
    end
    return true
end

function test_quantiles_make_dist_cdf()
    # Public entry point + equal-tail prior recipe: 1% below 0.5, 1% above 50.
    d = make_dist(InverseGamma, quantiles=[(0.01, 0.5), (0.99, 50.0)])
    ok_cdf = isapprox(cdf(d, 0.5), 0.01, atol=1e-6) &&
             isapprox(cdf(d, 50.0), 0.99, atol=1e-6)
    ok_exists = dist_exists(InverseGamma, quantiles=[(0.01, 0.5), (0.99, 50.0)]) === true
    return ok_cdf && ok_exists
end
