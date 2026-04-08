function test_mean_std_roundtrip()
    for (D, params_list) in [
        (Normal,    [(0.0, 1.0), (5.0, 2.0), (-3.0, 0.5)]),
        (Gamma,     [(1.0, 1.0), (5.0, 2.0), (10.0, 0.1)]),
        (LogNormal, [(1.0, 0.5), (3.0, 1.0)]),
    ]
        for (μ, σ) in params_list
            d = dist_from_mean_std(D, μ, σ)
            if !isapprox(mean(d), μ, atol=1e-8) || !isapprox(std(d), σ, rtol=1e-8)
                @info "Mismatch:", D, μ, σ, mean(d), std(d)
                return false
            end
        end
    end
    return true
end

function test_mean_cv_roundtrip()
    for (D, params_list) in [
        (Gamma,     [(2.0, 0.5), (5.0, 1.0), (1.0, 2.0)]),
        (LogNormal, [(2.0, 0.5), (1.0, 1.0)]),
        (Weibull,   [(3.0, 0.3), (1.0, 0.8)]),
    ]
        for (μ, cv) in params_list
            d = dist_from_mean_cv(D, μ, cv)
            actual_cv = std(d) / mean(d)
            if !isapprox(mean(d), μ, atol=1e-8) || !isapprox(actual_cv, cv, rtol=1e-8)
                @info "Mismatch:", D, μ, cv, mean(d), actual_cv
                return false
            end
        end
    end
    return true
end

function test_mean_scv_roundtrip()
    for (D, params_list) in [
        (Gamma,     [(2.0, 0.25), (5.0, 1.0), (1.0, 4.0)]),
        (LogNormal, [(2.0, 0.25), (1.0, 1.0)]),
    ]
        for (μ, scv) in params_list
            d = dist_from_mean_scv(D, μ, scv)
            actual_scv = var(d) / mean(d)^2
            if !isapprox(mean(d), μ, atol=1e-8) || !isapprox(actual_scv, scv, rtol=1e-8)
                @info "Mismatch:", D, μ, scv, mean(d), actual_scv
                return false
            end
        end
    end
    return true
end

function test_variants_consistency()
    μ = 3.0
    σ = 1.5
    v = σ^2
    cv = σ / μ
    scv = v / μ^2

    d_var = dist_from_mean_var(Gamma, μ, v)
    d_std = dist_from_mean_std(Gamma, μ, σ)
    d_cv  = dist_from_mean_cv(Gamma, μ, cv)
    d_scv = dist_from_mean_scv(Gamma, μ, scv)

    for d in [d_std, d_cv, d_scv]
        if !all(isapprox.(params(d_var), params(d), atol=1e-10))
            @info "Inconsistency:", params(d_var), params(d)
            return false
        end
    end
    return true
end

function test_exists_variants()
    # Valid cases
    exists_unique_dist_from_mean_std(Gamma, 2.0, 1.0) || return false
    exists_unique_dist_from_mean_cv(Gamma, 2.0, 0.5) || return false
    exists_unique_dist_from_mean_scv(Gamma, 2.0, 0.25) || return false

    # Invalid: negative mean for Gamma propagates through
    for f in [exists_unique_dist_from_mean_std, exists_unique_dist_from_mean_cv, exists_unique_dist_from_mean_scv]
        try
            f(Gamma, -1.0, 1.0)
            return false
        catch e
            e isa DomainError || return false
        end
    end
    return true
end
