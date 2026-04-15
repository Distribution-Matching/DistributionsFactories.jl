# --- @dist creates the right types ---

function test_dist_macro_bare_type()
    D = @dist Gamma
    return D === Gamma
end

function test_dist_macro_full_instance()
    d = @dist Normal(3.0, 2.0)
    return d isa Normal && isapprox(mean(d), 3.0, atol=1e-8)
end

function test_dist_macro_partial()
    p = @dist Gamma(3.0, _)
    return p isa PartialDist{Gamma} && p.params.α == 3.0 && ismissing(p.params.θ)
end

# --- PartialDist with dist_from_mean ---

function test_partial_dist_from_mean()
    d = dist_from_mean(@dist(Gamma(3.0, _)), 5.0)
    isapprox(mean(d), 5.0, atol=1e-8) || return false
    isapprox(params(d)[1], 3.0, atol=1e-8) || return false
    return true
end

function test_partial_normal_fix_sigma()
    d = dist_from_mean(@dist(Normal(_, 1.0)), 3.0)
    isapprox(mean(d), 3.0, atol=1e-8) || return false
    isapprox(params(d)[2], 1.0, atol=1e-8) || return false
    return true
end

# --- PartialDist with dist_from_mean_var ---

function test_partial_dist_from_mean_var()
    d = dist_from_mean_var(@dist(Normal(0.0, _)), 0.0, 4.0)
    isapprox(params(d)[1], 0.0, atol=1e-8) || return false
    isapprox(var(d), 4.0, rtol=1e-8) || return false
    return true
end

function test_partial_all_missing_delegates()
    # All params missing → same as type dispatch
    d = dist_from_mean_var(@dist(Gamma(_, _)), 5.0, 3.0)
    isapprox(mean(d), 5.0, atol=1e-8) || return false
    isapprox(var(d), 3.0, rtol=1e-8) || return false
    return true
end

function test_partial_beta()
    d = dist_from_mean(@dist(Beta(2.0, _)), 0.4)
    isapprox(mean(d), 0.4, atol=1e-6) || return false
    isapprox(params(d)[1], 2.0, atol=1e-8) || return false
    return true
end

# --- Full instance with dist_from_mean_var (LocationScale) ---

function test_instance_tdist()
    d = dist_from_mean_var(@dist(TDist(7)), 5.0, 2.0)
    isapprox(mean(d), 5.0, atol=1e-8) || return false
    isapprox(var(d), 2.0, rtol=1e-8) || return false
    return true
end

# --- dist_from_var and convenience wrappers with PartialDist ---

function test_partial_dist_from_var()
    # Logistic(μ=2, _): fix location, solve scale from variance
    d = dist_from_var(@dist(Logistic(2.0, _)), 22.3)
    isapprox(var(d), 22.3, rtol=1e-6) || return false
    isapprox(params(d)[1], 2.0, atol=1e-8) || return false
    return true
end

function test_partial_dist_from_std()
    d = dist_from_std(@dist(Logistic(2.0, _)), 3.0)
    isapprox(std(d), 3.0, rtol=1e-6) || return false
    isapprox(params(d)[1], 2.0, atol=1e-8) || return false
    return true
end

function test_partial_dist_from_mean_cv()
    # Gamma(α=4, θ=?): with α=4, CV = 1/√α = 0.5, so mean=5 → θ=5/4, var=5²*0.25=6.25
    d = dist_from_mean_cv(@dist(Gamma(4.0, _)), 5.0, 0.5)
    isapprox(mean(d), 5.0, atol=1e-6) || return false
    isapprox(std(d) / mean(d), 0.5, rtol=1e-6) || return false
    isapprox(params(d)[1], 4.0, atol=1e-8) || return false
    return true
end

function test_partial_dist_from_mean_std()
    d = dist_from_mean_std(@dist(Normal(_, 2.0)), 3.0, 2.0)
    isapprox(mean(d), 3.0, atol=1e-8) || return false
    isapprox(std(d), 2.0, rtol=1e-8) || return false
    return true
end

# --- Bare type with dist_from_mean_var ---

function test_type_via_macro()
    d = dist_from_mean_var(@dist(Gamma), 5.0, 3.0)
    isapprox(mean(d), 5.0, atol=1e-8) || return false
    isapprox(var(d), 3.0, rtol=1e-8) || return false
    return true
end
