# --- Moment-based ---

function test_make_dist_mean_var()
    d = make_dist(Gamma, mean=5.0, var=3.0)
    isapprox(mean(d), 5.0, atol=1e-8) || return false
    isapprox(var(d), 3.0, rtol=1e-8) || return false
    return true
end

function test_make_dist_mean_std()
    d = make_dist(Normal, mean=0.0, std=2.0)
    isapprox(std(d), 2.0, rtol=1e-8) || return false
    return true
end

function test_make_dist_mean_cv()
    d = make_dist(Gamma, mean=5.0, cv=0.5)
    isapprox(std(d) / mean(d), 0.5, rtol=1e-8) || return false
    return true
end

function test_make_dist_mean_scv()
    d = make_dist(Gamma, mean=5.0, scv=0.25)
    isapprox(var(d) / mean(d)^2, 0.25, rtol=1e-8) || return false
    return true
end

function test_make_dist_mean_only()
    d = make_dist(Exponential, mean=3.0)
    isapprox(mean(d), 3.0, atol=1e-8) || return false
    return true
end

function test_make_dist_var_only()
    d = make_dist(Exponential, var=4.0)
    isapprox(mean(d), 2.0, atol=1e-8) || return false
    return true
end

# --- Quantile-based ---

function test_make_dist_median()
    d = make_dist(Exponential, median=2.0)
    isapprox(median(d), 2.0, rtol=1e-8) || return false
    return true
end

function test_make_dist_q1_q3()
    d = make_dist(Normal, q1=10.0, q3=30.0)
    isapprox(quantile(d, 0.25), 10.0, rtol=1e-8) || return false
    isapprox(quantile(d, 0.75), 30.0, rtol=1e-8) || return false
    return true
end

function test_make_dist_mean_median()
    d = make_dist(Beta, mean=0.4, median=0.35)
    isapprox(mean(d), 0.4, atol=1e-4) || return false
    isapprox(median(d), 0.35, rtol=1e-4) || return false
    return true
end

function test_make_dist_median_iqr()
    d = make_dist(Normal, median=5.0, iqr=4.0)
    isapprox(median(d), 5.0, atol=1e-8) || return false
    iqr_actual = quantile(d, 0.75) - quantile(d, 0.25)
    isapprox(iqr_actual, 4.0, rtol=1e-8) || return false
    return true
end

function test_make_dist_quantiles()
    d = make_dist(Gamma, quantiles=[(0.1, 1.0), (0.9, 10.0)])
    isapprox(quantile(d, 0.1), 1.0, rtol=1e-6) || return false
    isapprox(quantile(d, 0.9), 10.0, rtol=1e-6) || return false
    return true
end

# --- With support ---

function test_make_dist_support_affine()
    d = make_dist(Beta, mean=3.5, var=0.5, support=2..7)
    isapprox(mean(d), 3.5, rtol=1e-8) || return false
    isapprox(minimum(d), 2.0, atol=1e-10) || return false
    isapprox(maximum(d), 7.0, atol=1e-10) || return false
    return true
end

function test_make_dist_support_shift()
    d = make_dist(Gamma, mean=8.0, var=3.0, support=3..Inf)
    isapprox(mean(d), 8.0, rtol=1e-8) || return false
    isapprox(minimum(d), 3.0, atol=1e-10) || return false
    return true
end

# --- With @dist ---

function test_make_dist_partial_mean()
    spec = @dist Gamma(3.0, _)
    d = make_dist(spec, mean=5.0)
    isapprox(mean(d), 5.0, atol=1e-8) || return false
    isapprox(params(d)[1], 3.0, atol=1e-8) || return false
    return true
end

function test_make_dist_partial_var()
    spec = @dist Logistic(2.0, _)
    d = make_dist(spec, var=22.3)
    isapprox(var(d), 22.3, rtol=1e-6) || return false
    isapprox(params(d)[1], 2.0, atol=1e-8) || return false
    return true
end

# --- dist_exists ---

function test_dist_exists_true()
    return dist_exists(Beta, mean=0.5, var=0.1)
end

# --- Mode-based ---

function test_make_dist_mode_rayleigh()
    d = make_dist(Rayleigh, mode=2.0)
    isapprox(mode(d), 2.0, atol=1e-8) || return false
    isapprox(params(d)[1], 2.0, atol=1e-8) || return false
    return true
end

function test_make_dist_mean_mode_gamma()
    d = make_dist(Gamma, mean=5.0, mode=3.0)
    isapprox(mean(d), 5.0, atol=1e-8) || return false
    isapprox(mode(d), 3.0, atol=1e-8) || return false
    return true
end

function test_make_dist_mean_mode_beta()
    d = make_dist(Beta, mean=0.4, mode=0.3)
    isapprox(mean(d), 0.4, atol=1e-6) || return false
    isapprox(mode(d), 0.3, atol=1e-6) || return false
    return true
end

function test_make_dist_mode_var_normal()
    d = make_dist(Normal, mode=3.0, var=4.0)
    isapprox(mean(d), 3.0, atol=1e-8) || return false
    isapprox(var(d), 4.0, rtol=1e-8) || return false
    return true
end

function test_make_dist_mode_var_gamma()
    d = make_dist(Gamma, mode=2.0, var=3.0)
    isapprox(mode(d), 2.0, atol=1e-6) || return false
    isapprox(var(d), 3.0, rtol=1e-6) || return false
    return true
end

# --- dist_exists ---

function test_dist_exists_false()
    # dist_exists is a pure predicate: infeasible moments → false, never throws.
    return dist_exists(Exponential, mean=2.5, var=1.5) == false
end
