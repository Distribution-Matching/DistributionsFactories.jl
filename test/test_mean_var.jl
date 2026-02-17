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

function test_mean_var_FDist()
    for d₁ ∈ 2:0.5:50
        for d₂ ∈ 4.5:0.5:50
            true_dist = FDist(d₁, d₂)
            m, v = mean(true_dist), var(true_dist)
            # @show d₁, d₂, m, v
            new_dist = dist_from_mean_var(FDist, m, v)
            if !all(isapprox.(params(true_dist), params(new_dist), atol = 1e-8))
                @info "Mismatch:", new_dist, true_dist
                return false
            end
        end
    end
    return true
end