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