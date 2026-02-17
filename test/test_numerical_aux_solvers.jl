using SpecialFunctions: loggamma

# Forward evaluation: compute C = f(γ) = γ · Γ(2/γ) / Γ(1/γ)²
_f(γ) = exp(log(γ) + loggamma(2 / γ) - 2 * loggamma(1 / γ))

function test_solve_beta_ratio_known_values()
    # γ = 1 → C = 1
    isapprox(DistributionsFactories.solve_beta_ratio(1.0), 1.0; atol=1e-10) || return false
    # γ = 0.5 → C = 3
    isapprox(DistributionsFactories.solve_beta_ratio(3.0), 0.5; atol=1e-8) || return false
    # γ = 2 → C = 2/π
    isapprox(DistributionsFactories.solve_beta_ratio(2 / π), 2.0; atol=1e-8) || return false
    # γ = 0.25 → C = 35
    isapprox(DistributionsFactories.solve_beta_ratio(35.0), 0.25; atol=1e-8) || return false
    return true
end

function test_solve_beta_ratio_roundtrip()
    for γ in [0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]
        C = _f(γ)
        γ_recovered = DistributionsFactories.solve_beta_ratio(C)
        if !isapprox(γ_recovered, γ; rtol=1e-8)
            @info "Round-trip mismatch" γ C γ_recovered
            return false
        end
    end
    return true
end

function test_solve_beta_ratio_edge_cases()
    # Very large γ → C close to 1/2
    C_near_half = _f(100.0)
    γ_rec = DistributionsFactories.solve_beta_ratio(C_near_half)
    isapprox(γ_rec, 100.0; rtol=1e-6) || return false

    # Very small γ → large C
    C_large = _f(0.05)
    γ_rec = DistributionsFactories.solve_beta_ratio(C_large)
    isapprox(γ_rec, 0.05; rtol=1e-6) || return false

    return true
end

function test_solve_beta_ratio_errors()
    try
        DistributionsFactories.solve_beta_ratio(0.5)
        return false
    catch e
        e isa DomainError || return false
    end
    try
        DistributionsFactories.solve_beta_ratio(0.3)
        return false
    catch e
        e isa DomainError || return false
    end
    try
        DistributionsFactories.solve_beta_ratio(-1.0)
        return false
    catch e
        e isa DomainError || return false
    end
    try
        DistributionsFactories.solve_beta_ratio(Inf)
        return false
    catch e
        e isa DomainError || return false
    end
    return true
end
