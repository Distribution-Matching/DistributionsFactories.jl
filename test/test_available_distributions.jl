function test_available_by_interval()
    # IntervalSets syntax
    real = available_distributions(-Inf..Inf)
    pos  = available_distributions(0..Inf)
    unit = available_distributions(0..1)
    Normal ∈ real && Gamma ∈ pos && Beta ∈ unit || return false

    # Bounded not yet supported
    isempty(available_distributions(2..7)) || return false
    return true
end

function test_available_by_realinterval()
    # Distributions.jl support()
    real = available_distributions(support(Normal()))
    pos  = available_distributions(support(Gamma(1,1)))
    unit = available_distributions(support(Beta(2,3)))
    Normal ∈ real && Gamma ∈ pos && Beta ∈ unit || return false

    # Should match IntervalSets results
    real == available_distributions(-Inf..Inf) || return false
    pos  == available_distributions(0..Inf)    || return false
    unit == available_distributions(0..1)      || return false
    return true
end

function test_available_by_range()
    bounded = available_distributions(0:10)
    Binomial ∈ bounded && DiscreteUniform ∈ bounded || return false
    return true
end

function test_available_with_mean_var()
    feasible = available_distributions(0..Inf, mean=5.0, var=3.0)
    Gamma ∈ feasible || return false
    # Exponential requires var=μ², so var=3 with μ=5 should exclude it
    Exponential ∉ feasible || return false
    return true
end

function test_available_with_mean_std()
    feasible = available_distributions(0..1, mean=0.5, std=0.1)
    Beta ∈ feasible || return false
    return true
end

function test_available_with_mean_cv()
    feasible = available_distributions(0..Inf, mean=5.0, cv=1.0)
    Exponential ∈ feasible || return false
    return true
end

function test_available_kwargs_only()
    feasible = available_distributions(mean=5.0, var=25.0)
    # Should include distributions from multiple supports
    Normal ∈ feasible && Gamma ∈ feasible && Exponential ∈ feasible || return false
    return true
end

function test_available_kwargs_requires_dispersion()
    # Calling with no kwargs at all is still rejected; mean-alone is now
    # supported (returns the candidates whose dist_from_mean works).
    try
        available_distributions()
        return false
    catch e
        return e isa ArgumentError
    end
end
