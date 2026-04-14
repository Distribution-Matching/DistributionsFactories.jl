# Single-parameter distributions: construct from mean alone

"""
    dist_from_mean(D, μ̄)

Construct a 1-parameter distribution `D` from its mean `μ̄` alone. Only supported for
distributions where all parameters are determined by the mean:
`Exponential`, `Poisson`, `Rayleigh`, `Chisq`, `Geometric`.
"""
function dist_from_mean end

function dist_from_mean(::Type{Exponential}, μ̄::Number)
    μ̄ > 0 || throw(DomainError(μ̄, "Exponential: μ̄ must be > 0"))
    return Exponential(μ̄)
end

function dist_from_mean(::Type{Poisson}, μ̄::Number)
    μ̄ > 0 || throw(DomainError(μ̄, "Poisson: μ̄ must be > 0"))
    return Poisson(μ̄)
end

function dist_from_mean(::Type{Rayleigh}, μ̄::Number)
    μ̄ > 0 || throw(DomainError(μ̄, "Rayleigh: μ̄ must be > 0"))
    σ = μ̄ / √(π / 2)
    return Rayleigh(σ)
end

function dist_from_mean(::Type{Chisq}, μ̄::Number)
    μ̄ > 0 || throw(DomainError(μ̄, "Chisq: μ̄ must be > 0"))
    isinteger(μ̄) || throw(DomainError(μ̄, "Chisq: μ̄ must be a positive integer"))
    return Chisq(μ̄)
end

function dist_from_mean(::Type{Geometric}, μ̄::Number)
    μ̄ > 0 || throw(DomainError(μ̄, "Geometric: μ̄ must be > 0"))
    p = 1 / (1 + μ̄)
    return Geometric(p)
end
