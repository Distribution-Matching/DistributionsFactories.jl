# Single-parameter distributions: construct from mean alone

"""
    dist_from_mean(D, μ)

Construct a 1-parameter distribution `D` from its mean alone. Only supported for
distributions where all parameters are determined by the mean:
`Exponential`, `Poisson`, `Rayleigh`, `Chisq`, `Geometric`.
"""
function dist_from_mean end

function dist_from_mean(::Type{Exponential}, μ::Number)
    μ > 0 || throw(DomainError(μ, "Exponential: μ must be > 0"))
    return Exponential(μ)
end

function dist_from_mean(::Type{Poisson}, μ::Number)
    μ > 0 || throw(DomainError(μ, "Poisson: μ must be > 0"))
    return Poisson(μ)
end

function dist_from_mean(::Type{Rayleigh}, μ::Number)
    μ > 0 || throw(DomainError(μ, "Rayleigh: μ must be > 0"))
    σ = μ / √(π / 2)
    return Rayleigh(σ)
end

function dist_from_mean(::Type{Chisq}, μ::Number)
    μ > 0 || throw(DomainError(μ, "Chisq: μ must be > 0"))
    isinteger(μ) || throw(DomainError(μ, "Chisq: μ must be a positive integer"))
    return Chisq(μ)
end

function dist_from_mean(::Type{Geometric}, μ::Number)
    μ > 0 || throw(DomainError(μ, "Geometric: μ must be > 0"))
    p = 1 / (1 + μ)
    return Geometric(p)
end
