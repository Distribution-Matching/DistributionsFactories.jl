# Single-parameter distributions: construct from mean alone

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
