# FoldedNormal: |X| where X ~ Normal(μ, σ).
#
# Not in Distributions.jl. Defined here as a first-class
# `ContinuousUnivariateDistribution` so it composes with the rest of
# Distributions.jl machinery (quantile, fit, sampling, etc.).
#
# Parameterization follows the parent Normal: `FoldedNormal(μ, σ)` means
# the parent has mean μ and standard deviation σ. Setting μ = 0 recovers
# the half-normal.

"""
    FoldedNormal(μ, σ)

The folded normal distribution: the distribution of `|X|` where
`X ~ Normal(μ, σ)`. Supported on `[0, ∞)`.

The mean and variance are
```
E[|X|]   = σ √(2/π) exp(-μ²/(2σ²)) + μ erf(μ / (σ√2))
Var[|X|] = μ² + σ² − E[|X|]²
```
"""
struct FoldedNormal{T<:Real} <: ContinuousUnivariateDistribution
    μ::T
    σ::T
    FoldedNormal{T}(μ::T, σ::T) where {T<:Real} = (σ > 0 || throw(DomainError(σ, "σ must be > 0")); new{T}(μ, σ))
end

FoldedNormal(μ::Real, σ::Real) = FoldedNormal{promote_type(typeof(μ), typeof(σ))}(promote(μ, σ)...)
FoldedNormal(σ::Real) = FoldedNormal(zero(σ), σ)
FoldedNormal() = FoldedNormal(0.0, 1.0)

Distributions.params(d::FoldedNormal) = (d.μ, d.σ)
Distributions.partype(::FoldedNormal{T}) where {T} = T
Distributions.minimum(::FoldedNormal) = 0.0
Distributions.maximum(::FoldedNormal{T}) where {T} = T(Inf)
Distributions.support(d::FoldedNormal) = Distributions.RealInterval(0.0, Inf)
Distributions.insupport(::FoldedNormal, x::Real) = x ≥ 0

function Distributions.mean(d::FoldedNormal)
    μ, σ = d.μ, d.σ
    return σ * √(2/π) * exp(-μ^2 / (2σ^2)) + μ * erf(μ / (σ * √2))
end

function Distributions.var(d::FoldedNormal)
    μ, σ = d.μ, d.σ
    m = mean(d)
    return μ^2 + σ^2 - m^2
end

Distributions.std(d::FoldedNormal) = √(var(d))

function Distributions.pdf(d::FoldedNormal, x::Real)
    x ≥ 0 || return zero(float(x))
    μ, σ = d.μ, d.σ
    c = 1 / (σ * √(2π))
    return c * (exp(-(x - μ)^2 / (2σ^2)) + exp(-(x + μ)^2 / (2σ^2)))
end

function Distributions.logpdf(d::FoldedNormal, x::Real)
    x ≥ 0 || return -Inf
    return log(pdf(d, x))
end

function Distributions.cdf(d::FoldedNormal, x::Real)
    x ≥ 0 || return zero(float(x))
    μ, σ = d.μ, d.σ
    return (erf((x - μ) / (σ * √2)) + erf((x + μ) / (σ * √2))) / 2
end

function Distributions.ccdf(d::FoldedNormal, x::Real)
    return 1 - cdf(d, x)
end

function Distributions.quantile(d::FoldedNormal, p::Real)
    0 ≤ p ≤ 1 || throw(DomainError(p, "p must be in [0,1]"))
    p == 0 && return 0.0
    p == 1 && return Inf
    # Bracket: median is upper-bounded by parent |μ| + 3σ
    hi = abs(d.μ) + 5 * d.σ
    while cdf(d, hi) < p
        hi *= 2
    end
    return find_zero(x -> cdf(d, x) - p, (0.0, hi))
end

function Base.rand(rng::Random.AbstractRNG, d::FoldedNormal)
    return abs(d.μ + d.σ * randn(rng))
end
