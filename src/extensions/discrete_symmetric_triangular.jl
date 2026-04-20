# DiscreteSymmetricTriangular: integer-valued symmetric triangular PMF.
#
# Parameters: center μ ∈ ℤ, half-width n ∈ ℕ₀.
# Support: {μ-n, μ-n+1, …, μ+n}.
# PMF:    P(μ + k) = (n + 1 − |k|) / (n + 1)²  for k ∈ {-n, …, n}.
# Moments (closed form):
#     E[X]   = μ
#     Var[X] = n(n+2)/6

"""
    DiscreteSymmetricTriangular(μ::Integer, n::Integer)

Integer-valued symmetric triangular distribution centered at `μ` with
half-width `n ≥ 0`. Support is `{μ-n, μ-n+1, …, μ+n}` and
`P(μ + k) = (n + 1 − |k|) / (n + 1)²`.

Mean is `μ`; variance is `n(n+2)/6`.
"""
struct DiscreteSymmetricTriangular <: DiscreteUnivariateDistribution
    μ::Int
    n::Int
    function DiscreteSymmetricTriangular(μ::Integer, n::Integer)
        n ≥ 0 || throw(DomainError(n, "half-width n must be ≥ 0"))
        return new(Int(μ), Int(n))
    end
end

Distributions.params(d::DiscreteSymmetricTriangular) = (d.μ, d.n)
Distributions.minimum(d::DiscreteSymmetricTriangular) = d.μ - d.n
Distributions.maximum(d::DiscreteSymmetricTriangular) = d.μ + d.n
Distributions.support(d::DiscreteSymmetricTriangular) = (d.μ - d.n):(d.μ + d.n)
Distributions.insupport(d::DiscreteSymmetricTriangular, x::Real) =
    isinteger(x) && (d.μ - d.n) ≤ x ≤ (d.μ + d.n)

Distributions.mean(d::DiscreteSymmetricTriangular) = float(d.μ)
Distributions.var(d::DiscreteSymmetricTriangular) = d.n * (d.n + 2) / 6
Distributions.std(d::DiscreteSymmetricTriangular) = √(var(d))
Distributions.mode(d::DiscreteSymmetricTriangular) = d.μ

function Distributions.pdf(d::DiscreteSymmetricTriangular, x::Real)
    insupport(d, x) || return 0.0
    k = Int(x) - d.μ
    return (d.n + 1 - abs(k)) / (d.n + 1)^2
end

Distributions.logpdf(d::DiscreteSymmetricTriangular, x::Real) =
    insupport(d, x) ? log(pdf(d, x)) : -Inf

function Distributions.cdf(d::DiscreteSymmetricTriangular, x::Real)
    x < d.μ - d.n && return 0.0
    x ≥ d.μ + d.n && return 1.0
    k = floor(Int, x) - d.μ                # supported k value at floor(x)
    # Sum P(μ+j) for j = -n .. k.
    # By symmetry, S(k) = S(-k-1) reflected; just compute directly.
    s = 0.0
    for j in -d.n:k
        s += (d.n + 1 - abs(j)) / (d.n + 1)^2
    end
    return s
end

function Distributions.quantile(d::DiscreteSymmetricTriangular, p::Real)
    0 ≤ p ≤ 1 || throw(DomainError(p, "p must be in [0,1]"))
    p == 0 && return d.μ - d.n
    p == 1 && return d.μ + d.n
    s = 0.0
    for k in -d.n:d.n
        s += (d.n + 1 - abs(k)) / (d.n + 1)^2
        s ≥ p && return d.μ + k
    end
    return d.μ + d.n
end

function Base.rand(rng::Random.AbstractRNG, d::DiscreteSymmetricTriangular)
    # X = ⌊U₁ + U₂⌋ over integer half-widths gives a symmetric triangular shape.
    # Implement via inverse CDF on U(0,1) — fine for moderate n.
    return quantile(d, rand(rng))
end
