# DiscreteTriangular: integer-valued asymmetric triangular PMF.
#
# Parameters: a ≤ c ≤ b (all integers). Support is {a, a+1, …, b}.
# PMF is two linear ramps meeting at the mode c:
#     P(k) ∝ (k − a + 1) / (c − a + 1)   for k ∈ [a, c]
#     P(k) ∝ (b − k + 1) / (b − c + 1)   for k ∈ [c, b]
# Both pieces equal 1 at k = c. Normalising constant: Z = (b − a + 2) / 2.

"""
    DiscreteTriangular(a::Integer, b::Integer, c::Integer)

Integer-valued (asymmetric) triangular distribution on `{a, a+1, …, b}`
with mode at `c`, where `a ≤ c ≤ b`. The PMF is two linear ramps
meeting at `c`. When `b − c == c − a`, the distribution is symmetric
(see [`DiscreteSymmetricTriangular`](@ref) for the closed-form
parameterization).
"""
struct DiscreteTriangular <: DiscreteUnivariateDistribution
    a::Int
    b::Int
    c::Int
    function DiscreteTriangular(a::Integer, b::Integer, c::Integer)
        a ≤ c ≤ b || throw(DomainError((a, b, c), "must satisfy a ≤ c ≤ b"))
        return new(Int(a), Int(b), Int(c))
    end
end

Distributions.params(d::DiscreteTriangular) = (d.a, d.b, d.c)
Distributions.minimum(d::DiscreteTriangular) = d.a
Distributions.maximum(d::DiscreteTriangular) = d.b
Distributions.support(d::DiscreteTriangular) = d.a:d.b
Distributions.insupport(d::DiscreteTriangular, x::Real) =
    isinteger(x) && d.a ≤ x ≤ d.b
Distributions.mode(d::DiscreteTriangular) = d.c

function Distributions.pdf(d::DiscreteTriangular, x::Real)
    insupport(d, x) || return 0.0
    k = Int(x)
    Z = (d.b - d.a + 2) / 2
    if k ≤ d.c
        return (k - d.a + 1) / (d.c - d.a + 1) / Z
    else
        return (d.b - k + 1) / (d.b - d.c + 1) / Z
    end
end

Distributions.logpdf(d::DiscreteTriangular, x::Real) =
    insupport(d, x) ? log(pdf(d, x)) : -Inf

function Distributions.cdf(d::DiscreteTriangular, x::Real)
    x < d.a && return 0.0
    x ≥ d.b && return 1.0
    s = 0.0
    for k in d.a:floor(Int, x)
        s += pdf(d, k)
    end
    return s
end

function Distributions.mean(d::DiscreteTriangular)
    s = 0.0
    for k in d.a:d.b
        s += k * pdf(d, k)
    end
    return s
end

function Distributions.var(d::DiscreteTriangular)
    m = mean(d)
    s = 0.0
    for k in d.a:d.b
        s += (k - m)^2 * pdf(d, k)
    end
    return s
end

Distributions.std(d::DiscreteTriangular) = √(var(d))

function Distributions.quantile(d::DiscreteTriangular, p::Real)
    0 ≤ p ≤ 1 || throw(DomainError(p, "p must be in [0,1]"))
    p == 0 && return d.a
    p == 1 && return d.b
    s = 0.0
    for k in d.a:d.b
        s += pdf(d, k)
        s ≥ p && return k
    end
    return d.b
end

function Base.rand(rng::Random.AbstractRNG, d::DiscreteTriangular)
    return quantile(d, rand(rng))
end
