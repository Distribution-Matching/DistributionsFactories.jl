using Distributions
using Roots
using Polynomials
using SpecialFunctions


"""
    dist_from_mean_var(D, μ̄, σ̄²)

Construct a distribution of type `D` with the given mean `μ̄` and variance `σ̄²`.

Dispatches on the distribution type (or instance for truncated/TDist).
Throws `DomainError` if no valid distribution exists for the given moments.
Use `exists_unique_dist_from_mean_var` to check feasibility before calling.

See also: [`dist_from_mean_std`](@ref), [`dist_from_mean_cv`](@ref),
[`exists_unique_dist_from_mean_var`](@ref)
"""
function dist_from_mean_var end

"""
    dist_from_mean_var(::Type{Beta}, μ̄, σ̄²)

Direct formula. Construct a `Beta(α, β)` distribution.
Requires `0 < μ̄ < 1` and `0 < σ̄² < μ̄(1-μ̄)`.
"""
function dist_from_mean_var(::Type{Beta}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(Beta, μ̄, σ̄²)
    S = (μ̄*(1-μ̄))/σ̄²-1
    α = μ̄*S
    β = (1-μ̄)*S
    return Beta(α,β)
end

"""
    dist_from_mean_var(::Type{Uniform}, μ̄, σ̄²)

Direct formula. Construct a `Uniform(a, b)` distribution. Any `μ̄ ∈ ℝ` and `σ̄² > 0`.
"""
function dist_from_mean_var(::Type{Uniform}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(Uniform, μ̄, σ̄²)
    diff = √(3*σ̄²)
    a = μ̄-diff
    b = μ̄+diff
    return Uniform(a,b)
end

"""
    dist_from_mean_var(::Type{Normal}, μ̄, σ̄²)

Direct formula. Construct a `Normal(μ, σ)` distribution. Any `μ̄ ∈ ℝ` and `σ̄² > 0`.
"""
function dist_from_mean_var(::Type{Normal}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(Normal, μ̄, σ̄²)
    return Normal(μ̄,√(σ̄²))
end

"""
    dist_from_mean_var(::Type{TDist}, μ̄, σ̄²)

Direct formula. Construct a standard `TDist(ν)` distribution.
Requires `μ̄ = 0` and `σ̄² > 1`.
"""
function dist_from_mean_var(::Type{TDist}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(TDist, μ̄, σ̄²)
    v=2*σ̄²/(σ̄²-1)
    return TDist(v)
end

"""
    dist_from_mean_var(d::TDist, μ̄, σ̄²)

Direct formula. Construct an affine-transformed `TDist` (location-scale) with
arbitrary mean and variance. The input `d` provides the degrees of freedom `ν`
(must be > 2). Returns `μ̄ + σ * d` as a `LocationScale` distribution.
"""
function dist_from_mean_var(d::TDist, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(d, μ̄, σ̄²)
    ν = dof(d)
    base_var = ν / (ν - 2)
    σ = √(σ̄² / base_var)
    return μ̄ + σ * d
end

"""
    dist_from_mean_var(::Type{Cauchy}, μ̄, σ̄²)

Infeasible. Always throws `DomainError` — the Cauchy distribution has no defined
mean or variance. Use quantile-based construction instead: [`dist_from_quantiles`](@ref).
"""
function dist_from_mean_var(::Type{Cauchy}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(Cauchy, μ̄, σ̄²)
end

"""
    dist_from_mean_var(::Type{Logistic}, μ̄, σ̄²)

Direct formula. Construct a `Logistic(μ, s)` distribution. Any `μ̄ ∈ ℝ` and `σ̄² > 0`.
"""
function dist_from_mean_var(::Type{Logistic}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(Logistic, μ̄, σ̄²)
    x̄=μ̄
    s=√(3*σ̄²/π^2)
    return Logistic(x̄,s)
end

"""
    dist_from_mean_var(::Type{Laplace}, μ̄, σ̄²)

Direct formula. Construct a `Laplace(μ, b)` distribution. Any `μ̄ ∈ ℝ` and `σ̄² > 0`.
"""
function dist_from_mean_var(::Type{Laplace}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(Laplace, μ̄, σ̄²)
    x̄=μ̄
    b=√(σ̄²/2)
    return Laplace(x̄,b)
end

"""
    dist_from_mean_var(::Type{LogNormal}, μ̄, σ̄²)

Direct formula. Construct a `LogNormal(μ_log, σ_log)` distribution.
Requires `μ̄ > 0` and `σ̄² > 0`.
"""
function dist_from_mean_var(::Type{LogNormal}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(LogNormal, μ̄, σ̄²)
    σ=√(log(σ̄²/μ̄^2+1))
    x̄=log(μ̄^2/√(σ̄²+μ̄^2))
    return LogNormal(x̄,σ)
end

"""
    dist_from_mean_var(::Type{Chisq}, μ̄, σ̄²)

Direct formula. Construct a `Chisq(k)` distribution.
Requires `μ̄ ∈ ℕ` and `σ̄² = 2μ̄` (1 DOF).
"""
function dist_from_mean_var(::Type{Chisq}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(Chisq, μ̄, σ̄²)
    return Chisq(μ̄)
end

"""
    dist_from_mean_var(::Type{Exponential}, μ̄, σ̄²)

Direct formula. Construct an `Exponential(μ)` distribution.
Requires `μ̄ > 0` and `σ̄² = μ̄²` (1 DOF).
"""
function dist_from_mean_var(::Type{Exponential}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(Exponential, μ̄, σ̄²)
    return Exponential(μ̄)
end

"""
    dist_from_mean_var(::Type{Gamma}, μ̄, σ̄²)

Direct formula. Construct a `Gamma(α, θ)` distribution.
Requires `μ̄ > 0` and `σ̄² > 0`.
"""
function dist_from_mean_var(::Type{Gamma}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(Gamma, μ̄, σ̄²)
    α = μ̄^2/σ̄²
    θ = σ̄²/μ̄
    return Gamma(α,θ)
end

"""
    dist_from_mean_var(::Type{Erlang}, μ̄, σ̄²)

Direct formula. Construct an `Erlang(k, θ)` distribution (Gamma with integer shape).
Requires `μ̄ > 0` and `σ̄² > 0`. Shape `k` is rounded to the nearest integer.
"""
function dist_from_mean_var(::Type{Erlang}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(Erlang, μ̄, σ̄²)
    k = round(Int, μ̄^2/σ̄²)
    θ = σ̄²/μ̄
    return Erlang(k, θ)
end


function ron_ashri_evt_approximation(μ̄::Number, σ̄²::Number, positiveSolution::Bool)
    μ_b  = BigFloat(μ̄)
    var_b = BigFloat(σ̄²)
    CV = √(var_b)/μ_b
    f(x) = x/beta(1/x,1/x)-(1+CV^2)/2
    if 0 < CV^2 < 1
        lowerBound = positiveSolution>0 ? 1/CV : min(-√(2π), -1/CV);
        upperBound = positiveSolution>0 ? (CV^2+1)/(2CV^2) : -2(1+CV^2)/(CV^2);
        return find_zero(f, (lowerBound, upperBound));
    elseif CV^2 == 1
        return positiveSolution>0 ? 1 : find_zero(f, -√(7));
    elseif CV^2 > 1
        lowerBound = positiveSolution>0 ? 0 : -2;
        upperBound = positiveSolution>0 ? 1 : -2(1+CV^2)/(CV^2);
        while true
            tmp = (upperBound+lowerBound)/2
            if f(tmp)<0
                upperBound = tmp;
            elseif f(tmp)>0
                lowerBound = tmp;
                break
            else
                return tmp;
            end
        end
        return find_zero(f, (lowerBound, upperBound));
    end
end

function ron_ashri_weibull_approximation(μ̄::Number, σ̄²::Number)
    k = ron_ashri_evt_approximation(μ̄,σ̄², true)
    λ = μ̄/gamma(1+1/k)
    return Weibull(k, λ)
end

function ron_ashri_frechet_approximation(μ̄::Number, σ̄²::Number)
    α=-1*ron_ashri_evt_approximation(μ̄,σ̄², false)
    s = μ̄/gamma(1-1/α)
    return Frechet(α, s)
end

"""
    dist_from_mean_var(::Type{Frechet}, μ̄, σ̄²)

Numerical (root-finding). Construct a `Frechet(α, s)` distribution.
Requires `μ̄ > 0` and `σ̄² > 0`. Solves the beta-ratio equation.
"""
function dist_from_mean_var(::Type{Frechet}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(Frechet, μ̄, σ̄²)
    return ron_ashri_frechet_approximation(μ̄,σ̄²)
end

"""
    dist_from_mean_var(::Type{Weibull}, μ̄, σ̄²)

Numerical (root-finding). Construct a `Weibull(k, λ)` distribution.
Requires `μ̄ > 0` and `σ̄² > 0`. Solves the beta-ratio equation.
"""
function dist_from_mean_var(::Type{Weibull}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(Weibull, μ̄, σ̄²)
    return ron_ashri_weibull_approximation(μ̄,σ̄²)
end

"""
    dist_from_mean_var(::Type{Gumbel}, μ̄, σ̄²)

Direct formula. Construct a `Gumbel(μ_loc, β)` distribution.
Any `μ̄ ∈ ℝ` and `σ̄² > 0`.
"""
function dist_from_mean_var(::Type{Gumbel}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(Gumbel, μ̄, σ̄²)
    β = √(6*σ̄²/π^2)
    x̄ = μ̄-β*Base.MathConstants.γ
    return Gumbel(x̄,β)
end

"""
    dist_from_mean_var(::Type{Chi}, μ̄, σ̄²)

Direct formula. Construct a `Chi(ν)` distribution. Requires `μ̄ > 0`.
Degrees of freedom: `ν = μ̄² + σ̄²`.
"""
function dist_from_mean_var(::Type{Chi}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(Chi, μ̄, σ̄²)
    return Chi(μ̄^2+σ̄²)
end

"""
    dist_from_mean_var(::Type{Rayleigh}, μ̄, σ̄²)

Direct formula. Construct a `Rayleigh(σ)` distribution.
Requires `μ̄ > 0` and `CV = √((4-π)/π)` (1 DOF).
"""
function dist_from_mean_var(::Type{Rayleigh}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(Rayleigh, μ̄, σ̄²)
    σ=√(2/π)*μ̄
    return Rayleigh(σ)
end

"""
    dist_from_mean_var(::Type{FDist}, μ̄, σ̄²)

Direct formula. Construct an `FDist(ν₁, ν₂)` distribution.
Requires `1 < μ̄ < 2` and `σ̄² > μ̄²(μ̄-1)/(2-μ̄)`.
"""
function dist_from_mean_var(::Type{FDist}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(FDist, μ̄, σ̄²)
    v₂ = 2*μ̄/(μ̄-1)
    v₁ = 2*μ̄^2*(v₂-2)/(σ̄²*(v₂-4)-2*μ̄^2)
    return FDist(v₁,v₂)
end

"""
    dist_from_mean_var(::Type{InverseGamma}, μ̄, σ̄²)

Direct formula. Construct an `InverseGamma(α, β)` distribution.
Requires `μ̄ > 0` and `σ̄² > 0`.
"""
function dist_from_mean_var(::Type{InverseGamma}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(InverseGamma, μ̄, σ̄²)
    α=(μ̄^2+2*σ̄²)/σ̄²
    β=μ̄*(α-1)
    return InverseGamma(α,β)
end

"""
    dist_from_mean_var(::Type{Binomial}, μ̄, σ̄²)

Direct formula. Construct a `Binomial(n, p)` distribution.
Requires `μ̄ > 0` and `σ̄² < μ̄`. Parameter `n` is rounded to the nearest integer.
"""
function dist_from_mean_var(::Type{Binomial}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(Binomial, μ̄, σ̄²)
    p=1-σ̄²/μ̄
    n=round(Int, μ̄/p)
    return Binomial(n,p)
end

"""
    dist_from_mean_var(::Type{Poisson}, μ̄, σ̄²)

Direct formula. Construct a `Poisson(μ)` distribution.
Requires `σ̄² = μ̄` (1 DOF).
"""
function dist_from_mean_var(::Type{Poisson}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(Poisson, μ̄, σ̄²)
    return Poisson(μ̄)
end

"""
    dist_from_mean_var(::Type{NegativeBinomial}, μ̄, σ̄²)

Direct formula. Construct a `NegativeBinomial(r, p)` distribution.
Requires `σ̄² > μ̄ > 0`.
"""
function dist_from_mean_var(::Type{NegativeBinomial}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(NegativeBinomial, μ̄, σ̄²)
    p=μ̄/σ̄²
    r=μ̄^2/(σ̄²-μ̄)
    return NegativeBinomial(r,p)
end


# --- Recently implemented ---

"""
    dist_from_mean_var(::Type{Pareto}, μ̄, σ̄²)

Direct formula. Construct a `Pareto(α, θ)` distribution.
Requires `μ̄ > 0` and `σ̄² > 0`.
Shape `α` is derived from the coefficient of variation: `α = 1 + √(1 + 1/CV²)`.
"""
function dist_from_mean_var(::Type{Pareto}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(Pareto, μ̄, σ̄²)
    CV² = σ̄² / μ̄^2
    α = 1 + √(1 + 1 / CV²)
    θ = μ̄ * (α - 1) / α
    return Pareto(α, θ)
end

# TODO: needs further testing and validation
"""
    dist_from_mean_var(::Type{FoldedNormal}, μ̄, σ̄²)

Numerical (2D Newton iteration). Construct the parent `Normal(μp, σp)` whose folded
version `|X|` has mean `μ̄` and variance `σ̄²`. Requires `μ̄ > 0` and `σ̄² > 0`.
"""
function dist_from_mean_var(::Type{FoldedNormal}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(FoldedNormal, μ̄, σ̄²)
    μp, σp = _solve_folded_normal(Float64(μ̄), Float64(σ̄²))
    return Normal(μp, σp)  # returns the parent Normal; user takes |X|
end

"""
    dist_from_mean_var(::Type{Geometric}, μ̄, σ̄²)

Direct formula. Construct a `Geometric(p)` distribution.
Requires `μ̄ > 0` and `σ̄² = μ̄(1+μ̄)` (1 DOF).
"""
function dist_from_mean_var(::Type{Geometric}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(Geometric, μ̄, σ̄²)
    p = 1 / (1 + μ̄)
    return Geometric(p)
end

# TODO: truncated distributions need further testing and validation
"""
    dist_from_mean_var(d::Truncated{<:Normal}, μ̄, σ̄²)

Numerical (2D Newton iteration). Construct a truncated Normal on `[lo, hi]`
(taken from `d`) with mean `μ̄` and variance `σ̄²`. Uses quadrature for moments.
"""
function dist_from_mean_var(d::Truncated{<:Normal}, μ̄::Number, σ̄²::Number)
    lo, hi = extrema(d)
    return _solve_truncated_mean_var(Normal, lo, hi, Float64(μ̄), Float64(σ̄²))
end

"""
    dist_from_mean_var(d::Truncated{<:Laplace}, μ̄, σ̄²)

Numerical (2D Newton iteration). Construct a truncated Laplace on `[lo, hi]`
(taken from `d`) with mean `μ̄` and variance `σ̄²`. Uses quadrature for moments.
"""
function dist_from_mean_var(d::Truncated{<:Laplace}, μ̄::Number, σ̄²::Number)
    lo, hi = extrema(d)
    return _solve_truncated_mean_var(Laplace, lo, hi, Float64(μ̄), Float64(σ̄²))
end

"""
    dist_from_mean_var(d::Truncated{<:Logistic}, μ̄, σ̄²)

Numerical (2D Newton iteration). Construct a truncated Logistic on `[lo, hi]`
(taken from `d`) with mean `μ̄` and variance `σ̄²`. Uses quadrature for moments.
"""
function dist_from_mean_var(d::Truncated{<:Logistic}, μ̄::Number, σ̄²::Number)
    lo, hi = extrema(d)
    return _solve_truncated_mean_var(Logistic, lo, hi, Float64(μ̄), Float64(σ̄²))
end

function dist_from_mean_var(::Type{TriangularDist}, μ̄::Number, σ̄²::Number)
    throw(ErrorException("TriangularDist: dist_from_mean_var not yet implemented"))
end

"""
    dist_from_mean_var(::Type{SymTriangularDist}, μ̄, σ̄²)

Direct formula. Construct a `SymTriangularDist(μ, s)` distribution.
Any `μ̄ ∈ ℝ` and `σ̄² > 0`. Scale is `s = √(6 σ̄²)`.
"""
function dist_from_mean_var(::Type{SymTriangularDist}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(SymTriangularDist, μ̄, σ̄²)
    s = √(6 * σ̄²)
    return SymTriangularDist(μ̄, s)
end

function dist_from_mean_var(::Type{DiscreteTriangular}, μ̄::Number, σ̄²::Number)
    throw(ErrorException("DiscreteTriangular: dist_from_mean_var not yet implemented"))
end

function dist_from_mean_var(::Type{DiscreteSymmetricTriangular}, μ̄::Number, σ̄²::Number)
    throw(ErrorException("DiscreteSymmetricTriangular: dist_from_mean_var not yet implemented"))
end

function dist_from_mean_var(::Type{TruncatedPoisson}, μ̄::Number, σ̄²::Number)
    throw(ErrorException("TruncatedPoisson: dist_from_mean_var not yet implemented"))
end

"""
    dist_from_mean_var(::Type{DiscreteUniform}, μ̄, σ̄²)

Direct formula. Construct a `DiscreteUniform(a, b)` distribution. Requires that
`n = b - a` resolves to a non-negative integer and `a` is an integer.
"""
function dist_from_mean_var(::Type{DiscreteUniform}, μ̄::Number, σ̄²::Number)
    exists_unique_dist_from_mean_var(DiscreteUniform, μ̄, σ̄²)
    n = round(Int, -1 + √(1 + 12 * σ̄²))
    a = round(Int, μ̄ - n / 2)
    b = a + n
    return DiscreteUniform(a, b)
end
