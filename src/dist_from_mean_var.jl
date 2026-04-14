using Distributions
using Roots
using Polynomials
using SpecialFunctions


"""
    dist_from_mean_var(D, μ, var)

Construct a distribution of type `D` with the given mean `μ` and variance `var`.

Dispatches on the distribution type (or instance for truncated/TDist).
Throws `DomainError` if no valid distribution exists for the given moments.
Use `exists_unique_dist_from_mean_var` to check feasibility before calling.

See also: [`dist_from_mean_std`](@ref), [`dist_from_mean_cv`](@ref),
[`exists_unique_dist_from_mean_var`](@ref)
"""
function dist_from_mean_var end

"""
    dist_from_mean_var(::Type{Beta}, μ, var)

Construct a `Beta(α, β)` distribution. Requires `0 < μ < 1` and `0 < var < μ(1-μ)`.
"""
function dist_from_mean_var(::Type{Beta}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Beta, μ, var)
    S = (μ*(1-μ))/var-1
    α = μ*S
    β = (1-μ)*S
    return Beta(α,β)
end

"""
    dist_from_mean_var(::Type{Uniform}, μ, var)

Construct a `Uniform(a, b)` distribution. Any `μ ∈ ℝ` and `var > 0`.
"""
function dist_from_mean_var(::Type{Uniform}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Uniform, μ, var)
    diff = √(3*var)
    a = μ-diff
    b = μ+diff
    return Uniform(a,b)
end

"""
    dist_from_mean_var(::Type{Normal}, μ, var)

Construct a `Normal(μ, σ)` distribution. Any `μ ∈ ℝ` and `var > 0`.
"""
function dist_from_mean_var(::Type{Normal}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Normal, μ, var)
    return Normal(μ,√(var))
end

"""
    dist_from_mean_var(::Type{TDist}, μ, var)

Construct a standard `TDist(ν)` distribution. Requires `μ = 0` and `var > 1`.
"""
function dist_from_mean_var(::Type{TDist}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(TDist, μ, var)
    v=2*var/(var-1)
    return TDist(v)
end

"""
    dist_from_mean_var(d::TDist, μ, var)

Construct an affine-transformed `TDist` (location-scale) with arbitrary mean and variance.
The input `d` provides the degrees of freedom `ν` (must be > 2).
Returns `μ + σ * d` as a `LocationScale` distribution.
"""
function dist_from_mean_var(d::TDist, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(d, μ, var)
    ν = dof(d)
    base_var = ν / (ν - 2)
    σ = √(var / base_var)
    return μ + σ * d
end

"""
    dist_from_mean_var(::Type{Cauchy}, μ, var)

Always throws `DomainError` — the Cauchy distribution has no defined mean or variance.
Use quantile-based construction instead: [`dist_from_quantiles`](@ref).
"""
function dist_from_mean_var(::Type{Cauchy}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Cauchy, μ, var)
end

"""
    dist_from_mean_var(::Type{Logistic}, μ, var)

Construct a `Logistic(μ, s)` distribution. Any `μ ∈ ℝ` and `var > 0`.
"""
function dist_from_mean_var(::Type{Logistic}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Logistic, μ, var)
    x̄=μ
    s=√(3*var/π^2)
    return Logistic(x̄,s)
end

"""
    dist_from_mean_var(::Type{Laplace}, μ, var)

Construct a `Laplace(μ, b)` distribution. Any `μ ∈ ℝ` and `var > 0`.
"""
function dist_from_mean_var(::Type{Laplace}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Laplace, μ, var)
    x̄=μ
    b=√(var/2)
    return Laplace(x̄,b)
end

"""
    dist_from_mean_var(::Type{LogNormal}, μ, var)

Construct a `LogNormal(μ_log, σ_log)` distribution. Requires `μ > 0` and `var > 0`.
"""
function dist_from_mean_var(::Type{LogNormal}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(LogNormal, μ, var)
    σ=√(log(var/μ^2+1))
    x̄=log(μ^2/√(var+μ^2))
    return LogNormal(x̄,σ)
end

"""
    dist_from_mean_var(::Type{Chisq}, μ, var)

Construct a `Chisq(k)` distribution. Requires `μ ∈ ℕ` and `var = 2μ` (1 DOF).
"""
function dist_from_mean_var(::Type{Chisq}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Chisq, μ, var)
    return Chisq(μ)
end

"""
    dist_from_mean_var(::Type{Exponential}, μ, var)

Construct an `Exponential(μ)` distribution. Requires `μ > 0` and `var = μ²` (1 DOF).
"""
function dist_from_mean_var(::Type{Exponential}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Exponential, μ, var)
    return Exponential(μ)
end

"""
    dist_from_mean_var(::Type{Gamma}, μ, var)

Construct a `Gamma(α, θ)` distribution. Requires `μ > 0` and `var > 0`.
"""
function dist_from_mean_var(::Type{Gamma}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Gamma, μ, var)
    α = μ^2/var
    θ = var/μ
    return Gamma(α,θ)
end

"""
    dist_from_mean_var(::Type{Erlang}, μ, var)

Construct an `Erlang(k, θ)` distribution (Gamma with integer shape).
Requires `μ > 0` and `var > 0`. Shape `k` is rounded to the nearest integer.
"""
function dist_from_mean_var(::Type{Erlang}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Erlang, μ, var)
    k = round(Int, μ^2/var)
    θ = var/μ
    return Erlang(k, θ)
end


function ron_ashri_evt_approximation(μ::Number, var::Number, positiveSolution::Bool)
    μ_b  = BigFloat(μ)
    var_b = BigFloat(var)
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

function ron_ashri_weibull_approximation(μ::Number, var::Number)
    k = ron_ashri_evt_approximation(μ,var, true)
    λ = μ/gamma(1+1/k)
    return Weibull(k, λ)
end

function ron_ashri_frechet_approximation(μ::Number, var::Number)
    α=-1*ron_ashri_evt_approximation(μ,var, false)
    s = μ/gamma(1-1/α)
    return Frechet(α, s)
end

"""
    dist_from_mean_var(::Type{Frechet}, μ, var)

Construct a `Frechet(α, s)` distribution. Requires `μ > 0` and `var > 0`.
Uses the Ron Ashri EVT approximation via the beta-ratio equation.
"""
function dist_from_mean_var(::Type{Frechet}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Frechet, μ, var)
    return ron_ashri_frechet_approximation(μ,var)
end

"""
    dist_from_mean_var(::Type{Weibull}, μ, var)

Construct a `Weibull(k, λ)` distribution. Requires `μ > 0` and `var > 0`.
Uses the Ron Ashri EVT approximation via the beta-ratio equation.
"""
function dist_from_mean_var(::Type{Weibull}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Weibull, μ, var)
    return ron_ashri_weibull_approximation(μ,var)
end

"""
    dist_from_mean_var(::Type{Gumbel}, μ, var)

Construct a `Gumbel(μ_loc, β)` distribution. Any `μ ∈ ℝ` and `var > 0`.
"""
function dist_from_mean_var(::Type{Gumbel}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Gumbel, μ, var)
    β = √(6*var/π^2)
    x̄ = μ-β*Base.MathConstants.γ
    return Gumbel(x̄,β)
end

"""
    dist_from_mean_var(::Type{Chi}, μ, var)

Construct a `Chi(ν)` distribution. Requires `μ > 0`. The degrees of freedom
are computed as `ν = μ² + var`.
"""
function dist_from_mean_var(::Type{Chi}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Chi, μ, var)
    return Chi(μ^2+var)
end

"""
    dist_from_mean_var(::Type{Rayleigh}, μ, var)

Construct a `Rayleigh(σ)` distribution. Requires `μ > 0` and `CV = √((4-π)/π)` (1 DOF).
"""
function dist_from_mean_var(::Type{Rayleigh}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Rayleigh, μ, var)
    σ=√(2/π)*μ
    return Rayleigh(σ)
end

"""
    dist_from_mean_var(::Type{FDist}, μ, var)

Construct an `FDist(ν₁, ν₂)` distribution. Requires `1 < μ < 2` and
`var > μ²(μ-1)/(2-μ)`.
"""
function dist_from_mean_var(::Type{FDist}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(FDist, μ, var)
    v₂ = 2*μ/(μ-1)
    v₁ = 2*μ^2*(v₂-2)/(var*(v₂-4)-2*μ^2)
    return FDist(v₁,v₂)
end

"""
    dist_from_mean_var(::Type{InverseGamma}, μ, var)

Construct an `InverseGamma(α, β)` distribution. Requires `μ > 0` and `var > 0`.
"""
function dist_from_mean_var(::Type{InverseGamma}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(InverseGamma, μ, var)
    α=(μ^2+2*var)/var
    β=μ*(α-1)
    return InverseGamma(α,β)
end

"""
    dist_from_mean_var(::Type{Binomial}, μ, var)

Construct a `Binomial(n, p)` distribution. Requires `μ > 0` and `var < μ`.
Parameter `n` is rounded to the nearest integer.
"""
function dist_from_mean_var(::Type{Binomial}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Binomial, μ, var)
    p=1-var/μ
    n=round(Int, μ/p)
    return Binomial(n,p)
end

"""
    dist_from_mean_var(::Type{Poisson}, μ, var)

Construct a `Poisson(μ)` distribution. Requires `var = μ` (1 DOF).
"""
function dist_from_mean_var(::Type{Poisson}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Poisson, μ, var)
    return Poisson(μ)
end

"""
    dist_from_mean_var(::Type{NegativeBinomial}, μ, var)

Construct a `NegativeBinomial(r, p)` distribution. Requires `var > μ > 0`.
"""
function dist_from_mean_var(::Type{NegativeBinomial}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(NegativeBinomial, μ, var)
    p=μ/var
    r=μ^2/(var-μ)
    return NegativeBinomial(r,p)
end


# --- Recently implemented ---

"""
    dist_from_mean_var(::Type{Pareto}, μ, var)

Construct a `Pareto(α, θ)` distribution. Requires `μ > 0` and `var > 0`.
Shape `α` is derived from the coefficient of variation: `α = 1 + √(1 + 1/CV²)`.
"""
function dist_from_mean_var(::Type{Pareto}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Pareto, μ, var)
    CV² = var / μ^2
    α = 1 + √(1 + 1 / CV²)
    θ = μ * (α - 1) / α
    return Pareto(α, θ)
end

# TODO: needs further testing and validation
"""
    dist_from_mean_var(::Type{FoldedNormal}, μ, var)

Construct the parent `Normal(μp, σp)` whose folded version `|X|` has the given
mean and variance. Uses 2D Newton iteration. Requires `μ > 0` and `var > 0`.
"""
function dist_from_mean_var(::Type{FoldedNormal}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(FoldedNormal, μ, var)
    μp, σp = _solve_folded_normal(Float64(μ), Float64(var))
    return Normal(μp, σp)  # returns the parent Normal; user takes |X|
end

"""
    dist_from_mean_var(::Type{Geometric}, μ, var)

Construct a `Geometric(p)` distribution. Requires `μ > 0` and `var = μ(1+μ)` (1 DOF).
"""
function dist_from_mean_var(::Type{Geometric}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Geometric, μ, var)
    p = 1 / (1 + μ)
    return Geometric(p)
end

# TODO: truncated distributions need further testing and validation
"""
    dist_from_mean_var(d::Truncated{<:Normal}, μ, var)

Construct a truncated Normal on `[lo, hi]` (taken from `d`) with the given mean
and variance. Uses 2D Newton iteration with quadrature-based moment computation.
"""
function dist_from_mean_var(d::Truncated{<:Normal}, μ::Number, var::Number)
    lo, hi = extrema(d)
    return _solve_truncated_mean_var(Normal, lo, hi, Float64(μ), Float64(var))
end

"""
    dist_from_mean_var(d::Truncated{<:Laplace}, μ, var)

Construct a truncated Laplace on `[lo, hi]` (taken from `d`) with the given mean
and variance. Uses 2D Newton iteration with quadrature-based moment computation.
"""
function dist_from_mean_var(d::Truncated{<:Laplace}, μ::Number, var::Number)
    lo, hi = extrema(d)
    return _solve_truncated_mean_var(Laplace, lo, hi, Float64(μ), Float64(var))
end

"""
    dist_from_mean_var(d::Truncated{<:Logistic}, μ, var)

Construct a truncated Logistic on `[lo, hi]` (taken from `d`) with the given mean
and variance. Uses 2D Newton iteration with quadrature-based moment computation.
"""
function dist_from_mean_var(d::Truncated{<:Logistic}, μ::Number, var::Number)
    lo, hi = extrema(d)
    return _solve_truncated_mean_var(Logistic, lo, hi, Float64(μ), Float64(var))
end

function dist_from_mean_var(::Type{TriangularDist}, μ::Number, var::Number)
    throw(ErrorException("TriangularDist: dist_from_mean_var not yet implemented"))
end

"""
    dist_from_mean_var(::Type{SymTriangularDist}, μ, var)

Construct a `SymTriangularDist(μ, s)` distribution. Any `μ ∈ ℝ` and `var > 0`.
Scale is `s = √(6 var)`.
"""
function dist_from_mean_var(::Type{SymTriangularDist}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(SymTriangularDist, μ, var)
    s = √(6 * var)
    return SymTriangularDist(μ, s)
end

function dist_from_mean_var(::Type{DiscreteTriangular}, μ::Number, var::Number)
    throw(ErrorException("DiscreteTriangular: dist_from_mean_var not yet implemented"))
end

function dist_from_mean_var(::Type{DiscreteSymmetricTriangular}, μ::Number, var::Number)
    throw(ErrorException("DiscreteSymmetricTriangular: dist_from_mean_var not yet implemented"))
end

function dist_from_mean_var(::Type{TruncatedPoisson}, μ::Number, var::Number)
    throw(ErrorException("TruncatedPoisson: dist_from_mean_var not yet implemented"))
end

"""
    dist_from_mean_var(::Type{DiscreteUniform}, μ, var)

Construct a `DiscreteUniform(a, b)` distribution. Requires that `n = b - a`
resolves to a non-negative integer and `a` is an integer.
"""
function dist_from_mean_var(::Type{DiscreteUniform}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(DiscreteUniform, μ, var)
    n = round(Int, -1 + √(1 + 12 * var))
    a = round(Int, μ - n / 2)
    b = a + n
    return DiscreteUniform(a, b)
end
