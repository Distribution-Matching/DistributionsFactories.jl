# `exists_dist_from_mean_var` methods live in exists_dist.jl (included earlier
# in DistributionsFactories.jl); this file owns the constructors below.


"""
    dist_from_mean_var(D, μ̄, σ̄²)

Construct a distribution of type `D` with the given mean `μ̄` and variance `σ̄²`.

Dispatches on the distribution type (or instance for truncated/TDist).
Throws `DomainError` with a reason when no valid distribution exists for the
given moments. Use [`exists_dist_from_mean_var`](@ref) for a non-throwing
`Bool` feasibility predicate.

See also: [`make_dist`](@ref), [`exists_dist_from_mean_var`](@ref).
"""
function dist_from_mean_var end

"""
    dist_from_mean_var(::Type{Beta}, μ̄, σ̄²)

Direct formula. Construct a `Beta(α, β)` distribution.
Requires `0 < μ̄ < 1` and `0 < σ̄² < μ̄(1-μ̄)`.
"""
function dist_from_mean_var(::Type{Beta}, μ̄::Number, σ̄²::Number)
    _require_dist_from_mean_var(Beta, μ̄, σ̄²)
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
    _require_dist_from_mean_var(Uniform, μ̄, σ̄²)
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
    _require_dist_from_mean_var(Normal, μ̄, σ̄²)
    return Normal(μ̄,√(σ̄²))
end

"""
    dist_from_mean_var(::Type{TDist}, μ̄, σ̄²)

Direct formula. Construct a standard `TDist(ν)` distribution.
Requires `μ̄ = 0` and `σ̄² > 1`.
"""
function dist_from_mean_var(::Type{TDist}, μ̄::Number, σ̄²::Number)
    _require_dist_from_mean_var(TDist, μ̄, σ̄²)
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
    _require_dist_from_mean_var(d, μ̄, σ̄²)
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
    _require_dist_from_mean_var(Cauchy, μ̄, σ̄²)
end

"""
    dist_from_mean_var(::Type{Logistic}, μ̄, σ̄²)

Direct formula. Construct a `Logistic(μ, s)` distribution. Any `μ̄ ∈ ℝ` and `σ̄² > 0`.
"""
function dist_from_mean_var(::Type{Logistic}, μ̄::Number, σ̄²::Number)
    _require_dist_from_mean_var(Logistic, μ̄, σ̄²)
    x̄=μ̄
    s=√(3*σ̄²/π^2)
    return Logistic(x̄,s)
end

"""
    dist_from_mean_var(::Type{Laplace}, μ̄, σ̄²)

Direct formula. Construct a `Laplace(μ, b)` distribution. Any `μ̄ ∈ ℝ` and `σ̄² > 0`.
"""
function dist_from_mean_var(::Type{Laplace}, μ̄::Number, σ̄²::Number)
    _require_dist_from_mean_var(Laplace, μ̄, σ̄²)
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
    _require_dist_from_mean_var(LogNormal, μ̄, σ̄²)
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
    _require_dist_from_mean_var(Chisq, μ̄, σ̄²)
    return Chisq(μ̄)
end

"""
    dist_from_mean_var(::Type{Exponential}, μ̄, σ̄²)

Direct formula. Construct an `Exponential(μ)` distribution.
Requires `μ̄ > 0` and `σ̄² = μ̄²` (1 DOF).
"""
function dist_from_mean_var(::Type{Exponential}, μ̄::Number, σ̄²::Number)
    _require_dist_from_mean_var(Exponential, μ̄, σ̄²)
    return Exponential(μ̄)
end

"""
    dist_from_mean_var(::Type{Gamma}, μ̄, σ̄²)

Direct formula. Construct a `Gamma(α, θ)` distribution.
Requires `μ̄ > 0` and `σ̄² > 0`.
"""
function dist_from_mean_var(::Type{Gamma}, μ̄::Number, σ̄²::Number)
    _require_dist_from_mean_var(Gamma, μ̄, σ̄²)
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
    _require_dist_from_mean_var(Erlang, μ̄, σ̄²)
    k = round(Int, μ̄^2/σ̄²)
    θ = σ̄²/μ̄
    return Erlang(k, θ)
end


function _solve_evt_shape(μ̄::Number, σ̄²::Number, positiveSolution::Bool)
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

function _weibull_from_mean_var(μ̄::Number, σ̄²::Number)
    k = _solve_evt_shape(μ̄,σ̄², true)
    λ = μ̄/gamma(1+1/k)
    return Weibull(k, λ)
end

function _frechet_from_mean_var(μ̄::Number, σ̄²::Number)
    α=-1*_solve_evt_shape(μ̄,σ̄², false)
    s = μ̄/gamma(1-1/α)
    return Frechet(α, s)
end

"""
    dist_from_mean_var(::Type{Frechet}, μ̄, σ̄²)

Numerical (root-finding). Construct a `Frechet(α, s)` distribution.
Requires `μ̄ > 0` and `σ̄² > 0`. Solves the beta-ratio equation.
"""
function dist_from_mean_var(::Type{Frechet}, μ̄::Number, σ̄²::Number)
    _require_dist_from_mean_var(Frechet, μ̄, σ̄²)
    return _frechet_from_mean_var(μ̄,σ̄²)
end

"""
    dist_from_mean_var(::Type{Weibull}, μ̄, σ̄²)

Numerical (root-finding). Construct a `Weibull(k, λ)` distribution.
Requires `μ̄ > 0` and `σ̄² > 0`. Solves the beta-ratio equation.
"""
function dist_from_mean_var(::Type{Weibull}, μ̄::Number, σ̄²::Number)
    _require_dist_from_mean_var(Weibull, μ̄, σ̄²)
    return _weibull_from_mean_var(μ̄,σ̄²)
end

"""
    dist_from_mean_var(::Type{Gumbel}, μ̄, σ̄²)

Direct formula. Construct a `Gumbel(μ_loc, β)` distribution.
Any `μ̄ ∈ ℝ` and `σ̄² > 0`.
"""
function dist_from_mean_var(::Type{Gumbel}, μ̄::Number, σ̄²::Number)
    _require_dist_from_mean_var(Gumbel, μ̄, σ̄²)
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
    _require_dist_from_mean_var(Chi, μ̄, σ̄²)
    return Chi(μ̄^2+σ̄²)
end

"""
    dist_from_mean_var(::Type{Rayleigh}, μ̄, σ̄²)

Direct formula. Construct a `Rayleigh(σ)` distribution.
Requires `μ̄ > 0` and `CV = √((4-π)/π)` (1 DOF).
"""
function dist_from_mean_var(::Type{Rayleigh}, μ̄::Number, σ̄²::Number)
    _require_dist_from_mean_var(Rayleigh, μ̄, σ̄²)
    σ=√(2/π)*μ̄
    return Rayleigh(σ)
end

"""
    dist_from_mean_var(::Type{FDist}, μ̄, σ̄²)

Direct formula. Construct an `FDist(ν₁, ν₂)` distribution.
Requires `1 < μ̄ < 2` and `σ̄² > μ̄²(μ̄-1)/(2-μ̄)`.
"""
function dist_from_mean_var(::Type{FDist}, μ̄::Number, σ̄²::Number)
    _require_dist_from_mean_var(FDist, μ̄, σ̄²)
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
    _require_dist_from_mean_var(InverseGamma, μ̄, σ̄²)
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
    _require_dist_from_mean_var(Binomial, μ̄, σ̄²)
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
    _require_dist_from_mean_var(Poisson, μ̄, σ̄²)
    return Poisson(μ̄)
end

"""
    dist_from_mean_var(::Type{NegativeBinomial}, μ̄, σ̄²)

Direct formula. Construct a `NegativeBinomial(r, p)` distribution.
Requires `σ̄² > μ̄ > 0`.
"""
function dist_from_mean_var(::Type{NegativeBinomial}, μ̄::Number, σ̄²::Number)
    _require_dist_from_mean_var(NegativeBinomial, μ̄, σ̄²)
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
    _require_dist_from_mean_var(Pareto, μ̄, σ̄²)
    CV² = σ̄² / μ̄^2
    α = 1 + √(1 + 1 / CV²)
    θ = μ̄ * (α - 1) / α
    return Pareto(α, θ)
end

"""
    dist_from_mean_var(::Type{FoldedNormal}, μ̄, σ̄²)

Numerical (2D Newton iteration). Construct a `FoldedNormal(μ, σ)` distribution
whose mean is `μ̄` and variance is `σ̄²`. Requires `μ̄ > 0` and `σ̄² > 0`.
"""
function dist_from_mean_var(::Type{FoldedNormal}, μ̄::Number, σ̄²::Number)
    _require_dist_from_mean_var(FoldedNormal, μ̄, σ̄²)
    μp, σp = _solve_folded_normal(Float64(μ̄), Float64(σ̄²))
    return FoldedNormal(μp, σp)
end

"""
    dist_from_mean_var(::Type{Geometric}, μ̄, σ̄²)

Direct formula. Construct a `Geometric(p)` distribution.
Requires `μ̄ > 0` and `σ̄² = μ̄(1+μ̄)` (1 DOF).
"""
function dist_from_mean_var(::Type{Geometric}, μ̄::Number, σ̄²::Number)
    _require_dist_from_mean_var(Geometric, μ̄, σ̄²)
    p = 1 / (1 + μ̄)
    return Geometric(p)
end

# Routes the truncated location-scale factory to the right solver based on
# whether bounds are finite. Two-sided uses the standardize-to-[-0.5, 0.5]
# unit solver; one-sided uses the standardize-to-[0, ∞) half-truncated
# solver. See `solvers.jl` for the canonical solver kernels.
function _dispatch_truncated_locscale_factory(::Type{D}, lo::Real, hi::Real,
                                              μ̄::Number, σ̄²::Number) where {D<:Distribution}
    if isfinite(lo) && isfinite(hi)
        return _solve_truncated_mean_var(D, lo, hi, Float64(μ̄), Float64(σ̄²))
    elseif isfinite(lo) && !isfinite(hi)
        return _solve_truncated_half_below(D, lo, Float64(μ̄), Float64(σ̄²))
    elseif !isfinite(lo) && isfinite(hi)
        return _solve_truncated_half_above(D, hi, Float64(μ̄), Float64(σ̄²))
    else
        # Both infinite: no truncation; fall back to the untruncated factory.
        return dist_from_mean_var(D, μ̄, σ̄²)
    end
end

"""
    dist_from_mean_var(d::Truncated{<:Normal}, μ̄, σ̄²)

Numerical. Construct a truncated Normal on `extrema(d)` with the given
moments. Two-sided uses the standardize-to-`[-0.5, 0.5]` 2D Newton solver;
one-sided uses the standardize-to-`[0, ∞)` half-truncated solver.
"""
function dist_from_mean_var(d::Truncated{<:Normal}, μ̄::Number, σ̄²::Number)
    _require_dist_from_mean_var(d, μ̄, σ̄²)
    lo, hi = extrema(d)
    return _dispatch_truncated_locscale_factory(Normal, lo, hi, μ̄, σ̄²)
end

"""
    dist_from_mean_var(d::Truncated{<:Laplace}, μ̄, σ̄²)

Numerical. As for `Truncated{<:Normal}`.
"""
function dist_from_mean_var(d::Truncated{<:Laplace}, μ̄::Number, σ̄²::Number)
    _require_dist_from_mean_var(d, μ̄, σ̄²)
    lo, hi = extrema(d)
    return _dispatch_truncated_locscale_factory(Laplace, lo, hi, μ̄, σ̄²)
end

"""
    dist_from_mean_var(d::Truncated{<:Logistic}, μ̄, σ̄²)

Numerical. As for `Truncated{<:Normal}`.
"""
function dist_from_mean_var(d::Truncated{<:Logistic}, μ̄::Number, σ̄²::Number)
    _require_dist_from_mean_var(d, μ̄, σ̄²)
    lo, hi = extrema(d)
    return _dispatch_truncated_locscale_factory(Logistic, lo, hi, μ̄, σ̄²)
end

function dist_from_mean_var(::Type{TriangularDist}, μ̄::Number, σ̄²::Number)
    throw(ArgumentError("TriangularDist: 3 parameters and only 2 moment constraints — supply `mode` as well (use `make_dist(TriangularDist, mean=…, var=…, mode=…)`)"))
end

"""
    dist_from_mean_var(::Type{SymTriangularDist}, μ̄, σ̄²)

Direct formula. Construct a `SymTriangularDist(μ, s)` distribution.
Any `μ̄ ∈ ℝ` and `σ̄² > 0`. Scale is `s = √(6 σ̄²)`.
"""
function dist_from_mean_var(::Type{SymTriangularDist}, μ̄::Number, σ̄²::Number)
    _require_dist_from_mean_var(SymTriangularDist, μ̄, σ̄²)
    s = √(6 * σ̄²)
    return SymTriangularDist(μ̄, s)
end

"""
    dist_from_mean_var(::Type{DiscreteTriangular}, μ̄, σ̄²)

Always throws: `DiscreteTriangular` has 3 integer parameters and mean+var
alone is underdetermined. Supply `mode` as well via
`make_dist(DiscreteTriangular, mean=…, var=…, mode=…)`.
"""
function dist_from_mean_var(::Type{DiscreteTriangular}, μ̄::Number, σ̄²::Number)
    throw(ArgumentError("DiscreteTriangular: mean+var alone is underdetermined; supply `mode` as well"))
end

"""
    dist_from_mean_var(::Type{DiscreteSymmetricTriangular}, μ̄, σ̄²)

Direct formula. Construct a `DiscreteSymmetricTriangular(μ, n)` distribution.
Requires `μ̄ ∈ ℤ` and `σ̄²` such that `n = -1 + √(1 + 6σ̄²)` is a non-negative
integer.
"""
function dist_from_mean_var(::Type{DiscreteSymmetricTriangular}, μ̄::Number, σ̄²::Number)
    _require_dist_from_mean_var(DiscreteSymmetricTriangular, μ̄, σ̄²)
    n = round(Int, -1 + √(1 + 6σ̄²))
    return DiscreteSymmetricTriangular(round(Int, μ̄), n)
end

"""
    dist_from_mean_var(d::Truncated{<:Poisson}, μ̄, σ̄²)

Numerical (1D root-finding). Construct a Poisson truncated to `extrema(d)` whose
mean equals `μ̄`. Variance is then determined; if `σ̄²` deviates significantly
from the resulting variance an `ArgumentError` is raised.
"""
function dist_from_mean_var(d::Truncated{<:Poisson}, μ̄::Number, σ̄²::Number)
    lo, hi = Float64(d.lower), Float64(d.upper)
    σ̄² > 0 || throw(DomainError(σ̄², "Poisson: σ̄² must be > 0"))
    (lo < μ̄ < hi) || throw(DomainError(μ̄, "Truncated Poisson: μ̄ must be in ($lo, $hi)"))
    td = _solve_truncated_poisson_mean(lo, hi, Float64(μ̄))
    achieved_var = _truncated_poisson_moments(td.untruncated.λ, lo, hi)[2]
    isapprox(achieved_var, σ̄²; rtol=1e-3) || throw(ArgumentError(
        "Truncated Poisson: variance is determined by mean on [$lo, $hi]. " *
        "Got achieved var $achieved_var, requested $σ̄²."))
    return td
end

"""
    dist_from_mean(d::Truncated{<:Poisson}, μ̄)

Numerical (1D root-finding). Construct a Poisson truncated to `extrema(d)` whose
mean equals `μ̄`. The Poisson rate λ is the single free parameter.
"""
function dist_from_mean(d::Truncated{<:Poisson}, μ̄::Number)
    lo, hi = Float64(d.lower), Float64(d.upper)
    (lo < μ̄ < hi) || throw(DomainError(μ̄, "Truncated Poisson: μ̄ must be in ($lo, $hi)"))
    return _solve_truncated_poisson_mean(lo, hi, Float64(μ̄))
end

"""
    dist_from_mean_var(::Type{DiscreteUniform}, μ̄, σ̄²)

Direct formula. Construct a `DiscreteUniform(a, b)` distribution. Requires that
`n = b - a` resolves to a non-negative integer and `a` is an integer.
"""
function dist_from_mean_var(::Type{DiscreteUniform}, μ̄::Number, σ̄²::Number)
    _require_dist_from_mean_var(DiscreteUniform, μ̄, σ̄²)
    n = round(Int, -1 + √(1 + 12 * σ̄²))
    a = round(Int, μ̄ - n / 2)
    b = a + n
    return DiscreteUniform(a, b)
end



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



# --- Variance-only construction (1-parameter distributions) ---

dist_from_var(D, σ̄²::Number) = dist_from_mean_var(D, _mean_from_var(D, σ̄²), σ̄²)

_mean_from_var(::Type{Exponential}, σ̄²::Number) = √σ̄²                         # σ̄² = μ̄²
_mean_from_var(::Type{Poisson}, σ̄²::Number) = σ̄²                               # σ̄² = μ̄
_mean_from_var(::Type{Chisq}, σ̄²::Number) = σ̄² / 2                             # σ̄² = 2μ̄
_mean_from_var(::Type{Rayleigh}, σ̄²::Number) = √(σ̄² * π / (4 - π))            # σ̄² = μ̄²(4-π)/π
_mean_from_var(::Type{Geometric}, σ̄²::Number) = (-1 + √(1 + 4σ̄²)) / 2        # σ̄² = μ̄(1+μ̄)
_mean_from_var(::Type{D}, σ̄²::Number) where {D<:Distribution} =
    throw(ErrorException("$D: dist_from_var not supported (mean is not determined by variance alone)"))

# --- Mode-based construction ---

function dist_from_mode end

# 1-parameter: mode determines everything
function dist_from_mode(::Type{Rayleigh}, m::Number)
    m > 0 || throw(DomainError(m, "Rayleigh: mode must be > 0"))
    return Rayleigh(m)  # mode = σ
end

# 2-parameter: mode + mean
function dist_from_mean_mode end

function dist_from_mean_mode(::Type{Gamma}, μ̄::Number, m::Number)
    μ̄ > 0 || throw(DomainError(μ̄, "Gamma: μ̄ must be > 0"))
    μ̄ > m || throw(DomainError("Gamma: mean must be > mode (mean=$μ̄, mode=$m)"))
    # mode = (α-1)θ, mean = αθ → θ = mean - mode, α = mean/θ
    θ = μ̄ - m
    α = μ̄ / θ
    return Gamma(α, θ)
end

function dist_from_mean_mode(::Type{Normal}, μ̄::Number, m::Number)
    # Normal: mode = mean, so this is only consistent if mode ≈ mean
    isapprox(μ̄, m, atol=1e-10) || throw(DomainError(
        "Normal: mode must equal mean (mean=$μ̄, mode=$m)"))
    throw(ArgumentError("Normal: mode=mean, need another constraint (var, std, etc.)"))
end

function dist_from_mean_mode(::Type{Beta}, μ̄::Number, m::Number)
    (0 < μ̄ < 1) || throw(DomainError(μ̄, "Beta: μ̄ must be in (0,1)"))
    (0 < m < 1) || throw(DomainError(m, "Beta: mode must be in (0,1)"))
    # mean = α/(α+β), mode = (α-1)/(α+β-2)
    # From these two equations:
    # α+β = α/μ̄  →  β = α(1-μ̄)/μ̄
    # α+β-2 = (α-1)/m  →  α/μ̄ - 2 = (α-1)/m
    # → α(1/μ̄ - 1/m) = 2 - 1/m  →  α = (2 - 1/m) / (1/μ̄ - 1/m)
    # → α = (2m - 1) / (m/μ̄ - 1) = μ̄(2m - 1) / (m - μ̄)
    (m != μ̄) || throw(DomainError("Beta: mode must differ from mean"))
    α = μ̄ * (2m - 1) / (m - μ̄)
    α > 1 || throw(DomainError("Beta: resulting α=$α must be > 1 for mode to exist"))
    β = α * (1 - μ̄) / μ̄
    β > 1 || throw(DomainError("Beta: resulting β=$β must be > 1 for mode to exist"))
    return Beta(α, β)
end

# 2-parameter: mode + var
function dist_from_mode_var end

function dist_from_mode_iqr end

function dist_from_mode_iqr(::Type{Normal}, m::Number, iqr::Number)
    # Normal: mode = μ, IQR = 2·z₀.₇₅·σ
    σ = iqr / (2 * quantile(Normal(), 0.75))
    return Normal(m, σ)
end

function dist_from_mode_iqr(::Type{Gamma}, m::Number, iqr::Number)
    m >= 0 || throw(DomainError(m, "Gamma: mode must be ≥ 0"))
    iqr > 0 || throw(DomainError(iqr, "Gamma: IQR must be > 0"))
    # mode = (α-1)θ, IQR = q75 - q25. Solve numerically for α.
    sol = find_zero(
        logα -> begin
            α = exp(logα) + 1  # ensure α > 1
            θ = m / (α - 1)
            d = Gamma(α, θ)
            quantile(d, 0.75) - quantile(d, 0.25) - iqr
        end,
        0.0
    )
    α = exp(sol) + 1
    θ = m / (α - 1)
    return Gamma(α, θ)
end

function dist_from_mode_iqr(::Type{Logistic}, m::Number, iqr::Number)
    # Logistic: mode = μ, IQR = 2·θ·ln(3)
    θ = iqr / (2 * log(3))
    return Logistic(m, θ)
end

function dist_from_mode_iqr(::Type{Laplace}, m::Number, iqr::Number)
    # Laplace: mode = μ, IQR = 2·b·ln(2)
    b = iqr / (2 * log(2))
    return Laplace(m, b)
end

function dist_from_mode_var(::Type{Normal}, m::Number, σ̄²::Number)
    # mode = μ for Normal
    return Normal(m, √σ̄²)
end

function dist_from_mode_quantile end

function dist_from_mode_quantile(::Type{Gamma}, m::Number, p::Number, q::Number)
    m >= 0 || throw(DomainError(m, "Gamma: mode must be ≥ 0"))
    (0 < p < 1) || throw(DomainError(p, "p must be in (0,1)"))
    q > 0 || throw(DomainError(q, "Gamma: quantile must be > 0"))
    # mode = (α-1)θ, quantile(Gamma(α,θ), p) = q
    # Solve for α: given α, θ = mode/(α-1), then check quantile
    sol = find_zero(
        logα -> begin
            α = exp(logα) + 1  # ensure α > 1
            θ = m / (α - 1)
            quantile(Gamma(α, θ), p) - q
        end,
        0.0  # initial guess: α ≈ 2
    )
    α = exp(sol) + 1
    θ = m / (α - 1)
    return Gamma(α, θ)
end

# --- Mean + variance + mode (3-parameter triangular families) ---

"""
    dist_from_mean_var_mode(D, μ̄, σ̄², mode)

Construct a 3-parameter distribution from a mean, variance, and mode. Used
by [`make_dist`](@ref) when `mean`, `var`, and `mode` are all supplied.
Currently implemented for `TriangularDist` (continuous, exact) and
`DiscreteTriangular` (integer, approximate).
"""
function dist_from_mean_var_mode end

"""
    dist_from_mean_var_mode(::Type{TriangularDist}, μ̄, σ̄², c)

Direct formula. Construct a `TriangularDist(a, b, c)` whose mean is `μ̄`,
variance `σ̄²`, and mode `c`. Solving:

    a + b = 3μ̄ - c
    ab    = ((3μ̄ - c)² + c² - c(3μ̄ - c) - 18σ̄²) / 3

with `a, b` the two roots of `t² - (a+b)t + ab = 0`. The discriminant must be
non-negative and the resulting `a ≤ c ≤ b`.
"""
function dist_from_mean_var_mode(::Type{TriangularDist}, μ̄::Number, σ̄²::Number, c::Number)
    σ̄² > 0 || throw(DomainError(σ̄², "TriangularDist: σ̄² must be > 0"))
    S = 3μ̄ - c
    ab = (S^2 + c^2 - c * S - 18σ̄²) / 3
    Δ = S^2 - 4 * ab
    Δ ≥ 0 || throw(DomainError((μ̄, σ̄², c),
        "TriangularDist: no real (a, b) for these moments and mode (discriminant=$Δ)"))
    a = (S - √Δ) / 2
    b = (S + √Δ) / 2
    (a ≤ c ≤ b) || throw(DomainError((μ̄, σ̄², c),
        "TriangularDist: solved (a=$a, b=$b) does not satisfy a ≤ c ≤ b"))
    return TriangularDist(a, b, c)
end

"""
    dist_from_mean_var_mode(::Type{DiscreteTriangular}, μ̄, σ̄², c)

Approximate. Solves the *continuous* triangular `(a, b)` for the requested
moments and mode `c`, rounds to integers, then searches a ±1 neighbourhood of
`(a, b)` to pick the integer combination whose `(mean, var)` minimise the
squared relative error against `(μ̄, σ̄²)`. The result will have mean and
variance close to `μ̄, σ̄²` but generally not matching exactly (3 integer
parameters vs. 3 continuous constraints).
"""
function dist_from_mean_var_mode(::Type{DiscreteTriangular}, μ̄::Number, σ̄²::Number, c::Number)
    cont = dist_from_mean_var_mode(TriangularDist, μ̄, σ̄², c)
    a0 = round(Int, minimum(cont))
    b0 = round(Int, maximum(cont))
    c_int = round(Int, c)

    best = nothing
    best_err = Inf
    for da in -1:1, db in -1:1
        a_try = min(a0 + da, c_int)
        b_try = max(b0 + db, c_int)
        a_try ≤ c_int ≤ b_try || continue
        d = DiscreteTriangular(a_try, b_try, c_int)
        # Squared relative error in mean and variance (both nonzero by feasibility).
        err = ((mean(d) - μ̄) / max(abs(μ̄), 1.0))^2 +
              ((var(d) - σ̄²) / max(σ̄², 1.0))^2
        if err < best_err
            best_err = err
            best = d
        end
    end
    return best
end

function dist_from_mode_var(::Type{Gamma}, m::Number, σ̄²::Number)
    m >= 0 || throw(DomainError(m, "Gamma: mode must be ≥ 0"))
    # mode = (α-1)θ, var = αθ²
    # From var: θ² = var/α → θ = √(var/α)
    # From mode: (α-1)√(var/α) = m → (α-1)²(var/α) = m²
    # → var(α² - 2α + 1)/α = m² → var·α - 2var + var/α = m²
    # Solve numerically
    sol = find_zero(
        α -> (α - 1) * √(σ̄² / α) - m,
        max(m^2 / σ̄² + 1, 1.1)  # initial guess
    )
    θ = √(σ̄² / sol)
    return Gamma(sol, θ)
end
