"""
    exists_dist_from_mean_var(D, μ̄, σ̄²)

Check whether a unique distribution of type `D` exists for mean `μ̄` and variance `σ̄²`.
Returns `true` if feasible, throws `DomainError` otherwise.
"""
function exists_dist_from_mean_var(disttype::Type{<:Distribution}, μ̄::Number, σ̄²::Number)
    error("$disttype distribution not supported")
end

function base_exists_dist_from_mean_var(disttype::Type{<:Distribution}, μ̄::Number, σ̄²::Number)
    if σ̄²≤0
        error("$disttype: the condition σ̄² > 0 is not satisfied")
    end
end


function exists_dist_from_mean_var(::Type{Beta}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(Beta, μ̄, σ̄²)

    if μ̄≤0
        throw(DomainError("Beta: the condition μ̄ > 0 is not satisfied"))
    elseif μ̄≥1
        throw(DomainError("Beta: the condition μ̄ < 1 is not satisfied"))
    end

    if σ̄²≥μ̄*(1-μ̄)
        throw(DomainError("Beta: the condition σ̄² < μ̄(1-μ̄) is not satisfied"))
    end
    return true
end

function exists_dist_from_mean_var(::Type{Uniform}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(Uniform, μ̄, σ̄²)
    return true
end


function exists_dist_from_mean_var(::Type{Normal}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(Normal, μ̄, σ̄²)
    return true
end


function exists_dist_from_mean_var(::Type{TDist}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(TDist, μ̄, σ̄²)
    if μ̄≠0
        throw(DomainError("TDist: the condition μ̄ = 0 is not satisfied"))
    end
    if σ̄²≤1
        throw(DomainError("TDist: the condition σ̄² > 1 is not satisfied"))
    end
    return true
end

function exists_dist_from_mean_var(d::TDist, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(TDist, μ̄, σ̄²)
    ν = dof(d)
    if ν ≤ 2
        throw(DomainError("TDist instance: the condition ν > 2 is not satisfied (ν=$ν). Variance is undefined for ν ≤ 2."))
    end
    return true
end

function exists_dist_from_mean_var(::Type{Cauchy}, μ̄::Number, σ̄²::Number)
    if !isnan(μ̄)
        throw(DomainError("Cauchy: distribution cannot have a defined μ̄"))
    end
    if !isnan(σ̄²)
        throw(DomainError("Cauchy: distribution cannot have a defined σ̄²"))
    end
    throw(DomainError("Cauchy: no unique distribution exists"))
end


function exists_dist_from_mean_var(::Type{Logistic}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(Logistic, μ̄, σ̄²)
    return true
end

function exists_dist_from_mean_var(::Type{Laplace}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(Laplace, μ̄, σ̄²)
    return true
end


function exists_dist_from_mean_var(::Type{LogNormal}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(LogNormal, μ̄, σ̄²)
    if μ̄≤0
        throw(DomainError("LogNormal: the condition μ̄ > 0 is not satisfied"))
    end
    return true
end


function exists_dist_from_mean_var(::Type{Chisq}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(Chisq, μ̄, σ̄²)
    if μ̄≤0
        throw(DomainError("Chisq: the condition μ̄ > 0 is not satisfied"))
    elseif !isinteger(μ̄)
        throw(DomainError("Chisq: the condition μ̄ ∈ ℕ is not satisfied"))
    end
    if σ̄² ≠ 2μ̄
        throw(DomainError("Chisq: the condition σ̄² = 2μ̄ is not satisfied"))
    end
    return true
end

function exists_dist_from_mean_var(::Type{Exponential}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(Exponential, μ̄, σ̄²)
    if μ̄≤0
        throw(DomainError("Exponential: the condition μ̄ > 0 is not satisfied"))
    end
    if !isapprox(σ̄², μ̄^2; rtol=1e-10)
        throw(DomainError("Exponential: the condition σ̄² = μ̄² is not satisfied"))
    end
    return true
end

function exists_dist_from_mean_var(::Type{Gamma}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(Gamma, μ̄, σ̄²)
    if μ̄≤0
        throw(DomainError("Gamma: the condition μ̄ > 0 is not satisfied"))
    end
    return true
end

function exists_dist_from_mean_var(::Type{Erlang}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(Erlang, μ̄, σ̄²)
    if μ̄≤0
        throw(DomainError("Erlang: the condition μ̄ > 0 is not satisfied"))
    end
    return true
end


function exists_dist_from_mean_var(::Type{Frechet}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(Frechet, μ̄, σ̄²)
    if μ̄≤0
        throw(DomainError("Frechet: the condition μ̄ > 0 is not satisfied"))
    end
    return true
end

function exists_dist_from_mean_var(::Type{Weibull}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(Weibull, μ̄, σ̄²)
    if μ̄≤0
        throw(DomainError("Weibull: the condition μ̄ > 0 is not satisfied"))
    end
    return true
end

function exists_dist_from_mean_var(::Type{Gumbel}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(Gumbel, μ̄, σ̄²)
    return true
end


function exists_dist_from_mean_var(::Type{Chi}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(Chi, μ̄, σ̄²)
    if isapprox(μ̄, √(2)*gamma((μ̄^2+σ̄²+1)/2)/gamma((μ̄^2+σ̄²)/2); rtol=1e-10, atol=1e-12)
        throw(DomainError("Chi: the condition μ̄=√(2)Γ((μ̄²+σ̄²+1)/2)/Γ((μ̄²+σ̄²)/2) is not satisfied"))
    end
    return true
end

function exists_dist_from_mean_var(::Type{Rayleigh}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(Rayleigh, μ̄, σ̄²)
    if √(σ̄²)/μ̄≠√((4-π)/π)
        throw(DomainError("Rayleigh: the condition CV = √((4-π)/π) is not satisfied"))
    end
    return true
end


function exists_dist_from_mean_var(::Type{FDist}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(FDist, μ̄, σ̄²)
    if μ̄≤1
        throw(DomainError("FDist: the condition μ̄ > 1 is not satisfied"))
    elseif μ̄≥2
        throw(DomainError("FDist: the condition μ̄ < 2 is not satisfied"))
    end
    if σ̄²≤μ̄^2*(μ̄-1)/(2-μ̄)
        throw(DomainError("FDist: the condition σ̄² > μ̄²(μ̄-1)/(2-μ̄) is not satisfied"))
    end
    return true
end


function exists_dist_from_mean_var(::Type{InverseGamma}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(InverseGamma, μ̄, σ̄²)
    return true
end


function exists_dist_from_mean_var(::Type{Binomial}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(Binomial, μ̄, σ̄²)
    if σ̄²≥μ̄
        throw(DomainError("Binomial: the condition μ̄ > σ̄² is not satisfied"))
    end
    if μ̄≤0
        throw(DomainError("Binomial: the condition μ̄ > 0 is not satisfied"))
    elseif !isapprox(μ̄^2/(μ̄-σ̄²), round(μ̄^2/(μ̄-σ̄²)); atol=1e-8)
        throw(DomainError("Binomial: the condition μ̄²/(μ̄-σ̄²) ∈ ℕ is not satisfied"))
    end
    return true
end


function exists_dist_from_mean_var(::Type{Poisson}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(Poisson, μ̄, σ̄²)
    if μ̄≠σ̄²
        throw(DomainError("Poisson: the condition μ̄ = σ̄² is not satisfied"))
    end
    return true
end


function exists_dist_from_mean_var(::Type{NegativeBinomial}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(NegativeBinomial, μ̄, σ̄²)
    if μ̄≥σ̄²
        throw(DomainError("NegativeBinomial: the condition μ̄ < σ̄² is not satisfied"))
    end
    return true
end


# --- Recently implemented ---

function exists_dist_from_mean_var(::Type{Pareto}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(Pareto, μ̄, σ̄²)
    if μ̄ ≤ 0
        throw(DomainError("Pareto: the condition μ̄ > 0 is not satisfied"))
    end
    CV² = σ̄² / μ̄^2
    α = 1 + √(1 + 1/CV²)
    if α ≤ 2
        throw(DomainError("Pareto: the condition α > 2 is not satisfied (variance is infinite for α ≤ 2)"))
    end
    return true
end

function exists_dist_from_mean_var(::Type{FoldedNormal}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(FoldedNormal, μ̄, σ̄²)
    if μ̄ ≤ 0
        throw(DomainError("FoldedNormal: the condition μ̄ > 0 is not satisfied"))
    end
    return true
end

function exists_dist_from_mean_var(::Type{Geometric}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(Geometric, μ̄, σ̄²)
    if μ̄ ≤ 0
        throw(DomainError("Geometric: the condition μ̄ > 0 is not satisfied"))
    end
    if !isapprox(σ̄², μ̄ * (1 + μ̄); rtol=1e-10)
        throw(DomainError("Geometric: the condition σ̄² = μ̄(1+μ̄) is not satisfied"))
    end
    return true
end

function exists_dist_from_mean_var(d::Truncated{<:Normal}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(Normal, μ̄, σ̄²)
    lo, hi = extrema(d)
    if μ̄ ≤ lo || μ̄ ≥ hi
        throw(DomainError("Truncated Normal: μ̄ must be in ($lo, $hi)"))
    end
    return true
end

function exists_dist_from_mean_var(d::Truncated{<:Laplace}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(Laplace, μ̄, σ̄²)
    lo, hi = extrema(d)
    if μ̄ ≤ lo || μ̄ ≥ hi
        throw(DomainError("Truncated Laplace: μ̄ must be in ($lo, $hi)"))
    end
    return true
end

function exists_dist_from_mean_var(d::Truncated{<:Logistic}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(Logistic, μ̄, σ̄²)
    lo, hi = extrema(d)
    if μ̄ ≤ lo || μ̄ ≥ hi
        throw(DomainError("Truncated Logistic: μ̄ must be in ($lo, $hi)"))
    end
    return true
end

function exists_dist_from_mean_var(::Type{TriangularDist}, μ̄::Number, σ̄²::Number)
    throw(ErrorException("TriangularDist: exists_dist_from_mean_var not yet implemented"))
end

function exists_dist_from_mean_var(::Type{SymTriangularDist}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(SymTriangularDist, μ̄, σ̄²)
    return true
end

function exists_dist_from_mean_var(::Type{DiscreteTriangular}, μ̄::Number, σ̄²::Number)
    throw(ErrorException("DiscreteTriangular: exists_dist_from_mean_var not yet implemented"))
end

function exists_dist_from_mean_var(::Type{DiscreteSymmetricTriangular}, μ̄::Number, σ̄²::Number)
    throw(ErrorException("DiscreteSymmetricTriangular: exists_dist_from_mean_var not yet implemented"))
end

function exists_dist_from_mean_var(::Type{TruncatedPoisson}, μ̄::Number, σ̄²::Number)
    throw(ErrorException("TruncatedPoisson: exists_dist_from_mean_var not yet implemented"))
end

function exists_dist_from_mean_var(::Type{DiscreteUniform}, μ̄::Number, σ̄²::Number)
    base_exists_dist_from_mean_var(DiscreteUniform, μ̄, σ̄²)
    n_raw = -1 + √(1 + 12 * σ̄²)
    if !isapprox(n_raw, round(n_raw); atol=1e-8) || round(n_raw) < 0
        throw(DomainError("DiscreteUniform: n = b - a must be a non-negative integer (got n ≈ $n_raw)"))
    end
    n = round(Int, n_raw)
    a_raw = μ̄ - n / 2
    if !isapprox(a_raw, round(a_raw); atol=1e-8)
        throw(DomainError("DiscreteUniform: lower bound a must be an integer (got a ≈ $a_raw)"))
    end
    return true
end





"""
    dist_from_mean_var(D, μ̄, σ̄²)

Construct a distribution of type `D` with the given mean `μ̄` and variance `σ̄²`.

Dispatches on the distribution type (or instance for truncated/TDist).
Throws `DomainError` if no valid distribution exists for the given moments.
Use `exists_dist_from_mean_var` to check feasibility before calling.

See also: [`dist_from_mean_std`](@ref), [`dist_from_mean_cv`](@ref),
[`exists_dist_from_mean_var`](@ref)
"""
function dist_from_mean_var end

"""
    dist_from_mean_var(::Type{Beta}, μ̄, σ̄²)

Direct formula. Construct a `Beta(α, β)` distribution.
Requires `0 < μ̄ < 1` and `0 < σ̄² < μ̄(1-μ̄)`.
"""
function dist_from_mean_var(::Type{Beta}, μ̄::Number, σ̄²::Number)
    exists_dist_from_mean_var(Beta, μ̄, σ̄²)
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
    exists_dist_from_mean_var(Uniform, μ̄, σ̄²)
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
    exists_dist_from_mean_var(Normal, μ̄, σ̄²)
    return Normal(μ̄,√(σ̄²))
end

"""
    dist_from_mean_var(::Type{TDist}, μ̄, σ̄²)

Direct formula. Construct a standard `TDist(ν)` distribution.
Requires `μ̄ = 0` and `σ̄² > 1`.
"""
function dist_from_mean_var(::Type{TDist}, μ̄::Number, σ̄²::Number)
    exists_dist_from_mean_var(TDist, μ̄, σ̄²)
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
    exists_dist_from_mean_var(d, μ̄, σ̄²)
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
    exists_dist_from_mean_var(Cauchy, μ̄, σ̄²)
end

"""
    dist_from_mean_var(::Type{Logistic}, μ̄, σ̄²)

Direct formula. Construct a `Logistic(μ, s)` distribution. Any `μ̄ ∈ ℝ` and `σ̄² > 0`.
"""
function dist_from_mean_var(::Type{Logistic}, μ̄::Number, σ̄²::Number)
    exists_dist_from_mean_var(Logistic, μ̄, σ̄²)
    x̄=μ̄
    s=√(3*σ̄²/π^2)
    return Logistic(x̄,s)
end

"""
    dist_from_mean_var(::Type{Laplace}, μ̄, σ̄²)

Direct formula. Construct a `Laplace(μ, b)` distribution. Any `μ̄ ∈ ℝ` and `σ̄² > 0`.
"""
function dist_from_mean_var(::Type{Laplace}, μ̄::Number, σ̄²::Number)
    exists_dist_from_mean_var(Laplace, μ̄, σ̄²)
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
    exists_dist_from_mean_var(LogNormal, μ̄, σ̄²)
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
    exists_dist_from_mean_var(Chisq, μ̄, σ̄²)
    return Chisq(μ̄)
end

"""
    dist_from_mean_var(::Type{Exponential}, μ̄, σ̄²)

Direct formula. Construct an `Exponential(μ)` distribution.
Requires `μ̄ > 0` and `σ̄² = μ̄²` (1 DOF).
"""
function dist_from_mean_var(::Type{Exponential}, μ̄::Number, σ̄²::Number)
    exists_dist_from_mean_var(Exponential, μ̄, σ̄²)
    return Exponential(μ̄)
end

"""
    dist_from_mean_var(::Type{Gamma}, μ̄, σ̄²)

Direct formula. Construct a `Gamma(α, θ)` distribution.
Requires `μ̄ > 0` and `σ̄² > 0`.
"""
function dist_from_mean_var(::Type{Gamma}, μ̄::Number, σ̄²::Number)
    exists_dist_from_mean_var(Gamma, μ̄, σ̄²)
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
    exists_dist_from_mean_var(Erlang, μ̄, σ̄²)
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
    exists_dist_from_mean_var(Frechet, μ̄, σ̄²)
    return _frechet_from_mean_var(μ̄,σ̄²)
end

"""
    dist_from_mean_var(::Type{Weibull}, μ̄, σ̄²)

Numerical (root-finding). Construct a `Weibull(k, λ)` distribution.
Requires `μ̄ > 0` and `σ̄² > 0`. Solves the beta-ratio equation.
"""
function dist_from_mean_var(::Type{Weibull}, μ̄::Number, σ̄²::Number)
    exists_dist_from_mean_var(Weibull, μ̄, σ̄²)
    return _weibull_from_mean_var(μ̄,σ̄²)
end

"""
    dist_from_mean_var(::Type{Gumbel}, μ̄, σ̄²)

Direct formula. Construct a `Gumbel(μ_loc, β)` distribution.
Any `μ̄ ∈ ℝ` and `σ̄² > 0`.
"""
function dist_from_mean_var(::Type{Gumbel}, μ̄::Number, σ̄²::Number)
    exists_dist_from_mean_var(Gumbel, μ̄, σ̄²)
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
    exists_dist_from_mean_var(Chi, μ̄, σ̄²)
    return Chi(μ̄^2+σ̄²)
end

"""
    dist_from_mean_var(::Type{Rayleigh}, μ̄, σ̄²)

Direct formula. Construct a `Rayleigh(σ)` distribution.
Requires `μ̄ > 0` and `CV = √((4-π)/π)` (1 DOF).
"""
function dist_from_mean_var(::Type{Rayleigh}, μ̄::Number, σ̄²::Number)
    exists_dist_from_mean_var(Rayleigh, μ̄, σ̄²)
    σ=√(2/π)*μ̄
    return Rayleigh(σ)
end

"""
    dist_from_mean_var(::Type{FDist}, μ̄, σ̄²)

Direct formula. Construct an `FDist(ν₁, ν₂)` distribution.
Requires `1 < μ̄ < 2` and `σ̄² > μ̄²(μ̄-1)/(2-μ̄)`.
"""
function dist_from_mean_var(::Type{FDist}, μ̄::Number, σ̄²::Number)
    exists_dist_from_mean_var(FDist, μ̄, σ̄²)
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
    exists_dist_from_mean_var(InverseGamma, μ̄, σ̄²)
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
    exists_dist_from_mean_var(Binomial, μ̄, σ̄²)
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
    exists_dist_from_mean_var(Poisson, μ̄, σ̄²)
    return Poisson(μ̄)
end

"""
    dist_from_mean_var(::Type{NegativeBinomial}, μ̄, σ̄²)

Direct formula. Construct a `NegativeBinomial(r, p)` distribution.
Requires `σ̄² > μ̄ > 0`.
"""
function dist_from_mean_var(::Type{NegativeBinomial}, μ̄::Number, σ̄²::Number)
    exists_dist_from_mean_var(NegativeBinomial, μ̄, σ̄²)
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
    exists_dist_from_mean_var(Pareto, μ̄, σ̄²)
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
    exists_dist_from_mean_var(FoldedNormal, μ̄, σ̄²)
    μp, σp = _solve_folded_normal(Float64(μ̄), Float64(σ̄²))
    return Normal(μp, σp)  # returns the parent Normal; user takes |X|
end

"""
    dist_from_mean_var(::Type{Geometric}, μ̄, σ̄²)

Direct formula. Construct a `Geometric(p)` distribution.
Requires `μ̄ > 0` and `σ̄² = μ̄(1+μ̄)` (1 DOF).
"""
function dist_from_mean_var(::Type{Geometric}, μ̄::Number, σ̄²::Number)
    exists_dist_from_mean_var(Geometric, μ̄, σ̄²)
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
    exists_dist_from_mean_var(SymTriangularDist, μ̄, σ̄²)
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
    exists_dist_from_mean_var(DiscreteUniform, μ̄, σ̄²)
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



# Variants that convert to mean+var and delegate

"""
    dist_from_mean_std(D, μ̄, σ̄)

Construct distribution `D` from mean `μ̄` and standard deviation `σ̄`.
Delegates to `dist_from_mean_var(D, μ̄, σ̄²)`.
"""
dist_from_mean_std(D, μ̄::Number, σ̄::Number) = dist_from_mean_var(D, μ̄, σ̄^2)

"""
    dist_from_mean_cv(D, μ̄, cv)

Construct distribution `D` from mean `μ̄` and coefficient of variation `cv = σ̄/μ̄`.
Delegates to `dist_from_mean_var(D, μ̄, (cv·μ̄)²)`.
"""
dist_from_mean_cv(D, μ̄::Number, cv::Number) = dist_from_mean_var(D, μ̄, (cv * μ̄)^2)

"""
    dist_from_mean_scv(D, μ̄, scv)

Construct distribution `D` from mean `μ̄` and squared coefficient of variation `scv = σ̄²/μ̄²`.
Delegates to `dist_from_mean_var(D, μ̄, scv·μ̄²)`.
"""
dist_from_mean_scv(D, μ̄::Number, scv::Number) = dist_from_mean_var(D, μ̄, scv * μ̄^2)

"""
    dist_from_mean_second_moment(D, μ̄, m2)

Construct distribution `D` from mean `μ̄` and second moment `m2 = E[X²]`.
Delegates to `dist_from_mean_var(D, μ̄, m2 - μ̄²)`.
"""
dist_from_mean_second_moment(D, μ̄::Number, m2::Number) = dist_from_mean_var(D, μ̄, m2 - μ̄^2)

"""
    exists_dist_from_mean_std(D, μ̄, σ̄)

Check feasibility for `dist_from_mean_std`. See [`exists_dist_from_mean_var`](@ref).
"""
exists_dist_from_mean_std(D, μ̄::Number, σ̄::Number) = exists_dist_from_mean_var(D, μ̄, σ̄^2)

"""
    exists_dist_from_mean_cv(D, μ̄, cv)

Check feasibility for `dist_from_mean_cv`. See [`exists_dist_from_mean_var`](@ref).
"""
exists_dist_from_mean_cv(D, μ̄::Number, cv::Number) = exists_dist_from_mean_var(D, μ̄, (cv * μ̄)^2)

"""
    exists_dist_from_mean_scv(D, μ̄, scv)

Check feasibility for `dist_from_mean_scv`. See [`exists_dist_from_mean_var`](@ref).
"""
exists_dist_from_mean_scv(D, μ̄::Number, scv::Number) = exists_dist_from_mean_var(D, μ̄, scv * μ̄^2)

"""
    exists_dist_from_mean_second_moment(D, μ̄, m2)

Check feasibility for `dist_from_mean_second_moment`. See [`exists_dist_from_mean_var`](@ref).
"""
exists_dist_from_mean_second_moment(D, μ̄::Number, m2::Number) = exists_dist_from_mean_var(D, μ̄, m2 - μ̄^2)

# Variants that only take a dispersion measure (for 1-parameter distributions where mean is determined)

"""
    dist_from_var(D, σ̄²)

Construct a 1-parameter distribution `D` from variance alone. Only supported for
distributions where the mean is uniquely determined by the variance:
`Exponential`, `Poisson`, `Chisq`, `Rayleigh`, `Geometric`.
"""
dist_from_var(D, σ̄²::Number) = dist_from_mean_var(D, _mean_from_var(D, σ̄²), σ̄²)

"""
    dist_from_std(D, σ̄)

Construct a 1-parameter distribution `D` from standard deviation alone.
Delegates to `dist_from_var(D, σ̄²)`.
"""
dist_from_std(D, σ̄::Number) = dist_from_var(D, σ̄^2)

# mean-from-var helpers for 1-parameter distributions where var determines mean
_mean_from_var(::Type{Exponential}, σ̄²::Number) = √σ̄²                         # σ̄² = μ̄²
_mean_from_var(::Type{Poisson}, σ̄²::Number) = σ̄²                               # σ̄² = μ̄
_mean_from_var(::Type{Chisq}, σ̄²::Number) = σ̄² / 2                             # σ̄² = 2μ̄
_mean_from_var(::Type{Rayleigh}, σ̄²::Number) = √(σ̄² * π / (4 - π))            # σ̄² = μ̄²(4-π)/π
_mean_from_var(::Type{Geometric}, σ̄²::Number) = (-1 + √(1 + 4σ̄²)) / 2        # σ̄² = μ̄(1+μ̄)
_mean_from_var(::Type{D}, σ̄²::Number) where {D<:Distribution} =
    throw(ErrorException("$D: dist_from_var not supported (mean is not determined by variance alone)"))

# Quantile convenience wrappers (delegate to dist_from_quantile / dist_from_quantiles)

"""
    dist_from_median(D, m)

Construct distribution `D` from its median. Equivalent to `dist_from_quantile(D, 0.5, m)`.
"""
dist_from_median(D, m::Number) = dist_from_quantile(D, 0.5, m)

"""
    dist_from_q1(D, q)

Construct distribution `D` from its first quartile. Equivalent to `dist_from_quantile(D, 0.25, q)`.
"""
dist_from_q1(D, q::Number) = dist_from_quantile(D, 0.25, q)

"""
    dist_from_q3(D, q)

Construct distribution `D` from its third quartile. Equivalent to `dist_from_quantile(D, 0.75, q)`.
"""
dist_from_q3(D, q::Number) = dist_from_quantile(D, 0.75, q)

"""
    dist_from_q1_q3(D, q1, q3)

Construct distribution `D` from its first and third quartiles.
"""
dist_from_q1_q3(D, q1::Number, q3::Number) = dist_from_quantiles(D, 0.25, q1, 0.75, q3)

"""
    dist_from_median_iqr(D, median, iqr)

Construct distribution `D` from its median and interquartile range.
"""
dist_from_median_iqr(D, median::Number, iqr::Number) = dist_from_quantiles(D, 0.25, median - iqr / 2, 0.75, median + iqr / 2)

"""
    dist_from_mean_median(D, μ̄, median)

Construct distribution `D` from its mean `μ̄` and median.
"""
dist_from_mean_median(D, μ̄::Number, median::Number) = dist_from_mean_quantile(D, μ̄, 0.5, median)
