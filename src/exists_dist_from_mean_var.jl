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
