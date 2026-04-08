function exists_unique_dist_from_mean_var(disttype::Type{<:Distribution}, mean::Number, var::Number)
    error("$disttype distribution not supported")
end

function base_exists_unique_dist_from_mean_var(disttype::Type{<:Distribution}, mean::Number, var::Number)
    if var≤0
        error("$disttype: the condition var > zero(var) is not satisfied")
    end
end


function exists_unique_dist_from_mean_var(::Type{Beta}, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(Beta, μ,var)

    if μ≤0
        throw(DomainError("Beta: the condition μ > zero(μ) is not satisfied"))
    elseif μ≥1
        throw(DomainError("Beta: the condition μ < one(μ) is not satisfied"))
    end
    
    if var≥μ*(1-μ)
        throw(DomainError("Beta: the condition var < μ*(1-μ) is not satisfied"))
    end
    return true
end

function exists_unique_dist_from_mean_var(::Type{Uniform}, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(Uniform, μ,var)
    return true
end


function exists_unique_dist_from_mean_var(::Type{Normal}, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(Normal, μ,var)
    return true
end


function exists_unique_dist_from_mean_var(::Type{TDist}, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(TDist, μ,var)
    if μ≠0
        throw(DomainError("TDist: the condition μ = zero(μ) is not satisfied"))
    end
    if var≤1
        throw(DomainError("TDist: the condition var > one(var) is not satisfied"))
    end
    return true
end

function exists_unique_dist_from_mean_var(d::TDist, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(TDist, μ,var)
    ν = dof(d)
    if ν ≤ 2
        throw(DomainError("TDist instance: the condition ν > 2 is not satisfied (ν=$ν). Variance is undefined for ν ≤ 2."))
    end
    return true
end

function exists_unique_dist_from_mean_var(::Type{Cauchy}, μ::Number, var::Number)
    if !isnan(μ)
        throw(DomainError("Cauchy: distribution cannot have a defined μ"))
    end
    if !isnan(var)
        throw(DomainError("Cauchy: distribution cannot have a defined var"))
    end
    throw(DomainError("Cauchy: no unique distribution exists"))
end


function exists_unique_dist_from_mean_var(::Type{Logistic}, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(Logistic, μ,var)
    return true
end

function exists_unique_dist_from_mean_var(::Type{Laplace}, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(Laplace, μ,var)
    return true
end


function exists_unique_dist_from_mean_var(::Type{LogNormal}, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(LogNormal, μ,var)
    if μ≤0
        throw(DomainError("LogNormal: the condition μ > zero(μ) is not satisfied"))
    end
    return true
end


function exists_unique_dist_from_mean_var(::Type{Chisq}, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(Chisq, μ,var)
    if μ≤0
        throw(DomainError("Chisq: the condition μ > zero(μ) is not satisfied"))
    elseif !isinteger(μ)
        throw(DomainError("Chisq: the condition μ ∈ ℕ is not satisfied"))
    end
    if var ≠ 2μ
        throw(DomainError("Chisq: the condition var=2μ is not satisfied"))
    end
    return true
end

function exists_unique_dist_from_mean_var(::Type{Exponential}, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(Exponential, μ,var)
    if μ≤0
        throw(DomainError("Exponential: the condition μ > zero(μ) is not satisfied"))
    end
    if !isapprox(var, μ^2; rtol=1e-10)
        throw(DomainError("Exponential: the condition var = μ² is not satisfied"))
    end
    return true
end

function exists_unique_dist_from_mean_var(::Type{Gamma}, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(Gamma, μ,var)
    if μ≤0
        throw(DomainError("Gamma: the condition μ > zero(μ) is not satisfied"))
    end
    return true
end

function exists_unique_dist_from_mean_var(::Type{Erlang}, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(Erlang, μ,var)
    if μ≤0
        throw(DomainError("Erlang: the condition μ > zero(μ) is not satisfied"))
    end
    return true
end


function exists_unique_dist_from_mean_var(::Type{Frechet}, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(Frechet, μ,var)
    if μ≤0
        throw(DomainError("Frechet: the condition μ > zero(μ) is not satisfied"))
    end
    return true
end

function exists_unique_dist_from_mean_var(::Type{Weibull}, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(Weibull, μ,var)
    if μ≤0
        throw(DomainError("Weibull: the condition μ > zero(μ) is not satisfied"))
    end
    return true
end

function exists_unique_dist_from_mean_var(::Type{Gumbel}, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(Gumbel, μ,var)
    return true
end


function exists_unique_dist_from_mean_var(::Type{Chi}, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(Chi, μ,var) 
    if isapprox(μ, √(2)*gamma((μ^2+var+1)/2)/gamma((μ^2+var)/2); rtol=1e-10, atol=1e-12)
        throw(DomainError("Chi: the condition μ=√(2)*gamma((μ^2+var+1)/2)/gamma((μ^2+var)/2) is not satisfied (if μ^2+var>10 then letting var=1/2(1-μ^2) will generate a valid distribution)"))
    end
    return true
end

function exists_unique_dist_from_mean_var(::Type{Rayleigh}, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(Rayleigh, μ,var)
    if √(var)/μ≠√((4-π)/π)
        throw(DomainError("Rayleigh: the condition CV = √((4-π)/π) is not satisfied"))
    end
    return true
end


function exists_unique_dist_from_mean_var(::Type{FDist}, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(FDist, μ,var)
    if μ≤1
        throw(DomainError("FDist: the condition μ > one(μ) is not satisfied"))
    elseif μ≥2
        throw(DomainError("FDist: the condition μ < 2 is not satisfied"))
    end
    if var≤μ^2*(μ-1)/(2-μ)
        throw(DomainError("FDist: the condition var > μ^2*(μ-1)/(2-μ) is not satisfied"))
    end
    return true
end


function exists_unique_dist_from_mean_var(::Type{InverseGamma}, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(InverseGamma, μ,var)
    return true
end


function exists_unique_dist_from_mean_var(::Type{Binomial}, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(Binomial, μ,var)
    if var≥μ
        throw(DomainError("Binomial: the condition μ > var is not satisfied"))
    end
    if μ≤0
        throw(DomainError("Binomial: the condition μ > zero(μ) is not satisfied"))
    elseif !isinteger(μ^2/(μ-var))
        throw(DomainError("Binomial: the condition μ^2/(μ-var) ∈ ℕ is not satisfied"))
    end
    return true 
end


function exists_unique_dist_from_mean_var(::Type{Poisson}, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(Poisson, μ,var)
    if μ≠var
        throw(DomainError("Poisson: the condition μ=var is not satisfied"))
    end
    return true 
end


function exists_unique_dist_from_mean_var(::Type{NegativeBinomial}, μ::Number, var::Number)
    base_exists_unique_dist_from_mean_var(NegativeBinomial, μ,var)
    if μ≥var
        throw(DomainError("NegativeBinomial: the condition μ<var is not satisfied"))
    end
    return true
end


# --- Not yet implemented ---

function exists_unique_dist_from_mean_var(::Type{Pareto}, μ::Number, var::Number)
    throw(ErrorException("Pareto: exists_unique_dist_from_mean_var not yet implemented"))
end

function exists_unique_dist_from_mean_var(::Type{FoldedNormal}, μ::Number, var::Number)
    throw(ErrorException("FoldedNormal: exists_unique_dist_from_mean_var not yet implemented"))
end

function exists_unique_dist_from_mean_var(::Type{Geometric}, μ::Number, var::Number)
    throw(ErrorException("Geometric: exists_unique_dist_from_mean_var not yet implemented"))
end

function exists_unique_dist_from_mean_var(d::Truncated{<:Normal}, μ::Number, var::Number)
    throw(ErrorException("Truncated Normal: exists_unique_dist_from_mean_var not yet implemented"))
end

function exists_unique_dist_from_mean_var(d::Truncated{<:Laplace}, μ::Number, var::Number)
    throw(ErrorException("Truncated Laplace: exists_unique_dist_from_mean_var not yet implemented"))
end

function exists_unique_dist_from_mean_var(d::Truncated{<:Logistic}, μ::Number, var::Number)
    throw(ErrorException("Truncated Logistic: exists_unique_dist_from_mean_var not yet implemented"))
end

function exists_unique_dist_from_mean_var(::Type{TriangularDist}, μ::Number, var::Number)
    throw(ErrorException("TriangularDist: exists_unique_dist_from_mean_var not yet implemented"))
end

function exists_unique_dist_from_mean_var(::Type{SymTriangularDist}, μ::Number, var::Number)
    throw(ErrorException("SymTriangularDist: exists_unique_dist_from_mean_var not yet implemented"))
end

function exists_unique_dist_from_mean_var(::Type{DiscreteTriangular}, μ::Number, var::Number)
    throw(ErrorException("DiscreteTriangular: exists_unique_dist_from_mean_var not yet implemented"))
end

function exists_unique_dist_from_mean_var(::Type{DiscreteSymmetricTriangular}, μ::Number, var::Number)
    throw(ErrorException("DiscreteSymmetricTriangular: exists_unique_dist_from_mean_var not yet implemented"))
end

function exists_unique_dist_from_mean_var(::Type{TruncatedPoisson}, μ::Number, var::Number)
    throw(ErrorException("TruncatedPoisson: exists_unique_dist_from_mean_var not yet implemented"))
end

function exists_unique_dist_from_mean_var(::Type{DiscreteUniform}, μ::Number, var::Number)
    throw(ErrorException("DiscreteUniform: exists_unique_dist_from_mean_var not yet implemented"))
end