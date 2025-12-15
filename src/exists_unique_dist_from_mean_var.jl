function exists_unique_dist_from_mean_var(disttype::Type{<:Distribution}, mean::Number, var::Number)
    error("$disttype distribution not supported")
end

function exists_unique_dist_from_mean_std(disttype::Type{<:Distribution}, μ::Number, σ::Number)
    return exists_unique_dist_from_mean_var(disttype, μ, σ^2)
end

function exists_unique_dist_from_mean_var(::Type{Beta}, μ::Number, var::Number)
    if μ≤0
        throw(DomainError("Beta: the condition μ > zero(μ) is not satisfied"))
    elseif μ≥1
        throw(DomainError("Beta: the condition μ < one(μ) is not satisfied"))
    end
    if var≤0
        throw(DomainError("Beta: the condition var > zero(var) is not satisfied"))
    elseif var≥μ*(1-μ)
        throw(DomainError("Beta: the condition var < μ*(1-μ) is not satisfied"))
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

function exists_unique_dist_from_mean_var(::Type{Chi}, μ::Number, var::Number) 
    if var≤-μ^2
        throw(DomainError("Chi: the condition var > -μ^2 is not satisfied"))
    end
    if isapprox(μ, √(2)*gamma((μ^2+var+1)/2)/gamma((μ^2+var)/2); rtol=1e-10, atol=1e-12)
        throw(DomainError("Chi: the condition μ=√(2)*gamma((μ^2+var+1)/2)/gamma((μ^2+var)/2) is not satisfied (if μ^2+var>10 then letting var=1/2(1-μ^2) will generate a valid distribution)"))
    end
    return true
end

function exists_unique_dist_from_mean_var(::Type{Chisq}, μ::Number, var::Number)
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

function exists_unique_dist_from_mean_var(::Type{Erlang}, μ::Number, var::Number)
    if μ≤0
        throw(DomainError("Erlang: the condition μ > zero(μ) is not satisfied"))
    end
    if var≤0
        throw(DomainError("Erlang: the condition var > zero(var) is not satisfied"))
    elseif !isinteger(μ^2/var)
        throw(DomainError("Erlang: the condition μ^2/var ∈ ℕ is not satisfied"))
    end
    return true
end

function exists_unique_dist_from_mean_var(::Type{Exponential}, μ::Number, var::Number)
    if μ≤0
        throw(DomainError("Exponential: the condition μ > zero(μ) is not satisfied"))
    end
    if var≤0
        throw(DomainError("Exponential: the condition var > zero(var) is not satisfied"))
    elseif !isapprox(var, μ^2; rtol=1e-10)
        throw(DomainError("Exponential: the condition var = μ² is not satisfied"))
    end
    return true
end

function exists_unique_dist_from_mean_var(::Type{FDist}, μ::Number, var::Number) 
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

function exists_unique_dist_from_mean_var(::Type{Frechet}, μ::Number, var::Number)
    throw(DomainError("Frechet: no unique Frechet distribution for any mean, var pair (3 parameter definition)"))
end

function exists_unique_dist_from_mean_var(::Type{Gamma}, μ::Number, var::Number)
    if μ≤0
        throw(DomainError("Gamma: the condition μ > zero(μ) is not satisfied"))
    end
    if var≤0
        throw(DomainError("Gamma: the condition var > zero(var) is not satisfied"))
    end
    return true
end

function exists_unique_dist_from_mean_var(::Type{Gumbel}, μ::Number, var::Number)
    if var≤0
        throw(DomainError("Gumbel: the condition var > zero(var) is not satisfied"))
    end
    return true
end

function exists_unique_dist_from_mean_var(::Type{InverseGamma}, μ::Number, var::Number)
    if var≤0
        throw(DomainError("InverseGamma: the condition var > zero(var) is not satisfied"))
    end
    return true
end

function exists_unique_dist_from_mean_var(::Type{Laplace}, μ::Number, var::Number)
    if var≤0
        throw(DomainError("Laplace: the condition var > zero(var) is not satisfied"))
    end
    return true
end

function exists_unique_dist_from_mean_var(::Type{Logistic}, μ::Number, var::Number)
    if var≤0
        throw(DomainError("Logistic: the condition var > zero(var) is not satisfied"))
    end
    return true
end

function exists_unique_dist_from_mean_var(::Type{LogNormal}, μ::Number, var::Number)
    if var≤0
        throw(DomainError("LogNormal: the condition var > zero(var) is not satisfied"))
    end
    return true
end

function exists_unique_dist_from_mean_var(::Type{Normal}, μ::Number, var::Number)
    if var≤0
        throw(DomainError("Normal: the condition var > zero(var) is not satisfied"))
    end
    return true
end

function exists_unique_dist_from_mean_var(::Type{Rayleigh}, μ::Number, var::Number)
    if var≤0
        throw(DomainError("Rayleigh: the condition var > zero(var) is not satisfied"))
    elseif var≠(2-π)*μ^2/π
        throw(DomainError("Rayleigh: the condition var = (2-π)*μ^2/π is not satisfied"))
    end
    return true
end

function exists_unique_dist_from_mean_var(::Type{TDist}, μ::Number, var::Number)
    if μ≠0
        throw(DomainError("TDist: the condition μ = zero(μ) is not satisfied"))
    end
    if var≤1
        throw(DomainError("Rayleigh: the condition var > one(var) is not satisfied"))
    end
    return true
end

function exists_unique_dist_from_mean_var(::Type{TriangularDist}, μ::Number, var::Number)
    throw(DomainError("TriangularDist: no unique TriangularDist distribution for any mean, var pair (3 parameter definition)"))
end

function exists_unique_dist_from_mean_var(::Type{Uniform}, μ::Number, var::Number)
    if var≤0
        throw(DomainError("Uniform: the condition var > zero(var) is not satisfied"))
    end
    return true
end

function exists_unique_dist_from_mean_var(::Type{Weibull}, μ::Number, var::Number)
    if var≤0
        throw(DomainError("Weibull: the condition var > zero(var) is not satisfied"))
    end
    return true
end