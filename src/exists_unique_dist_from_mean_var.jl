using Distributions
using SpecialFunctions

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
    elseif var≥μ(1-μ)
        throw(DomainError("Beta: the condition var < μ(1-μ) is not satisfied"))
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

function exists_unique_dist_from_mean_var(::Type{Chisq}, μ::Number, var::Number) end

function exists_unique_dist_from_mean_var(::Type{Erlang}, μ::Number, var::Number) end

function exists_unique_dist_from_mean_var(::Type{Exponential}, μ::Number, var::Number) end

function exists_unique_dist_from_mean_var(::Type{FDist}, μ::Number, var::Number) end

function exists_unique_dist_from_mean_var(::Type{Frechet}, μ::Number, var::Number) end

function exists_unique_dist_from_mean_var(::Type{Gamma}, μ::Number, var::Number) end

function exists_unique_dist_from_mean_var(::Type{Gumbel}, μ::Number, var::Number) end

function exists_unique_dist_from_mean_var(::Type{InverseGamma}, μ::Number, var::Number) end

function exists_unique_dist_from_mean_var(::Type{Laplace}, μ::Number, var::Number) end

function exists_unique_dist_from_mean_var(::Type{Logistic}, μ::Number, var::Number) end

function exists_unique_dist_from_mean_var(::Type{LogNormal}, μ::Number, var::Number) end

function exists_unique_dist_from_mean_var(::Type{Normal}, μ::Number, var::Number) end

function exists_unique_dist_from_mean_var(::Type{Rayleigh}, μ::Number, var::Number) end

function exists_unique_dist_from_mean_var(::Type{TDist}, μ::Number, var::Number) end

function exists_unique_dist_from_mean_var(::Type{TriangularDist}, μ::Number, var::Number) end

function exists_unique_dist_from_mean_var(::Type{Uniform}, μ::Number, var::Number) end

function exists_unique_dist_from_mean_var(::Type{Weibull}, μ::Number, var::Number) end