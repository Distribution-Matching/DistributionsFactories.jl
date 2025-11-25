using Distributions

function dist_from_mean_var(disttype::Type{<:Distribution}, mean::Number, var::Number) end

function dist_from_mean_std(disttype::Type{<:Distribution}, μ::Number, σ::Number)
    return dist_from_mean_var(disttype, μ, σ^2)
end

function dist_from_mean_var(::Type{Arcsine}, μ::Number, var::Number)
    b = √(2*var)+μ
    a = 2*μ-b
    return Arcsine(a,b)
end

function dist_from_mean_var(::Type{Gamma}, μ::Number, var::Number)
    θ = var/μ
    α = 1/var
    return Gamma(α,θ)
end

function dist_from_mean_var(::Type{Beta}, μ::Number, var::Number)
    S = (μ*(1-μ))/var-1
    α = μ*S
    β = (1-μ)*S
    return Beta(α,β)
end

function dist_from_mean_var(::Type{BetaPrime}, μ::Number, var::Number)
    S = (μ*(μ+1))/var+1
    α = μ*S
    β = S+1
    return BetaPrime(α,β)
end

function dist_from_mean_var(::Type{Biweight}, μ::Number, var::Number) end

function dist_from_mean_var(::Type{Normal}, μ::Number, var::Number) 
    return Normal(mean, var)
end

function dist_from_mean_var(::Type{Uniform}, μ::Number, var::Number)
    b = √(3*var)+μ
    a = 2*μ-b
    return Uniform(a,b)
end
