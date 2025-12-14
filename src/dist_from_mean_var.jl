using Distributions
function dist_from_mean_var(::Type{Beta}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Beta, μ, var)
    S = (μ*(1-μ))/var-1
    α = μ*S
    β = (1-μ)*S
    return Beta(α,β)
end

function dist_from_mean_var(::Type{Cauchy}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Cauchy, μ, var)
end

function dist_from_mean_var(::Type{Chi}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Chi, μ, var)
    return Chis(μ^2+var)
end

function dist_from_mean_var(::Type{Chisq}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Chisq, μ, var)
    return Chisq(μ)
end

function dist_from_mean_var(::Type{Erlang}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Erlang, μ, var)
    θ = var/μ
    k = round(Int, μ^2/var)
    return Erlang(k, θ)
end

function dist_from_mean_var(::Type{Exponential}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Exponential, μ, var)
    return Exponential(μ)
end

function dist_from_mean_var(::Type{FDist}, μ::Number, var::Number)
    # exists_unique_dist_from_mean_var(FDist, μ, var)
    d_2 = 2*μ/(μ-1)
    d_1 = 2*μ^2*(d_2-2)/(var*(d_2-4)-2*μ^2)
    return FDist(d_1,d_2)
end