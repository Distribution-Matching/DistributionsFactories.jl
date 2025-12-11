using Distributions
include("exists_unique_dist_from_mean_var.jl")

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