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
    exists_unique_dist_from_mean_var(FDist, μ, var)
    v₂ = 2*μ/(μ-1)
    v₁ = 2*μ^2*(v₂-2)/(var*(v₂-4)-2*μ^2)
    return FDist(v₁,v₂)
end

function dist_from_mean_var(::Type{Frechet}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Frechet, μ, var)
end

function dist_from_mean_var(::Type{Gamma}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Gamma, μ, var)
    α = μ^2/var
    θ = var/μ
    return Gamma(α,θ)
end

function dist_from_mean_var(::Type{Gumbel}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Gumbel, μ, var)
    β = √(6*var/π^2)
    x̄ = μ-β*Base.MathConstants.γ
    return Gumbel(x̄,β)
end

function dist_from_mean_var(::Type{InverseGamma}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(InverseGamma, μ, var)
    α=(μ^2+2*var)/var
    β=μ*(α-1)
    return InverseGamma(α,β)
end

function dist_from_mean_var(::Type{Laplace}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Laplace, μ, var)
    x̄=μ
    b=√(var/2)
    return Laplace(x̄,b)
end

function dist_from_mean_var(::Type{Logistic}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Logistic, μ, var)
    x̄=μ
    s=√(3*var/π^2)
    return Logistic(x̄,s)
end

function dist_from_mean_var(::Type{LogNormal}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(LogNormal, μ, var)
    σ=log((var+μ^2)/μ^2)
    x̄=(2*log(μ)-σ^2)/2
    return LogNormal(x̄,σ)
end

function dist_from_mean_var(::Type{LogNormal}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Normal, μ, var)
    σ=√(var)
    return Normal(μ,σ)
end

function dist_from_mean_var(::Type{Rayleigh}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Rayleigh, μ, var)
    σ=√(2/π)*μ
    return Rayleigh(σ)
end

function dist_from_mean_var(::Type{TDist}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(TDist, μ, var)
    v=2*var/(var-1)
    return TDist(v)
end

function dist_from_mean_var(::Type{TriangularDist}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(TriangularDist, μ, var)
end

function dist_from_mean_var(::Type{Uniform}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Uniform, μ, var)
    diff = √(3*var)
    a = μ-diff
    b = μ+diff
    return Uniform(a,b)
end

# Doesnt work, need to test further
function dist_from_mean_var(::Type{Weibull}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Weibull, μ, var)
    f(k)=μ^2*(gamma(1+2/k)/gamma(1+1/k)^1-1)-var
    k = find_zero(f,μ^2-var)
    λ = μ/gamma(1+1/k)
    return Weibull()
end