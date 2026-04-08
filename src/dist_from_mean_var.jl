using Distributions
using Roots
using Polynomials
using SpecialFunctions



function dist_from_mean_var(::Type{Beta}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Beta, μ, var)
    S = (μ*(1-μ))/var-1
    α = μ*S
    β = (1-μ)*S
    return Beta(α,β)
end

function dist_from_mean_var(::Type{Uniform}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Uniform, μ, var)
    diff = √(3*var)
    a = μ-diff
    b = μ+diff
    return Uniform(a,b)
end


function dist_from_mean_var(::Type{Normal}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Normal, μ, var)
    return Normal(μ,√(var))
end


function dist_from_mean_var(::Type{TDist}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(TDist, μ, var)
    v=2*var/(var-1)
    return TDist(v)
end

function dist_from_mean_var(d::TDist, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(d, μ, var)
    ν = dof(d)
    base_var = ν / (ν - 2)
    σ = √(var / base_var)
    return μ + σ * d
end

function dist_from_mean_var(::Type{Cauchy}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Cauchy, μ, var)
end


function dist_from_mean_var(::Type{Logistic}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Logistic, μ, var)
    x̄=μ
    s=√(3*var/π^2)
    return Logistic(x̄,s)
end

function dist_from_mean_var(::Type{Laplace}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Laplace, μ, var)
    x̄=μ
    b=√(var/2)
    return Laplace(x̄,b)
end


function dist_from_mean_var(::Type{LogNormal}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(LogNormal, μ, var)
    σ=√(log(var/μ^2+1))
    x̄=log(μ^2/√(var+μ^2))
    return LogNormal(x̄,σ)
end


function dist_from_mean_var(::Type{Chisq}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Chisq, μ, var)
    return Chisq(μ)
end

function dist_from_mean_var(::Type{Exponential}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Exponential, μ, var)
    return Exponential(μ)
end

function dist_from_mean_var(::Type{Gamma}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Gamma, μ, var)
    α = μ^2/var
    θ = var/μ
    return Gamma(α,θ)
end

function dist_from_mean_var(::Type{Erlang}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Erlang, μ, var)
    k = round(Int, μ^2/var)
    θ = var/μ
    return Erlang(k, θ)
end


function ron_ashri_evt_approximation(μ::Number, var::Number, positiveSolution::Bool)
    μ_b  = BigFloat(μ)
    var_b = BigFloat(var)
    CV = √(var_b)/μ_b
    f(x) = x/beta(1/x,1/x)-(1+CV^2)/2
    if 0 < CV^2 < 1
        lowerBound = positiveSolution>0 ? 1/CV : min(-√(2π), -1/CV);
        upperBound = positiveSolution>0 ? (CV^2+1)/(2CV^2) : -2(1+CV^2)/(CV^2);
        # println(f(lowerBound), f(upperBound))
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

function ron_ashri_weibull_approximation(μ::Number, var::Number)
    k = ron_ashri_evt_approximation(μ,var, true)
    λ = μ/gamma(1+1/k)
    return Weibull(k, λ)
end

function ron_ashri_frechet_approximation(μ::Number, var::Number)
    α=-1*ron_ashri_evt_approximation(μ,var, false)
    s = μ/gamma(1-1/α)
    return Frechet(α, s)
end

function dist_from_mean_var(::Type{Frechet}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Frechet, μ, var)
    return ron_ashri_frechet_approximation(μ,var)
end

function dist_from_mean_var(::Type{Weibull}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Weibull, μ, var)
    return ron_ashri_weibull_approximation(μ,var)
end

function dist_from_mean_var(::Type{Gumbel}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Gumbel, μ, var)
    β = √(6*var/π^2)
    x̄ = μ-β*Base.MathConstants.γ
    return Gumbel(x̄,β)
end


function dist_from_mean_var(::Type{Chi}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Chi, μ, var)
    return Chi(μ^2+var)
end

function dist_from_mean_var(::Type{Rayleigh}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Rayleigh, μ, var)
    σ=√(2/π)*μ
    return Rayleigh(σ)
end


function dist_from_mean_var(::Type{FDist}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(FDist, μ, var)
    v₂ = 2*μ/(μ-1)
    v₁ = 2*μ^2*(v₂-2)/(var*(v₂-4)-2*μ^2)
    return FDist(v₁,v₂)
end


function dist_from_mean_var(::Type{InverseGamma}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(InverseGamma, μ, var)
    α=(μ^2+2*var)/var
    β=μ*(α-1)
    return InverseGamma(α,β)
end


function dist_from_mean_var(::Type{Binomial}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Binomial, μ, var)
    p=1-var/μ
    n=round(Int, μ/p)
    return Binomial(n,p)
end


function dist_from_mean_var(::Type{Poisson}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(Poisson, μ, var)
    return Poisson(μ)
end


function dist_from_mean_var(::Type{NegativeBinomial}, μ::Number, var::Number)
    exists_unique_dist_from_mean_var(NegativeBinomial, μ, var)
    p=μ/var
    r=μ^2/(var-μ)
    return NegativeBinomial(r,p)
end


# --- Not yet implemented ---

function dist_from_mean_var(::Type{Pareto}, μ::Number, var::Number)
    throw(ErrorException("Pareto: dist_from_mean_var not yet implemented"))
end

function dist_from_mean_var(::Type{FoldedNormal}, μ::Number, var::Number)
    throw(ErrorException("FoldedNormal: dist_from_mean_var not yet implemented"))
end

function dist_from_mean_var(::Type{Geometric}, μ::Number, var::Number)
    throw(ErrorException("Geometric: dist_from_mean_var not yet implemented"))
end

function dist_from_mean_var(d::Truncated{<:Normal}, μ::Number, var::Number)
    throw(ErrorException("Truncated Normal: dist_from_mean_var not yet implemented"))
end

function dist_from_mean_var(d::Truncated{<:Laplace}, μ::Number, var::Number)
    throw(ErrorException("Truncated Laplace: dist_from_mean_var not yet implemented"))
end

function dist_from_mean_var(d::Truncated{<:Logistic}, μ::Number, var::Number)
    throw(ErrorException("Truncated Logistic: dist_from_mean_var not yet implemented"))
end

function dist_from_mean_var(::Type{TriangularDist}, μ::Number, var::Number)
    throw(ErrorException("TriangularDist: dist_from_mean_var not yet implemented"))
end

function dist_from_mean_var(::Type{SymTriangularDist}, μ::Number, var::Number)
    throw(ErrorException("SymTriangularDist: dist_from_mean_var not yet implemented"))
end

function dist_from_mean_var(::Type{DiscreteTriangular}, μ::Number, var::Number)
    throw(ErrorException("DiscreteTriangular: dist_from_mean_var not yet implemented"))
end

function dist_from_mean_var(::Type{DiscreteSymmetricTriangular}, μ::Number, var::Number)
    throw(ErrorException("DiscreteSymmetricTriangular: dist_from_mean_var not yet implemented"))
end

function dist_from_mean_var(::Type{TruncatedPoisson}, μ::Number, var::Number)
    throw(ErrorException("TruncatedPoisson: dist_from_mean_var not yet implemented"))
end

function dist_from_mean_var(::Type{DiscreteUniform}, μ::Number, var::Number)
    throw(ErrorException("DiscreteUniform: dist_from_mean_var not yet implemented"))
end