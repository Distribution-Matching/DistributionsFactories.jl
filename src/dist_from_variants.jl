# Variants that convert to mean+var and delegate

dist_from_mean_std(D, μ::Number, σ::Number) = dist_from_mean_var(D, μ, σ^2)
dist_from_mean_cv(D, μ::Number, cv::Number) = dist_from_mean_var(D, μ, (cv * μ)^2)
dist_from_mean_scv(D, μ::Number, scv::Number) = dist_from_mean_var(D, μ, scv * μ^2)
dist_from_mean_second_moment(D, μ::Number, m2::Number) = dist_from_mean_var(D, μ, m2 - μ^2)

exists_unique_dist_from_mean_std(D, μ::Number, σ::Number) = exists_unique_dist_from_mean_var(D, μ, σ^2)
exists_unique_dist_from_mean_cv(D, μ::Number, cv::Number) = exists_unique_dist_from_mean_var(D, μ, (cv * μ)^2)
exists_unique_dist_from_mean_scv(D, μ::Number, scv::Number) = exists_unique_dist_from_mean_var(D, μ, scv * μ^2)
exists_unique_dist_from_mean_second_moment(D, μ::Number, m2::Number) = exists_unique_dist_from_mean_var(D, μ, m2 - μ^2)

# Variants that only take a dispersion measure (for 1-parameter distributions where mean is determined)

dist_from_var(D, var::Number) = dist_from_mean_var(D, _mean_from_var(D, var), var)
dist_from_std(D, σ::Number) = dist_from_var(D, σ^2)

# mean-from-var helpers for 1-parameter distributions where var determines mean
_mean_from_var(::Type{Exponential}, var::Number) = √var                       # var = μ²
_mean_from_var(::Type{Poisson}, var::Number) = var                             # var = μ
_mean_from_var(::Type{Chisq}, var::Number) = var / 2                           # var = 2μ
_mean_from_var(::Type{Rayleigh}, var::Number) = √(var * π / (4 - π))          # var = μ²(4-π)/π
_mean_from_var(::Type{D}, var::Number) where {D<:Distribution} =
    throw(ErrorException("$D: dist_from_var not supported (mean is not determined by variance alone)"))

# Quantile convenience wrappers (delegate to dist_from_quantile / dist_from_quantiles)

dist_from_median(D, m::Number) = dist_from_quantile(D, 0.5, m)
dist_from_q1(D, q::Number) = dist_from_quantile(D, 0.25, q)
dist_from_q3(D, q::Number) = dist_from_quantile(D, 0.75, q)
dist_from_q1_q3(D, q1::Number, q3::Number) = dist_from_quantiles(D, 0.25, q1, 0.75, q3)
dist_from_median_iqr(D, median::Number, iqr::Number) = dist_from_quantiles(D, 0.25, median - iqr / 2, 0.75, median + iqr / 2)
dist_from_mean_median(D, μ::Number, median::Number) = dist_from_mean_quantile(D, μ, 0.5, median)
