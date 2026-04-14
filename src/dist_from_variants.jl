# Variants that convert to mean+var and delegate

"""
    dist_from_mean_std(D, μ̄, σ̄)

Construct distribution `D` from mean `μ̄` and standard deviation `σ̄`.
Delegates to `dist_from_mean_var(D, μ̄, σ̄²)`.
"""
dist_from_mean_std(D, μ̄::Number, σ̄::Number) = dist_from_mean_var(D, μ̄, σ̄^2)

"""
    dist_from_mean_cv(D, μ̄, cv)

Construct distribution `D` from mean `μ̄` and coefficient of variation `cv = σ̄/μ̄`.
Delegates to `dist_from_mean_var(D, μ̄, (cv·μ̄)²)`.
"""
dist_from_mean_cv(D, μ̄::Number, cv::Number) = dist_from_mean_var(D, μ̄, (cv * μ̄)^2)

"""
    dist_from_mean_scv(D, μ̄, scv)

Construct distribution `D` from mean `μ̄` and squared coefficient of variation `scv = σ̄²/μ̄²`.
Delegates to `dist_from_mean_var(D, μ̄, scv·μ̄²)`.
"""
dist_from_mean_scv(D, μ̄::Number, scv::Number) = dist_from_mean_var(D, μ̄, scv * μ̄^2)

"""
    dist_from_mean_second_moment(D, μ̄, m2)

Construct distribution `D` from mean `μ̄` and second moment `m2 = E[X²]`.
Delegates to `dist_from_mean_var(D, μ̄, m2 - μ̄²)`.
"""
dist_from_mean_second_moment(D, μ̄::Number, m2::Number) = dist_from_mean_var(D, μ̄, m2 - μ̄^2)

"""
    exists_unique_dist_from_mean_std(D, μ̄, σ̄)

Check feasibility for `dist_from_mean_std`. See [`exists_unique_dist_from_mean_var`](@ref).
"""
exists_unique_dist_from_mean_std(D, μ̄::Number, σ̄::Number) = exists_unique_dist_from_mean_var(D, μ̄, σ̄^2)

"""
    exists_unique_dist_from_mean_cv(D, μ̄, cv)

Check feasibility for `dist_from_mean_cv`. See [`exists_unique_dist_from_mean_var`](@ref).
"""
exists_unique_dist_from_mean_cv(D, μ̄::Number, cv::Number) = exists_unique_dist_from_mean_var(D, μ̄, (cv * μ̄)^2)

"""
    exists_unique_dist_from_mean_scv(D, μ̄, scv)

Check feasibility for `dist_from_mean_scv`. See [`exists_unique_dist_from_mean_var`](@ref).
"""
exists_unique_dist_from_mean_scv(D, μ̄::Number, scv::Number) = exists_unique_dist_from_mean_var(D, μ̄, scv * μ̄^2)

"""
    exists_unique_dist_from_mean_second_moment(D, μ̄, m2)

Check feasibility for `dist_from_mean_second_moment`. See [`exists_unique_dist_from_mean_var`](@ref).
"""
exists_unique_dist_from_mean_second_moment(D, μ̄::Number, m2::Number) = exists_unique_dist_from_mean_var(D, μ̄, m2 - μ̄^2)

# Variants that only take a dispersion measure (for 1-parameter distributions where mean is determined)

"""
    dist_from_var(D, σ̄²)

Construct a 1-parameter distribution `D` from variance alone. Only supported for
distributions where the mean is uniquely determined by the variance:
`Exponential`, `Poisson`, `Chisq`, `Rayleigh`, `Geometric`.
"""
dist_from_var(D, σ̄²::Number) = dist_from_mean_var(D, _mean_from_var(D, σ̄²), σ̄²)

"""
    dist_from_std(D, σ̄)

Construct a 1-parameter distribution `D` from standard deviation alone.
Delegates to `dist_from_var(D, σ̄²)`.
"""
dist_from_std(D, σ̄::Number) = dist_from_var(D, σ̄^2)

# mean-from-var helpers for 1-parameter distributions where var determines mean
_mean_from_var(::Type{Exponential}, σ̄²::Number) = √σ̄²                         # σ̄² = μ̄²
_mean_from_var(::Type{Poisson}, σ̄²::Number) = σ̄²                               # σ̄² = μ̄
_mean_from_var(::Type{Chisq}, σ̄²::Number) = σ̄² / 2                             # σ̄² = 2μ̄
_mean_from_var(::Type{Rayleigh}, σ̄²::Number) = √(σ̄² * π / (4 - π))            # σ̄² = μ̄²(4-π)/π
_mean_from_var(::Type{Geometric}, σ̄²::Number) = (-1 + √(1 + 4σ̄²)) / 2        # σ̄² = μ̄(1+μ̄)
_mean_from_var(::Type{D}, σ̄²::Number) where {D<:Distribution} =
    throw(ErrorException("$D: dist_from_var not supported (mean is not determined by variance alone)"))

# Quantile convenience wrappers (delegate to dist_from_quantile / dist_from_quantiles)

"""
    dist_from_median(D, m)

Construct distribution `D` from its median. Equivalent to `dist_from_quantile(D, 0.5, m)`.
"""
dist_from_median(D, m::Number) = dist_from_quantile(D, 0.5, m)

"""
    dist_from_q1(D, q)

Construct distribution `D` from its first quartile. Equivalent to `dist_from_quantile(D, 0.25, q)`.
"""
dist_from_q1(D, q::Number) = dist_from_quantile(D, 0.25, q)

"""
    dist_from_q3(D, q)

Construct distribution `D` from its third quartile. Equivalent to `dist_from_quantile(D, 0.75, q)`.
"""
dist_from_q3(D, q::Number) = dist_from_quantile(D, 0.75, q)

"""
    dist_from_q1_q3(D, q1, q3)

Construct distribution `D` from its first and third quartiles.
"""
dist_from_q1_q3(D, q1::Number, q3::Number) = dist_from_quantiles(D, 0.25, q1, 0.75, q3)

"""
    dist_from_median_iqr(D, median, iqr)

Construct distribution `D` from its median and interquartile range.
"""
dist_from_median_iqr(D, median::Number, iqr::Number) = dist_from_quantiles(D, 0.25, median - iqr / 2, 0.75, median + iqr / 2)

"""
    dist_from_mean_median(D, μ̄, median)

Construct distribution `D` from its mean `μ̄` and median.
"""
dist_from_mean_median(D, μ̄::Number, median::Number) = dist_from_mean_quantile(D, μ̄, 0.5, median)
