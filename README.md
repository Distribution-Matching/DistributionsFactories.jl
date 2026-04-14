# DistributionsFactories.jl

A Julia package for creating probability distributions parameterised by moments (mean, variance, etc.) or quantiles. Builds on [Distributions.jl](https://github.com/JuliaStats/Distributions.jl).

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Ron-Ash/DistributionsFactories.jl.git")
```

## Quick start

```julia
using DistributionsFactories
using Distributions

# Construct a Gamma from mean and variance
d = dist_from_mean_var(Gamma, 2.5, 1.5)
# Gamma{Float64}(α=4.166666666666667, θ=0.6)

mean(d), var(d)
# (2.5, 1.5)
```

When no valid distribution exists for the given moments, an exception is thrown:

```julia
dist_from_mean_var(Exponential, 2.5, 1.5)
# ERROR: DomainError with "Exponential: the condition var = μ² is not satisfied"
```

Check feasibility before constructing:

```julia
exists_unique_dist_from_mean_var(Beta, 0.5, 0.1)
# true
```

## Moment-based construction

### `dist_from_mean_var(D, μ, var)`

The primary interface. Constructs a distribution of type `D` with the given mean and variance.

### Convenience wrappers

| Function | Specification |
|----------|--------------|
| `dist_from_mean_std(D, μ, σ)` | Mean and standard deviation |
| `dist_from_mean_cv(D, μ, cv)` | Mean and coefficient of variation |
| `dist_from_mean_scv(D, μ, scv)` | Mean and squared coefficient of variation |
| `dist_from_mean_second_moment(D, μ, m2)` | Mean and second moment E[X^2] |
| `dist_from_var(D, var)` | Variance only (1-parameter distributions) |
| `dist_from_std(D, σ)` | Standard deviation only (1-parameter distributions) |
| `dist_from_mean(D, μ)` | Mean only (1-parameter distributions) |

### Example: convenience wrappers

```julia
# Gamma from mean and coefficient of variation
d = dist_from_mean_cv(Gamma, 5.0, 0.5)
mean(d), std(d) / mean(d)
# (5.0, 0.5)

# Exponential from variance alone (mean is determined)
d = dist_from_var(Exponential, 4.0)
mean(d)
# 2.0
```

## Quantile-based construction

### Single quantile (1-parameter distributions)

```julia
# Exponential with median = 2.0
d = dist_from_median(Exponential, 2.0)
median(d)
# 2.0
```

### Two quantiles (2-parameter distributions)

```julia
# Normal from Q1 and Q3
d = dist_from_q1_q3(Normal, 10.0, 30.0)
quantile(d, 0.25), quantile(d, 0.75)
# (10.0, 30.0)

# Gamma from two arbitrary quantiles
d = dist_from_quantiles(Gamma, 0.1, 1.0, 0.9, 10.0)
quantile(d, 0.1), quantile(d, 0.9)
# (1.0, 10.0)
```

### Hybrid: mean + quantile

```julia
# Beta with mean 0.4 and median 0.35
d = dist_from_mean_median(Beta, 0.4, 0.35)
mean(d), median(d)
# (0.4, 0.35)
```

| Function | Specification |
|----------|--------------|
| `dist_from_quantile(D, p, q)` | Single quantile |
| `dist_from_quantiles(D, p1, q1, p2, q2)` | Two quantiles |
| `dist_from_median(D, m)` | Median |
| `dist_from_q1(D, q)` | First quartile |
| `dist_from_q3(D, q)` | Third quartile |
| `dist_from_q1_q3(D, q1, q3)` | First and third quartiles |
| `dist_from_median_iqr(D, median, iqr)` | Median and IQR |
| `dist_from_mean_quantile(D, μ, p, q)` | Mean and a quantile |
| `dist_from_mean_median(D, μ, median)` | Mean and median |

## Supported distributions

### Continuous (real line)

`Normal`, `TDist`, `Cauchy` (quantile only), `Logistic`, `Laplace`

### Continuous (positive support)

`Gamma`, `Erlang`, `Exponential`, `LogNormal`, `Weibull`, `Frechet`, `Gumbel`,
`Chi`, `Chisq`, `Rayleigh`, `FDist`, `InverseGamma`, `Pareto`, `FoldedNormal`

### Continuous (bounded)

`Beta`, `Uniform`, `SymTriangularDist`

### Truncated

`Truncated{Normal}`, `Truncated{Laplace}`, `Truncated{Logistic}` (pass an instance with bounds)

### Discrete

`Binomial`, `Poisson`, `NegativeBinomial`, `Geometric`, `DiscreteUniform`

### Not yet implemented

`TriangularDist`, `DiscreteTriangular`, `DiscreteSymmetricTriangular`, `TruncatedPoisson`

## Algorithms

Most distributions use direct closed-form formulas. For certain distributions, specialised algorithms are used:

- **Weibull / Frechet** -- numerical solution of the beta-ratio equation (Ron Ashri EVT approximation)
- **FoldedNormal** -- 2D Newton iteration for parent Normal parameters
- **Truncated Normal / Laplace / Logistic** -- 2D Newton iteration with quadrature-based moment computation

## Reference

Ashri, R., Moka, S., & Nazarathy, Y. (2026). DistributionsFactories.jl: Moment and Quantile Parameterised Probability Distributions. *Journal of Statistical Software* (submitted).
