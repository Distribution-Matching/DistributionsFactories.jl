# DistributionsFactories.jl

A Julia package for creating probability distributions parameterised by moments (mean, variance, etc.) or quantiles, on standard or arbitrary supports. Builds on [Distributions.jl](https://github.com/JuliaStats/Distributions.jl).

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
mean(d), var(d)
# (2.5, 1.5)

# What distributions are feasible for these moments?
available_distributions(0..Inf, mean=2.5, var=1.5)
# [Gamma, Erlang, LogNormal, Weibull, Frechet, ...]
```

When no valid distribution exists for the given moments, an exception is thrown:

```julia
dist_from_mean_var(Exponential, 2.5, 1.5)
# ERROR: DomainError with "Exponential: the condition σ̄² = μ̄² is not satisfied"
```

Check feasibility before constructing:

```julia
exists_dist_from_mean_var(Beta, 0.5, 0.1)
# true
```

## Moment-based construction

### `dist_from_mean_var(D, μ̄, σ̄²)`

The primary interface. Constructs a distribution of type `D` with mean `μ̄` and variance `σ̄²`.

### Convenience wrappers

| Function | Specification |
|----------|--------------|
| `dist_from_mean_std(D, μ̄, σ̄)` | Mean and standard deviation |
| `dist_from_mean_cv(D, μ̄, cv)` | Mean and coefficient of variation |
| `dist_from_mean_scv(D, μ̄, scv)` | Mean and squared coefficient of variation |
| `dist_from_mean_second_moment(D, μ̄, m2)` | Mean and second moment E[X^2] |
| `dist_from_var(D, σ̄²)` | Variance only (1-parameter distributions) |
| `dist_from_std(D, σ̄)` | Standard deviation only (1-parameter distributions) |
| `dist_from_mean(D, μ̄)` | Mean only (1-parameter distributions) |

```julia
# Gamma from mean and coefficient of variation
d = dist_from_mean_cv(Gamma, 5.0, 0.5)

# Exponential from variance alone (mean is determined)
d = dist_from_var(Exponential, 4.0)
```

## Distributions on arbitrary supports

By default, distributions are constructed on their natural support (e.g. Beta on [0,1], Gamma on [0,∞)). Pass a support interval as a fourth argument to place the distribution on a different domain.

The function automatically determines whether to use an **affine transform** or **truncation**:

- **Affine**: when the requested support has the same shape as the natural one (e.g. Beta on [2,7], Gamma on [3,∞)). Returns a `LocationScale` wrapper.
- **Truncation**: when the requested support restricts the natural one (e.g. Normal on [0,1], Gamma on [0,10]). Returns a `Truncated` wrapper.

### Continuous examples

```julia
# Beta scaled from [0,1] to [2,7]
d = dist_from_mean_var(Beta, 3.5, 0.5, 2..7)
minimum(d), maximum(d)   # (2.0, 7.0)

# Gamma shifted from [0,∞) to [3,∞)
d = dist_from_mean_var(Gamma, 8.0, 3.0, 3..Inf)
minimum(d)   # 3.0

# Flipped Gamma on (-∞,10]
d = dist_from_mean_var(Gamma, 5.0, 3.0, -Inf..10)
maximum(d)   # 10.0

# Normal truncated to [0,1]
d = dist_from_mean_var(Normal, 0.5, 0.04, 0..1)
minimum(d), maximum(d)   # (0.0, 1.0)
```

### Discrete examples

```julia
# Binomial shifted from {0,...,5} to {10,...,15}
d = dist_from_mean_var(Binomial, 12.0, 1.2, 10:15)
minimum(d), maximum(d)   # (10, 15)
```

Supports are specified using [IntervalSets.jl](https://github.com/JuliaMath/IntervalSets.jl) syntax (`a..b`) or Distributions.jl's `RealInterval` (e.g. `support(some_dist)`), and integer ranges (`a:b`) for discrete distributions.

## Quantile-based construction

### Single quantile (1-parameter distributions)

```julia
# Exponential with median = 2.0
d = dist_from_median(Exponential, 2.0)
```

### Two quantiles (2-parameter distributions)

```julia
# Normal from Q1 and Q3
d = dist_from_q1_q3(Normal, 10.0, 30.0)

# Gamma from two arbitrary quantiles
d = dist_from_quantiles(Gamma, 0.1, 1.0, 0.9, 10.0)
```

### Hybrid: mean + quantile

```julia
# Beta with mean 0.4 and median 0.35
d = dist_from_mean_median(Beta, 0.4, 0.35)
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
| `dist_from_mean_quantile(D, μ̄, p, q)` | Mean and a quantile |
| `dist_from_mean_median(D, μ̄, median)` | Mean and median |

## Discovering available distributions

Use `available_distributions` to find which distribution types are available for a given support and/or moment specification.

```julia
# What distributions live on [0,∞)?
available_distributions(0..Inf)
# [Gamma, Erlang, Exponential, LogNormal, Weibull, ...]

# Which of those are feasible for μ̄=5, σ̄²=3?
available_distributions(0..Inf, mean=5.0, var=3.0)
# [Gamma, Erlang, LogNormal, Weibull, Frechet, ...]

# Use the support of an existing distribution
available_distributions(support(Beta(2, 3)), mean=0.5, std=0.2)

# Search across all supports
available_distributions(mean=5.0, cv=1.0)
```

Moment specification via named arguments: `mean`, `var`, `std`, `cv`, `scv`, `second_moment`.

## Supported distributions

### Continuous (real line)

`Normal`, `TDist`, `Cauchy` (quantile only), `Logistic`, `Laplace`, `Gumbel`, `SymTriangularDist`

### Continuous (positive support)

`Gamma`, `Erlang`, `Exponential`, `LogNormal`, `Weibull`, `Frechet`,
`Chi`, `Chisq`, `Rayleigh`, `FDist`, `InverseGamma`, `Pareto`, `FoldedNormal`

### Continuous (bounded)

`Beta`, `Uniform`

### Discrete

`Binomial`, `Poisson`, `NegativeBinomial`, `Geometric`, `DiscreteUniform`

### Not yet implemented

`TriangularDist`, `DiscreteTriangular`, `DiscreteSymmetricTriangular`, `TruncatedPoisson`

## Algorithms

Most distributions use direct closed-form formulas. For certain distributions, specialised algorithms are used:

| Distribution | Method |
|---|---|
| Beta, Normal, Gamma, Logistic, Laplace, ... | Direct formula |
| Erlang, Binomial, DiscreteUniform | Direct formula (with integer rounding) |
| Weibull, Frechet | Numerical root-finding (beta-ratio equation) |
| FoldedNormal | 2D Newton iteration |
| Truncated Normal, Laplace, Logistic, Gamma | 2D Newton iteration with quadrature |
| Quantile matching (Gamma, Beta) | Numerical root-finding |
| Quantile matching (Normal, Laplace, Logistic, Cauchy, Gumbel) | Direct formula (location-scale) |
| Arbitrary support (affine) | Moment transform + `LocationScale` wrapper |

## Reference

Ashri, R., Moka, S., & Nazarathy, Y. (2026). DistributionsFactories.jl: Moment and Quantile Parameterised Probability Distributions. *Journal of Statistical Software* (submitted).
