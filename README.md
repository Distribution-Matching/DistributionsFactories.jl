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
# ERROR: DomainError: "Exponential: the condition σ̄² = μ̄² is not satisfied"
```

Check feasibility before constructing:

```julia
exists_dist_from_mean_var(Beta, 0.5, 0.1)
# true
```

## Distributions covered

The following table lists all distributions covered by the package, matching Table 1 of the accompanying paper. The **Status** column indicates implementation state: **Yes** = fully implemented, **No** = not yet implemented. All distributions marked "Yes" support `dist_from_mean_var` construction and feasibility checking via `exists_dist_from_mean_var`.

### Continuous distributions

| # | Distribution | Support | Free params | Implemented | Method | Notes |
|---|---|---|---|---|---|---|
| 1 | Student's T | (-∞, ∞) | 3 | Yes | Direct formula | Type `TDist`; instance form for location-scale |
| 2 | Normal | (-∞, ∞) | 2 | Yes | Direct formula | Special case of Student's T (ν → ∞) |
| 3 | Cauchy | (-∞, ∞) | 2 | Quantile only | — | Moments undefined; use `dist_from_quantiles` |
| 4 | Laplace | (-∞, ∞) | 2 | Yes | Direct formula | |
| 5 | Logistic | (-∞, ∞) | 2 | Yes | Direct formula | |
| 6 | Gamma | [0, ∞) | 2 | Yes | Direct formula | |
| 7 | Erlang | [0, ∞) | 2 | Yes | Direct formula | Special case of Gamma (integer shape) |
| 8 | Exponential | [0, ∞) | 1 | Yes | Direct formula | Special case of Gamma (α = 1); 1 DOF |
| 9 | Chi-squared | [0, ∞) | 1 | Yes | Direct formula | Special case of Gamma; 1 DOF |
| 10 | Log-normal | [0, ∞) | 2 | Yes | Direct formula | |
| 11 | Weibull | [0, ∞) | 2 | Yes | Numerical (root-finding) | Beta-ratio equation |
| 12 | Frechet | [0, ∞) | 2 | Yes | Numerical (root-finding) | Beta-ratio equation |
| 13 | Gumbel | (-∞, ∞) | 2 | Yes | Direct formula | |
| 14 | F distribution | [0, ∞) | 2 | Yes | Direct formula | Requires 1 < μ̄ < 2 |
| 15 | Inverse Gamma | [0, ∞) | 2 | Yes | Direct formula | |
| 16 | Half-truncated Normal | [0, ∞) | 2 | No | | |
| 17 | Half-truncated Laplace | [0, ∞) | 2 | No | | |
| 18 | Half-truncated Logistic | [0, ∞) | 2 | No | | |
| 19 | Folded Normal | [0, ∞) | 1 | Yes | Numerical (2D Newton) | Returns parent Normal |
| 20 | Pareto | [0, ∞) | 1 | Yes | Direct formula | |
| 21 | Chi | [0, ∞) | 1 | Yes | Direct formula | |
| 22 | Rayleigh | [0, ∞) | 1 | Yes | Direct formula | Special case of Chi (ν = 2); 1 DOF |
| 23 | Triangular | [0, 1] | 1 | No | | |
| 24 | Symmetric Triangular | [0, 1] | 0 | Yes | Direct formula | Type `SymTriangularDist`; adjusts support |
| 25 | Beta | [0, 1] | 2 | Yes | Direct formula | |
| 26 | Uniform | [0, 1] | 0 | Yes | Direct formula | Special case of Beta; adjusts support |
| 27 | Doubly-truncated Normal | [0, 1] | 2 | No | | |
| 28 | Doubly-truncated Laplace | [0, 1] | 2 | No | | |
| 29 | Doubly-truncated Logistic | [0, 1] | 2 | No | | |

### Discrete distributions

| # | Distribution | Support | Free params | Implemented | Method | Notes |
|---|---|---|---|---|---|---|
| 30 | Discrete Triangular | {0, …, n} | 1 | No | | |
| 31 | Discrete Sym. Triangular | {0, …, n} | 0 | No | | |
| 32 | Binomial | {0, …, n} | 1 | Yes | Direct formula | |
| 33 | Truncated Poisson | {0, …, n} | 1 | No | | |
| 34 | Discrete Uniform | {0, …, n} | 0 | Yes | Direct formula | Adjusts support |
| 35 | Negative Binomial | {0, 1, 2, …} | 2 | Yes | Direct formula | |
| 36 | Geometric | {0, 1, 2, …} | 1 | Yes | Direct formula | Special case of Neg. Binomial; 1 DOF |
| 37 | Poisson | {0, 1, 2, …} | 1 | Yes | Direct formula | 1 DOF |

### Experimental / work in progress

The following have initial code but are not fully validated:

- **Truncated Normal on [a,b]** — 2D Newton solver with quadrature. Passes basic tests.
- **Truncated Gamma on [a,b]** — same solver. Passes basic tests.
- **Truncated Laplace on [a,b]** — same solver. Not yet tested end-to-end.
- **Truncated Logistic on [a,b]** — same solver. Not yet tested end-to-end.

## Moment-based construction

### `dist_from_mean_var(D, μ̄, σ̄²)`

The primary interface. Constructs a distribution of type `D` with mean `μ̄` and variance `σ̄²`.

### Convenience wrappers

| Function | Specification |
|----------|--------------|
| `dist_from_mean_std(D, μ̄, σ̄)` | Mean and standard deviation |
| `dist_from_mean_cv(D, μ̄, cv)` | Mean and coefficient of variation |
| `dist_from_mean_scv(D, μ̄, scv)` | Mean and squared coefficient of variation |
| `dist_from_mean_second_moment(D, μ̄, m2)` | Mean and second moment E[X²] |
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

By default, distributions are constructed on their natural support (e.g. Beta on [0,1], Gamma on [0,∞)). Use `dist_from_mean_var_on_support` to place a distribution on a different domain.

The function automatically determines whether to use an **affine transform** or **truncation**:

- **Affine**: when the requested support has the same shape as the natural one (e.g. Beta on [2,7], Gamma on [3,∞)). Returns a `LocationScale` wrapper.
- **Truncation**: when the requested support restricts the natural one (e.g. Normal on [0,1], Gamma on [0,10]). Returns a `Truncated` wrapper.

### Continuous examples

```julia
# Beta scaled from [0,1] to [2,7]
d = dist_from_mean_var_on_support(Beta, 3.5, 0.5, support=2..7)
minimum(d), maximum(d)   # (2.0, 7.0)

# Gamma shifted from [0,∞) to [3,∞)
d = dist_from_mean_var_on_support(Gamma, 8.0, 3.0, support=3..Inf)
minimum(d)   # 3.0

# Flipped Gamma on (-∞,10]
d = dist_from_mean_var_on_support(Gamma, 5.0, 3.0, support=-Inf..10)
maximum(d)   # 10.0

# Normal truncated to [0,1] (experimental)
d = dist_from_mean_var_on_support(Normal, 0.5, 0.04, support=0..1)
minimum(d), maximum(d)   # (0.0, 1.0)
```

### Discrete examples

```julia
# Binomial shifted from {0,...,5} to {10,...,15}
d = dist_from_mean_var_on_support(Binomial, 12.0, 1.2, support=10:15)
minimum(d), maximum(d)   # (10, 15)
```

Supports are specified using [IntervalSets.jl](https://github.com/JuliaMath/IntervalSets.jl) syntax (`a..b`), Distributions.jl's `RealInterval` (e.g. `support(some_dist)`), or integer ranges (`a:b`) for discrete distributions.

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

# Which of those are feasible for μ̄=5, σ̄²=3?
available_distributions(0..Inf, mean=5.0, var=3.0)

# Use the support of an existing distribution
available_distributions(support(Beta(2, 3)), mean=0.5, std=0.2)

# Search across all supports
available_distributions(mean=5.0, cv=1.0)
```

Moment specification via named arguments: `mean`, `var`, `std`, `cv`, `scv`, `second_moment`.

## Algorithms

| Distribution | Method |
|---|---|
| Normal, Beta, Gamma, Logistic, Laplace, Gumbel, LogNormal, ... | Direct formula |
| Erlang, Binomial, DiscreteUniform | Direct formula (with integer rounding) |
| Weibull, Frechet | Numerical root-finding (beta-ratio equation) |
| Folded Normal | Numerical (2D Newton iteration) |
| Truncated Normal, Gamma (experimental) | Numerical (2D Newton + quadrature) |
| Quantile matching (Gamma, Beta) | Numerical root-finding |
| Quantile matching (Normal, Laplace, Logistic, Cauchy, Gumbel) | Direct formula (location-scale) |
| Arbitrary support (affine) | Moment transform + `LocationScale` wrapper |

## Reference

Ashri, R., Moka, S., & Nazarathy, Y. (2026). DistributionsFactories.jl: Moment and Quantile Parameterised Probability Distributions. *Journal of Statistical Software* (submitted).
