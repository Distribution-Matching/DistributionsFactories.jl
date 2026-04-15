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
mean(d), var(d)   # (2.5, 1.5)

# Or use the @dist macro for partial specification
d = dist_from_mean(@dist(Gamma(3.0, _)), 5.0)   # fix α=3, solve θ from mean

# What distributions are feasible for these moments on [0,∞)?
available_distributions(0..Inf, mean=2.5, var=1.5)
# [Gamma, Erlang, LogNormal, Weibull, Frechet, ...]

# Check feasibility before constructing
exists_dist_from_mean_var(Beta, 0.5, 0.1)   # true

# When no valid distribution exists, an exception is thrown
dist_from_mean_var(Exponential, 2.5, 1.5)
# ERROR: DomainError: "Exponential: the condition σ̄² = μ̄² is not satisfied"
```

## Extended examples

```julia
using DistributionsFactories, Distributions

# --- Moment-based construction ---

# Gamma from mean and coefficient of variation
d = dist_from_mean_cv(Gamma, 10.0, 0.5)
mean(d), std(d) / mean(d)   # (10.0, 0.5)

# Beta on [2, 7] with target mean and variance
d = dist_from_mean_var_on_support(Beta, 3.5, 0.5, support=2..7)
mean(d), var(d), minimum(d), maximum(d)   # (3.5, 0.5, 2.0, 7.0)

# Gamma shifted to [3, ∞)
d = dist_from_mean_var_on_support(Gamma, 8.0, 3.0, support=3..Inf)
mean(d), minimum(d)   # (8.0, 3.0)

# --- Partial specification with @dist ---

# Fix shape α=3, solve scale θ from mean
d = dist_from_mean(@dist(Gamma(3.0, _)), 5.0)
params(d)   # (3.0, 1.667)

# Fix location μ=0, solve scale σ from variance
d = dist_from_mean_var(@dist(Normal(0.0, _)), 0.0, 4.0)
params(d)   # (0.0, 2.0)

# TDist with fixed degrees of freedom, arbitrary mean and variance
d = dist_from_mean_var(@dist(TDist(7)), 5.0, 2.0)
mean(d), var(d)   # (5.0, 2.0)

# --- Quantile-based construction ---

# Normal from first and third quartiles
d = dist_from_q1_q3(Normal, 10.0, 30.0)
quantile(d, 0.25), quantile(d, 0.75)   # (10.0, 30.0)

# Beta with mean 0.4 and median 0.35
d = dist_from_mean_median(Beta, 0.4, 0.35)
mean(d), median(d)   # (0.4, 0.35)

# --- Discovery ---

# What distributions can represent a positive r.v. with mean=5, var=25?
available_distributions(0..Inf, mean=5.0, var=25.0)
# [Gamma, Erlang, Exponential, LogNormal, Weibull, ...]

# Search across all supports
available_distributions(mean=5.0, cv=1.0)
```

## The `@dist` macro

The `@dist` macro creates distribution specifications where `_` marks parameters
to be solved from moment constraints. The result is a `PartialDist` that can be
passed to any `dist_from_*` function.

```julia
# Create a partial spec — α is fixed, θ is to be solved
spec = @dist Gamma(3.0, _)
fixed_params(spec)   # (α = 3.0,)
free_params(spec)    # (:θ,)
```

Use it with any `dist_from_*` function:

```julia
# Solve the free parameter from different moment specifications
dist_from_mean(@dist(Gamma(3.0, _)), 5.0)           # solve θ from mean
dist_from_var(@dist(Logistic(2.0, _)), 22.3)         # solve θ from variance
dist_from_std(@dist(Logistic(2.0, _)), 4.0)          # solve θ from std
dist_from_mean_var(@dist(Normal(0.0, _)), 0.0, 4.0)  # solve σ from mean+var
dist_from_mean_cv(@dist(Gamma(4.0, _)), 5.0, 0.5)    # solve θ from mean+cv
dist_from_mean_std(@dist(Normal(_, 2.0)), 3.0, 2.0)  # solve μ from mean+std
dist_from_mean(@dist(Beta(2.0, _)), 0.4)             # solve β from mean
```

The macro also handles bare types and full instances:

```julia
# Bare type — same as passing the type directly
dist_from_mean_var(@dist(Gamma), 5.0, 3.0)

# Full instance — wraps in LocationScale if moments differ
dist_from_mean_var(@dist(TDist(7)), 5.0, 2.0)
```

## Distributions covered

The following table lists all distributions covered by the package. The **Status** column indicates implementation state: **Yes** = fully implemented, **No** = not yet implemented. All distributions marked "Yes" support `dist_from_mean_var` construction and feasibility checking via `exists_dist_from_mean_var`.

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

The primary interface. Constructs a distribution of type `D` with mean `μ̄` and variance `σ̄²`. The first argument `D` is usually a type (e.g. `Gamma`, `Beta`), but can also be a distribution **instance** when some parameters are fixed:

```julia
# Type: standard TDist with μ=0, determines ν from variance
d = dist_from_mean_var(TDist, 0.0, 3.0)

# Instance: TDist with fixed ν=7, arbitrary mean and variance via location-scale
d = dist_from_mean_var(TDist(7), 5.0, 2.0)
# Returns a LocationScale wrapping TDist(7)
mean(d), var(d)   # (5.0, 2.0)
```

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

## Distributions on arbitrary supports

By default, distributions are constructed on their natural support (e.g. Beta on [0,1], Gamma on [0,∞)). Use `dist_from_mean_var_on_support` to place a distribution on a different domain.

The function automatically determines whether to use an **affine transform** or **truncation**:

- **Affine**: when the requested support has the same shape as the natural one (e.g. Beta on [2,7], Gamma on [3,∞)). Returns a `LocationScale` wrapper.
- **Truncation**: when the requested support restricts the natural one (e.g. Normal on [0,1], Gamma on [0,10]). Returns a `Truncated` wrapper.

```julia
# Beta scaled from [0,1] to [2,7]
d = dist_from_mean_var_on_support(Beta, 3.5, 0.5, support=2..7)

# Gamma shifted from [0,∞) to [3,∞)
d = dist_from_mean_var_on_support(Gamma, 8.0, 3.0, support=3..Inf)

# Flipped Gamma on (-∞,10]
d = dist_from_mean_var_on_support(Gamma, 5.0, 3.0, support=-Inf..10)

# Binomial shifted from {0,...,5} to {10,...,15}
d = dist_from_mean_var_on_support(Binomial, 12.0, 1.2, support=10:15)
```

Supports are specified using [IntervalSets.jl](https://github.com/JuliaMath/IntervalSets.jl) syntax (`a..b`), Distributions.jl's `RealInterval` (e.g. `support(some_dist)`), or integer ranges (`a:b`) for discrete distributions.

## Quantile-based construction

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
available_distributions(0..Inf)                          # by support
available_distributions(0..Inf, mean=5.0, var=3.0)       # filtered by moments
available_distributions(support(Beta(2, 3)), mean=0.5)   # from existing distribution
available_distributions(mean=5.0, cv=1.0)                # across all supports
```

Moment keywords: `mean`, `var`, `std`, `cv`, `scv`, `second_moment`.

## Algorithms

| Distribution | Method |
|---|---|
| Normal, Beta, Gamma, Logistic, Laplace, Gumbel, LogNormal, ... | Direct formula |
| Erlang, Binomial, DiscreteUniform | Direct formula (with integer rounding) |
| Weibull, Frechet | Numerical root-finding (beta-ratio equation) |
| Folded Normal | Numerical (2D Newton iteration) |
| Truncated Normal, Gamma (experimental) | Numerical (2D Newton + quadrature) |
| Partial specification (`@dist Gamma(3, _) ...`) | Numerical root-finding |
| Quantile matching (Gamma, Beta) | Numerical root-finding |
| Quantile matching (Normal, Laplace, Logistic, Cauchy, Gumbel) | Direct formula (location-scale) |
| Arbitrary support (affine) | Moment transform + `LocationScale` wrapper |
