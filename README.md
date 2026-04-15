# DistributionsFactories.jl

A Julia package for creating probability distributions parameterised by moments (mean, variance, etc.) or quantiles, on standard or arbitrary supports. Builds on [Distributions.jl](https://github.com/JuliaStats/Distributions.jl).

The `..` operator from [IntervalSets.jl](https://github.com/JuliaMath/IntervalSets.jl) is re-exported, so no separate `using IntervalSets` is needed.

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

# Fix some parameters, solve the rest with @dist
spec = @dist Gamma(3.0, _)
d = dist_from_mean(spec, 5.0)   # fix α=3, solve θ from mean

# What distributions are feasible for these moments on [0,∞)?
available_distributions(0..Inf, mean=2.5, var=1.5)
# [Gamma, Erlang, LogNormal, Weibull, Frechet, ...]

# Check feasibility
exists_dist_from_mean_var(Beta, 0.5, 0.1)   # true

# Infeasible moments throw an error
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

# Exponential from variance alone (1 DOF: mean is determined)
d = dist_from_var(Exponential, 4.0)
mean(d)   # 2.0

# --- Arbitrary supports ---

# Beta scaled from [0,1] to [2, 7]
d = dist_from_mean_var_on_support(Beta, 3.5, 0.5, support=2..7)
mean(d), minimum(d), maximum(d)   # (3.5, 2.0, 7.0)

# Gamma shifted to [3, ∞)
d = dist_from_mean_var_on_support(Gamma, 8.0, 3.0, support=3..Inf)
mean(d), minimum(d)   # (8.0, 3.0)

# Flipped Gamma on (-∞, 10]
d = dist_from_mean_var_on_support(Gamma, 5.0, 3.0, support=-Inf..10)
mean(d), maximum(d)   # (5.0, 10.0)

# Binomial shifted from {0,...,5} to {10,...,15}
d = dist_from_mean_var_on_support(Binomial, 12.0, 1.2, support=10:15)
minimum(d), maximum(d)   # (10, 15)

# --- Partial specification with @dist ---

# Fix shape, solve scale from mean
spec = @dist Gamma(3.0, _)
d = dist_from_mean(spec, 5.0)
params(d)   # (3.0, 1.667)

# Fix location, solve scale from variance
spec = @dist Logistic(2.0, _)
d = dist_from_var(spec, 22.3)
params(d)   # (2.0, 2.604)

# TDist with fixed degrees of freedom, arbitrary mean and variance
spec = @dist TDist(7)
d = dist_from_mean_var(spec, 5.0, 2.0)
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

# Use the support of an existing distribution
available_distributions(support(Beta(2, 3)), mean=0.5, std=0.2)

# Search across all supports
available_distributions(mean=5.0, cv=1.0)
```

## The `@dist` macro and `PartialDist`

The `@dist` macro creates distribution specifications where `_` marks parameters
to be solved from moment constraints. The result is a `PartialDist` that can be
passed to any `dist_from_*` function.

```julia
# Create a partial spec — α is fixed, θ is to be solved
spec = @dist Gamma(3.0, _)
fixed_params(spec)   # (α = 3.0,)
free_params(spec)    # (:θ,)
```

Use it with any `dist_from_*` function — the free parameter is solved numerically:

```julia
spec = @dist Gamma(3.0, _)              # fix α=3, θ to be solved
d = dist_from_mean(spec, 5.0)           # solve θ from mean

spec = @dist Logistic(2.0, _)           # fix μ=2, θ to be solved
d = dist_from_var(spec, 22.3)           # solve θ from variance
d = dist_from_std(spec, 4.0)            # solve θ from std

spec = @dist Normal(0.0, _)             # fix μ=0, σ to be solved
d = dist_from_mean_var(spec, 0.0, 4.0)  # solve σ from mean+var

spec = @dist Gamma(4.0, _)              # fix α=4, θ to be solved
d = dist_from_mean_cv(spec, 5.0, 0.5)   # solve θ from mean+cv

spec = @dist Normal(_, 2.0)             # fix σ=2, μ to be solved
d = dist_from_mean_std(spec, 3.0, 2.0)  # solve μ from mean+std

spec = @dist Beta(2.0, _)               # fix α=2, β to be solved
d = dist_from_mean(spec, 0.4)           # solve β from mean
```

The macro also handles bare types and full instances:

```julia
spec = @dist Gamma                       # bare type
d = dist_from_mean_var(spec, 5.0, 3.0)   # same as dist_from_mean_var(Gamma, ...)

spec = @dist TDist(7)                    # full instance
d = dist_from_mean_var(spec, 5.0, 2.0)   # wraps in LocationScale
```

The `PartialDist` type encodes which parameters are fixed vs free in the type system:

```julia
typeof(@dist Gamma(3.0, _))
# PartialDist{Gamma, @NamedTuple{α::Float64, θ::Missing}}
```

## Moment-based construction

### `dist_from_mean_var(D, μ̄, σ̄²)`

The primary interface. Constructs a distribution of type `D` with mean `μ̄` and variance `σ̄²`. The first argument `D` can be a type, a distribution instance, or a `PartialDist` (via `@dist`).

```julia
# Type: standard TDist with μ=0, determines ν from variance
d = dist_from_mean_var(TDist, 0.0, 3.0)

# Instance: TDist with fixed ν=7, arbitrary mean and variance via location-scale
d = dist_from_mean_var(TDist(7), 5.0, 2.0)

# PartialDist: fix one parameter, solve the other
spec = @dist Normal(0.0, _)
d = dist_from_mean_var(spec, 0.0, 4.0)
```

### Convenience wrappers

All of these accept a type, instance, or `PartialDist` as the first argument.

| Function | Specification |
|----------|--------------|
| `dist_from_mean_var(D, μ̄, σ̄²)` | Mean and variance |
| `dist_from_mean_std(D, μ̄, σ̄)` | Mean and standard deviation |
| `dist_from_mean_cv(D, μ̄, cv)` | Mean and coefficient of variation |
| `dist_from_mean_scv(D, μ̄, scv)` | Mean and squared coefficient of variation |
| `dist_from_mean_second_moment(D, μ̄, m2)` | Mean and second moment E[X²] |
| `dist_from_mean(D, μ̄)` | Mean only |
| `dist_from_var(D, σ̄²)` | Variance only |
| `dist_from_std(D, σ̄)` | Standard deviation only |

### Feasibility checking

Each `dist_from_*` function has a corresponding `exists_dist_from_*` that returns `true` or throws `DomainError`:

| Function | Checks feasibility for |
|----------|----------------------|
| `exists_dist_from_mean_var(D, μ̄, σ̄²)` | `dist_from_mean_var` |
| `exists_dist_from_mean_std(D, μ̄, σ̄)` | `dist_from_mean_std` |
| `exists_dist_from_mean_cv(D, μ̄, cv)` | `dist_from_mean_cv` |
| `exists_dist_from_mean_scv(D, μ̄, scv)` | `dist_from_mean_scv` |
| `exists_dist_from_mean_second_moment(D, μ̄, m2)` | `dist_from_mean_second_moment` |

## Distributions on arbitrary supports

Use `dist_from_mean_var_on_support` to place a distribution on a non-standard domain.

The function automatically determines whether to use an **affine transform** or **truncation**:

- **Affine**: when the requested support has the same shape as the natural one (e.g. Beta on [2,7], Gamma on [3,∞)). Returns a `LocationScale` wrapper.
- **Truncation** (experimental): when the requested support restricts the natural one (e.g. Normal on [0,1], Gamma on [0,10]). Returns a `Truncated` wrapper. Uses a 2D Newton solver with quadrature — works but not yet fully validated for all cases.

```julia
dist_from_mean_var_on_support(Beta, 3.5, 0.5, support=2..7)          # affine
dist_from_mean_var_on_support(Gamma, 8.0, 3.0, support=3..Inf)       # affine shift
dist_from_mean_var_on_support(Gamma, 5.0, 3.0, support=-Inf..10)     # affine flip
dist_from_mean_var_on_support(Normal, 0.5, 0.04, support=0..1)       # truncation
dist_from_mean_var_on_support(Binomial, 12.0, 1.2, support=10:15)    # discrete shift
```

Supports are specified using IntervalSets syntax (`a..b`), Distributions.jl's `RealInterval` (e.g. `support(some_dist)`), or integer ranges (`a:b`).

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
available_distributions(0:10, mean=5.0, var=2.0)         # discrete
available_distributions(mean=5.0, cv=1.0)                # across all supports
```

Moment keywords: `mean`, `var`, `std`, `cv`, `scv`, `second_moment`.

## Distributions covered

The following table lists all distributions covered by the package. **Yes** = fully implemented, **No** = not yet implemented.

### Continuous distributions

| # | Distribution | Support | Free params | Implemented | Method | Notes |
|---|---|---|---|---|---|---|
| 1 | Student's T | (-∞, ∞) | 3 | Yes | Direct formula | Type `TDist`; instance form for location-scale |
| 2 | Normal | (-∞, ∞) | 2 | Yes | Direct formula | |
| 3 | Cauchy | (-∞, ∞) | 2 | Quantile only | — | Moments undefined; use `dist_from_quantiles` |
| 4 | Laplace | (-∞, ∞) | 2 | Yes | Direct formula | |
| 5 | Logistic | (-∞, ∞) | 2 | Yes | Direct formula | |
| 6 | Gumbel | (-∞, ∞) | 2 | Yes | Direct formula | |
| 7 | Gamma | [0, ∞) | 2 | Yes | Direct formula | |
| 8 | Erlang | [0, ∞) | 2 | Yes | Direct formula | Special case of Gamma (integer shape) |
| 9 | Exponential | [0, ∞) | 1 | Yes | Direct formula | Special case of Gamma; 1 DOF |
| 10 | Chi-squared | [0, ∞) | 1 | Yes | Direct formula | Special case of Gamma; 1 DOF |
| 11 | Log-normal | [0, ∞) | 2 | Yes | Direct formula | |
| 12 | Weibull | [0, ∞) | 2 | Yes | Numerical (root-finding) | Beta-ratio equation |
| 13 | Frechet | [0, ∞) | 2 | Yes | Numerical (root-finding) | Beta-ratio equation |
| 14 | F distribution | [0, ∞) | 2 | Yes | Direct formula | Requires 1 < μ̄ < 2 |
| 15 | Inverse Gamma | [0, ∞) | 2 | Yes | Direct formula | |
| 16 | Folded Normal | [0, ∞) | 1 | Yes | Numerical (2D Newton) | Returns parent Normal |
| 17 | Pareto | [0, ∞) | 1 | Yes | Direct formula | |
| 18 | Chi | [0, ∞) | 1 | Yes | Direct formula | |
| 19 | Rayleigh | [0, ∞) | 1 | Yes | Direct formula | Special case of Chi; 1 DOF |
| 20 | Half-truncated Normal | [0, ∞) | 2 | No | | |
| 21 | Half-truncated Laplace | [0, ∞) | 2 | No | | |
| 22 | Half-truncated Logistic | [0, ∞) | 2 | No | | |
| 23 | Beta | [0, 1] | 2 | Yes | Direct formula | |
| 24 | Uniform | [0, 1] | 0 | Yes | Direct formula | Adjusts support |
| 25 | Symmetric Triangular | [0, 1] | 0 | Yes | Direct formula | Type `SymTriangularDist`; adjusts support |
| 26 | Triangular | [0, 1] | 1 | No | | |
| 27 | Doubly-truncated Normal | [0, 1] | 2 | No | | |
| 28 | Doubly-truncated Laplace | [0, 1] | 2 | No | | |
| 29 | Doubly-truncated Logistic | [0, 1] | 2 | No | | |

### Discrete distributions

| # | Distribution | Support | Free params | Implemented | Method | Notes |
|---|---|---|---|---|---|---|
| 30 | Binomial | {0, …, n} | 1 | Yes | Direct formula | |
| 31 | Discrete Uniform | {0, …, n} | 0 | Yes | Direct formula | Adjusts support |
| 32 | Negative Binomial | {0, 1, 2, …} | 2 | Yes | Direct formula | |
| 33 | Geometric | {0, 1, 2, …} | 1 | Yes | Direct formula | Special case of Neg. Binomial; 1 DOF |
| 34 | Poisson | {0, 1, 2, …} | 1 | Yes | Direct formula | 1 DOF |
| 35 | Discrete Triangular | {0, …, n} | 1 | No | | |
| 36 | Discrete Sym. Triangular | {0, …, n} | 0 | No | | |
| 37 | Truncated Poisson | {0, …, n} | 1 | No | | |

### Experimental / work in progress

The following have initial code but are not fully validated:

- **Truncated Normal on [a,b]** — 2D Newton solver with quadrature. Passes basic tests.
- **Truncated Gamma on [a,b]** — same solver. Passes basic tests.
- **Truncated Laplace on [a,b]** — same solver. Not yet tested end-to-end.
- **Truncated Logistic on [a,b]** — same solver. Not yet tested end-to-end.

## Algorithms

| Method | Used by |
|---|---|
| Direct formula | Normal, Beta, Gamma, Logistic, Laplace, Gumbel, LogNormal, Pareto, InverseGamma, FDist, Chi, Rayleigh, Chisq, Exponential, Erlang, SymTriangularDist, Uniform, Binomial, DiscreteUniform, Poisson, NegativeBinomial, Geometric |
| Numerical root-finding | Weibull, Frechet (beta-ratio equation) |
| 2D Newton iteration | Folded Normal, Truncated distributions (experimental) |
| Numerical root-finding | `PartialDist` — solving free parameters from moments |
| Location-scale formula | Quantile matching for Normal, Laplace, Logistic, Cauchy, Gumbel |
| Numerical root-finding | Quantile matching for Gamma, Beta |
| Moment transform + `LocationScale` | Arbitrary support (affine) |
