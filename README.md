# DistributionsFactories.jl

A Julia package for creating probability distributions parameterised by moments (mean, variance, etc.) or quantiles, on standard or arbitrary supports. Builds on [Distributions.jl](https://github.com/JuliaStats/Distributions.jl).

The `..` operator from [IntervalSets.jl](https://github.com/JuliaMath/IntervalSets.jl) is re-exported, so no separate `using IntervalSets` is needed.

## Installation

Once registered in the Julia General registry:

```julia
using Pkg
Pkg.add("DistributionsFactories")
```

Or directly from the repository:

```julia
using Pkg
Pkg.add(url="https://github.com/Distribution-Matching/DistributionsFactories.jl.git")
```

## Quick start

```julia
using DistributionsFactories
using Distributions

# Construct a Gamma from mean and variance
d = make_dist(Gamma, mean=5.0, var=3.0)
mean(d), var(d)   # (5.0, 3.0)

# From mean and coefficient of variation
d = make_dist(Gamma, mean=10.0, cv=0.5)

# From a single quantile
d = make_dist(Exponential, median=2.0)

# From two quantiles
d = make_dist(Normal, q1=10.0, q3=30.0)

# Check feasibility
dist_exists(Beta, mean=0.5, var=0.1)   # true

# What distributions are feasible for these moments on [0,∞)?
available_distributions(0..Inf, mean=5.0, var=3.0)
# [Gamma, Erlang, LogNormal, Weibull, Frechet, ...]
```

## Extended examples

```julia
using DistributionsFactories, Distributions

# --- Moment-based construction ---

d = make_dist(Gamma, mean=5.0, var=3.0)
d = make_dist(Gamma, mean=5.0, std=1.7)
d = make_dist(Gamma, mean=5.0, cv=0.5)
d = make_dist(Gamma, mean=5.0, scv=0.25)
d = make_dist(Exponential, mean=3.0)           # 1 parameter: mean determines all
d = make_dist(Exponential, var=4.0)             # 1 parameter: var determines all

# --- Quantile-based construction ---

d = make_dist(Exponential, median=2.0)
d = make_dist(Normal, q1=10.0, q3=30.0)
d = make_dist(Normal, median=5.0, iqr=4.0)
d = make_dist(Beta, mean=0.4, median=0.35)
d = make_dist(Gamma, quantiles=[(0.1, 1.0), (0.9, 10.0)])

# --- Mode-based construction ---

d = make_dist(Rayleigh, mode=2.0)                  # 1 parameter: mode determines all
d = make_dist(Gamma, mean=5.0, mode=3.0)           # mean + mode
d = make_dist(Beta, mean=0.4, mode=0.3)            # mean + mode
d = make_dist(Normal, mode=3.0, var=4.0)           # mode + variance
d = make_dist(Gamma, mode=2.0, var=3.0)            # mode + variance (numerical)
d = make_dist(Gamma, mode=2.0, median=3.0)         # mode + median
d = make_dist(Gamma, mode=2.0, q3=5.0)             # mode + quantile
d = make_dist(Normal, mode=3.0, iqr=4.0)           # mode + IQR
d = make_dist(Gamma, mode=2.0, iqr=3.0)            # mode + IQR (numerical)

# --- Arbitrary supports ---

d = make_dist(Beta, mean=3.5, var=0.5, support=2..7)       # Beta on [2,7]
d = make_dist(Gamma, mean=8.0, var=3.0, support=3..Inf)    # Gamma on [3,∞)
d = make_dist(Gamma, mean=5.0, var=3.0, support=-Inf..10)  # flipped Gamma on (-∞,10]
d = make_dist(Binomial, mean=12.0, var=1.2, support=10:15) # Binomial on {10,...,15}

# --- Partial specification with @dist ---

spec = @dist Gamma(3.0, _)            # fix α=3, θ to be solved
d = make_dist(spec, mean=5.0)         # solve θ from mean

spec = @dist Logistic(2.0, _)         # fix μ=2, θ to be solved
d = make_dist(spec, var=22.3)         # solve θ from variance

spec = @dist TDist(7)                 # full instance, ν=7 degrees of freedom
d = make_dist(spec, mean=5.0, var=2.0)  # LocationScale wrap

# --- Discovery ---

available_distributions(0..Inf, mean=5.0, var=25.0)
available_distributions(support(Beta(2, 3)), mean=0.5, std=0.2)
available_distributions(mean=5.0, cv=1.0)
```

## API reference

### `make_dist(D; kwargs...)`

The primary interface. Constructs a distribution from any combination of moment
and quantile keywords. `D` can be a type, instance, or `DistSpec` (via `@dist`).

| Keyword | Meaning |
|---------|---------|
| `mean` | Target mean μ̄ |
| `var` | Target variance σ̄² |
| `std` | Target standard deviation σ̄ |
| `cv` | Target coefficient of variation σ̄/μ̄ |
| `scv` | Target squared coefficient of variation σ̄²/μ̄² |
| `second_moment` | Target second moment E[X²] |
| `mode` | Target mode (peak of PDF/PMF) |
| `median` | Target median (0.5 quantile) |
| `q1` | Target first quartile (0.25 quantile) |
| `q3` | Target third quartile (0.75 quantile) |
| `iqr` | Target interquartile range (with `median`) |
| `quantiles` | Target quantiles as vector of `(p, q)` tuples |
| `support` | Target support interval (`a..b`, `a..Inf`, `a:b`) |

### `dist_exists(D; kwargs...)`

Pure predicate — returns `true` or `false`; never throws for infeasible moments.
Use `make_dist` to actually construct (that one throws `DomainError` with a
reason explaining why).

```julia
dist_exists(Beta, mean=0.5, var=0.1)           # true
dist_exists(Exponential, mean=2.5, var=1.5)    # false (Exponential needs σ̄² = μ̄²)
```

### `available_distributions(support; kwargs...)`

Discover which distribution types are feasible for a given support and/or moments.

```julia
available_distributions(0..Inf)                    # all positive distributions
available_distributions(0..Inf, mean=5.0, var=3.0) # that can fit these moments
available_distributions(mean=5.0, cv=1.0)          # across all supports
```

### The `@dist` macro and `DistSpec`

The `@dist` macro creates distribution specifications. Use `_` as a placeholder for
free parameters (to be determined). The macro's only job is to parse the expression
— all solving is done by `make_dist` or the lower-level `dist_from_*` functions.

**Three forms:**

```julia
# 1. Bare type — returns the type itself
@dist Gamma           # → Gamma

# 2. Partial — underscore marks parameters to solve
@dist Gamma(3.0, _)   # → DistSpec{Gamma}  (α=3.0 fixed, θ=missing)
@dist Normal(_, 1.0)  # → DistSpec{Normal} (μ=missing, σ=1.0 fixed)

# 3. Full instance — returns the distribution as-is
@dist TDist(7)        # → TDist(7)
@dist Normal(3.0, 2.0) # → Normal(3.0, 2.0)
```

**Inspect a `DistSpec`:**

```julia
spec = @dist Gamma(3.0, _)
fixed_params(spec)   # (α = 3.0,)
free_params(spec)    # (:θ,)
typeof(spec)         # DistSpec{Gamma, @NamedTuple{α::Float64, θ::Missing}}
```

**Pass to `make_dist` to solve the free parameter(s):**

```julia
spec = @dist Gamma(3.0, _)
d = make_dist(spec, mean=5.0)          # solve θ from mean
mean(d), var(d)                         # (5.0, 8.33)
params(d)                               # (3.0, 1.667)

d = make_dist(spec, var=3.0)           # solve θ from variance
mean(d), var(d)                         # (3.0, 3.0)

spec = @dist Logistic(2.0, _)
d = make_dist(spec, std=4.0)           # fix location, solve scale from std
mean(d), std(d)                         # (2.0, 4.0)

spec = @dist Normal(_, 1.0)
d = make_dist(spec, mean=3.0)          # fix σ, solve μ from mean
mean(d), std(d)                         # (3.0, 1.0)

spec = @dist Beta(2.0, _)
d = make_dist(spec, mean=0.4)          # fix α, solve β from mean
mean(d), var(d)                         # (0.4, 0.04)
params(d)                               # (2.0, 3.0)
```

**Full instances wrap in LocationScale when moments are given:**

```julia
spec = @dist TDist(7)
d = make_dist(spec, mean=5.0, var=2.0) # wraps TDist(7) in LocationScale
mean(d), var(d)                         # (5.0, 2.0)
```

**Note on syntax:** When used standalone, no parentheses are needed
(`@dist Gamma(3.0, _)`). When used inline as a function argument,
parentheses are required so Julia parses it correctly:

```julia
# Standalone — no parens
spec = @dist Gamma(3.0, _)
d = make_dist(spec, mean=5.0)

# Inline — parens required
d = make_dist(@dist(Gamma(3.0, _)), mean=5.0)
```

### Internal functions

The `make_dist` API is built on these internal functions. They are not exported but
can be accessed via `DistributionsFactories.dist_from_mean_var(...)` etc.:

| Function | Used for |
|----------|----------|
| `dist_from_mean_var(D, μ̄, σ̄²)` | Per-distribution moment matching (25+ methods) |
| `dist_from_mean(D, μ̄)` | 1-parameter distributions (Exponential, Poisson, etc.) |
| `dist_from_var(D, σ̄²)` | 1-parameter distributions where mean is determined by var |
| `dist_from_quantile(D, p, q)` | Single quantile matching |
| `dist_from_quantiles(D, p1, q1, p2, q2)` | Two-quantile matching |
| `dist_from_mean_quantile(D, μ̄, p, q)` | Mean + quantile matching |
| `exists_dist_from_mean_var(D, μ̄, σ̄²)` | Feasibility checking |

All conversions (std→var, cv→var, median→quantile, etc.) are handled by `make_dist`
internally. These functions also accept `DistSpec` as the first argument.

## Distributions covered

The following table lists all distributions covered by the package. **Yes** = fully implemented, **No** = not yet implemented.

### Continuous distributions

| # | Distribution | Support | Free params | Implemented | Method | Notes |
|---|---|---|---|---|---|---|
| 1 | Student's T | (-∞, ∞) | 3 | Yes | Direct formula | Type `TDist`; instance form for location-scale |
| 2 | Normal | (-∞, ∞) | 2 | Yes | Direct formula | |
| 3 | Cauchy | (-∞, ∞) | 2 | Quantile only | — | Moments undefined; use quantile matching |
| 4 | Laplace | (-∞, ∞) | 2 | Yes | Direct formula | |
| 5 | Logistic | (-∞, ∞) | 2 | Yes | Direct formula | |
| 6 | Gumbel | (-∞, ∞) | 2 | Yes | Direct formula | |
| 7 | Gamma | [0, ∞) | 2 | Yes | Direct formula | |
| 8 | Erlang | [0, ∞) | 2 | Yes | Direct formula | Special case of Gamma (integer shape) |
| 9 | Exponential | [0, ∞) | 1 | Yes | Direct formula | Special case of Gamma; 1 free parameter |
| 10 | Chi-squared | [0, ∞) | 1 | Yes | Direct formula | Special case of Gamma; 1 free parameter |
| 11 | Log-normal | [0, ∞) | 2 | Yes | Direct formula | |
| 12 | Weibull | [0, ∞) | 2 | Yes | Numerical (root-finding) | Beta-ratio equation |
| 13 | Frechet | [0, ∞) | 2 | Yes | Numerical (root-finding) | Beta-ratio equation |
| 14 | F distribution | [0, ∞) | 2 | Yes | Direct formula | Requires 1 < μ̄ < 2 |
| 15 | Inverse Gamma | [0, ∞) | 2 | Yes | Direct formula | |
| 16 | Folded Normal | [0, ∞) | 1 | Yes | Numerical (2D Newton) | Returns parent Normal |
| 17 | Pareto | [0, ∞) | 1 | Yes | Direct formula | |
| 18 | Chi | [0, ∞) | 1 | Yes | Direct formula | |
| 19 | Rayleigh | [0, ∞) | 1 | Yes | Direct formula | Special case of Chi; 1 free parameter |
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
| 33 | Geometric | {0, 1, 2, …} | 1 | Yes | Direct formula | Special case of Neg. Binomial; 1 free parameter |
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
| Numerical root-finding | `DistSpec` — solving free parameters from moments |
| Location-scale formula | Quantile matching for Normal, Laplace, Logistic, Cauchy, Gumbel |
| Numerical root-finding | Quantile matching for Gamma, Beta |
| Moment transform + `LocationScale` | Arbitrary support (affine) |
