# DistributionsFactories.jl

A Julia package for creating probability distributions parameterised by moments (mean, variance, etc.), quantiles, or mode, on standard or arbitrary supports.

## Installation

```julia
using Pkg
Pkg.add("DistributionsFactories")  # once registered, or:
# Pkg.add(url="https://github.com/Distribution-Matching/DistributionsFactories.jl.git")
```

## Quick start

```julia
using DistributionsFactories
using Distributions

# Construct a Gamma from mean and variance
d = make_dist(Gamma, mean=5.0, var=3.0)

# From mode and median
d = make_dist(Gamma, mode=2.0, median=3.0)

# Fix some parameters with @dist, solve the rest
spec = @dist Gamma(3.0, _)
d = make_dist(spec, mean=5.0)

# Discover feasible distributions
available_distributions(0..Inf, mean=5.0, var=3.0)

# Check feasibility
dist_exists(Beta, mean=0.5, var=0.1)
```

## Key features

- **Unified API**: `make_dist(D; mean=..., var=..., mode=..., median=..., support=...)`
- **30+ distributions** with moment-matching (direct formula or numerical)
- **Arbitrary supports** via affine transform or truncation: `make_dist(Beta, mean=3.5, var=0.5, support=2..7)`
- **Partial specification** with `@dist` macro: `@dist Gamma(3.0, _)` fixes shape, solves scale
- **Mode-based construction**: `make_dist(Gamma, mode=2.0, iqr=3.0)`
- **Quantile matching**: `make_dist(Normal, q1=10.0, q3=30.0)`
- **Discovery**: `available_distributions(0..Inf, mean=5.0, var=3.0)`

## Keyword arguments for `make_dist`

| Keyword | Meaning |
|---------|---------|
| `mean` | Target mean |
| `var` | Target variance |
| `std` | Target standard deviation |
| `cv` | Target coefficient of variation |
| `scv` | Target squared coefficient of variation |
| `second_moment` | Target second moment E[X^2] |
| `mode` | Target mode (peak of PDF/PMF) |
| `median` | Target median |
| `q1` | Target first quartile |
| `q3` | Target third quartile |
| `iqr` | Target interquartile range (with `median` or `mode`) |
| `quantiles` | Target quantiles as vector of `(p, q)` tuples |
| `support` | Target support interval (`a..b`, `a..Inf`, `a:b`) |

See [API Reference](@ref) for full details.
