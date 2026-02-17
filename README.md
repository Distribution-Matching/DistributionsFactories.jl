# DistributionsFactories.jl

A Julia package for creating probability distributions parameterised by moments (mean, variance, etc.).

## Installation

```julia
using Pkg
Pkg.add("DistributionsFactories")
```

## Usage

```julia
using DistributionsFactories
using Distributions
```

The main interface is `dist_from_mean_var`, which constructs a distribution from a given mean and variance:

```julia
julia> d = dist_from_mean_var(Gamma, 2.5, 1.5)
Gamma{Float64}(α=4.166666666666667, θ=0.6)

julia> mean(d), var(d)
(2.5, 1.5)
```

When no valid distribution exists for the given moments, an exception is thrown:

```julia
julia> dist_from_mean_var(Exponential, 2.5, 1.5)
ERROR: DomainError with "Exponential: the condition var = μ² is not satisfied"
```

You can check ahead of time whether a unique distribution exists using `exists_unique_dist_from_mean_var`:

```julia
julia> exists_unique_dist_from_mean_var(Beta, 0.5, 0.1)
true
```

## Supported distributions

### Continuous

`Beta`, `Chi`, `Chisq`, `Erlang`, `Exponential`, `FDist`, `Frechet`, `Gamma`, `Gumbel`, `InverseGamma`, `Laplace`, `Logistic`, `LogNormal`, `Normal`, `Rayleigh`, `TDist`, `Uniform`, `Weibull`

### Discrete

`Binomial`, `NegativeBinomial`, `Poisson`

### Weibull method selection

The Weibull `dist_from_mean_var` supports a `method` keyword argument:

```julia
dist_from_mean_var(Weibull, 5.0, 4.0; method = :simple_bracketing)  # default
dist_from_mean_var(Weibull, 5.0, 4.0; method = :oscar_garcia)
```

- `:simple_bracketing` (default) -- solves via `solve_beta_ratio`
- `:oscar_garcia` -- Oscar Garcia polynomial approximation as initial guess, refined with root finding

## Algorithms

In many cases the moment-to-parameter mapping is a direct formula. For certain distributions, specialised algorithms are used:

- **Weibull** -- numerical solution of the beta-ratio equation, or Oscar Garcia polynomial approximation with root refinement
- **Frechet** -- Victor Nawa & Saralees Nadarajah approximation with root refinement

## Coming soon

### Distributions

- `Bernoulli`, `Geometric`
- `TriangularDist`
- Truncated `Normal`, truncated `Laplace`

### Functions

- `dist_from_mean`, `dist_from_mean_std`, `dist_from_mean_cv`, `dist_from_mean_scv`
- `dist_from_var`, `dist_from_std`, `dist_from_cv`, `dist_from_scv`
- `dist_from_quantile`, `dist_from_quantiles`, `dist_from_median`, `dist_from_q1_q3`, `dist_from_iqr`
- `dist_from_median_iqr`
- `dist_from_mean_quantile`, `dist_from_mean_median`
- Partial parameter syntax, e.g. `dist_from_mean((Gamma, (α=4.167,)), 1.5)`
