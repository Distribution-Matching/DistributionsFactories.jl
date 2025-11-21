# DistributionsFactories.jl
A Julia package for creating probability distributions parametarised by moments (mean, variance, etc.), or other measures  (median, constrained domains, etc.). 

```
using DistributionsFactories
using Distributions
```

The general interface uses functions like `dist_from_mean_var`, `dist_from_mean`, `dist_from_quantiles` etc. For example: 

```
julia> d = dist_from_mean_var(Gamma, 2.5, 1.5)
Gamma{Float64}(α=4.166666666666667, θ=0.6)

julia> mean(d), var(d)
(2.5, 1.5)
```

The functions work when the distributional family can be parameterized as such, and otherwise an exception is thrown. For example:

```
julia> d = dist_from_mean_var(Exponential, 2.5, 1.5)
ERROR: The Exponential distribution has a single parameter
```

The complete list of functions is:

* `dist_from_mean`
* `dist_from_var`
* `dist_from_std`
* `dist_from_cv`
* `dist_from_scv`
* `dist_from_mean_var`, `dist_from_mean_std`, `dist_from_mean_cv`, `dist_from_mean_scv`
* `dist_from_quantile`
* `dist_from_qunatiles`
* `dist_from_median`
* `dist_from_q1_q3`
* `dist_from_iqr`
* `dist_from_median_iqr`
* `dist_from_mean_quantile`
* `dist_from_mean_median`

The first argument to each of these functions is a distribution type such as `Gamma` in the example above. Another option for the first argument is a 2-tuple with the first element being a distribution type, and the second element a named tuple of partial parameters. For example,

```
julia> d = dist_from_mean((Gamma,(α=4.166666666666667;),  1.5)
Gamma{Float64}(α=4.166666666666667, θ=0.6)
```


The supported (non-truncated) continuous distributions are: `Beta`, `Cauchy`, `Chi`, `Chisq`, `Erlang`, `Exponential`, `FDist`, `Frechet`, `Gamma`, `Gumbel`, `InverseGamma`, `Laplace`, `Logistic`, `LogNormal`, `Normal`, `Rayleigh`, `TDist`, `TriangularDist`, `Uniform`, `Weibull`.

The supported discrete distributions are: `Bernoulli`, `Binomial`, `Geometric`, `NegativeBinomial`, `Poisson`.

The supported truncated distributions are: `Normal` and `Laplace`.

In certain cases 

