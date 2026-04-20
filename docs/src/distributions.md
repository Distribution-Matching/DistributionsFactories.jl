# Distributions Covered

## Continuous distributions

| # | Distribution | Support | Free params | Implemented | Method | Notes |
|---|---|---|---|---|---|---|
| 1 | Student's T | (-Inf, Inf) | 3 | Yes | Direct formula | Type `TDist`; instance form for location-scale |
| 2 | Normal | (-Inf, Inf) | 2 | Yes | Direct formula | |
| 3 | Cauchy | (-Inf, Inf) | 2 | Quantile only | -- | Moments undefined; use quantile matching |
| 4 | Laplace | (-Inf, Inf) | 2 | Yes | Direct formula | |
| 5 | Logistic | (-Inf, Inf) | 2 | Yes | Direct formula | |
| 6 | Gumbel | (-Inf, Inf) | 2 | Yes | Direct formula | |
| 7 | Gamma | [0, Inf) | 2 | Yes | Direct formula | |
| 8 | Erlang | [0, Inf) | 2 | Yes | Direct formula | Special case of Gamma (integer shape) |
| 9 | Exponential | [0, Inf) | 1 | Yes | Direct formula | Special case of Gamma; 1 free parameter |
| 10 | Chi-squared | [0, Inf) | 1 | Yes | Direct formula | Special case of Gamma; 1 free parameter |
| 11 | Log-normal | [0, Inf) | 2 | Yes | Direct formula | |
| 12 | Weibull | [0, Inf) | 2 | Yes | Numerical (root-finding) | Beta-ratio equation |
| 13 | Frechet | [0, Inf) | 2 | Yes | Numerical (root-finding) | Beta-ratio equation |
| 14 | F distribution | [0, Inf) | 2 | Yes | Direct formula | Requires 1 < mean < 2 |
| 15 | Inverse Gamma | [0, Inf) | 2 | Yes | Direct formula | |
| 16 | Folded Normal | [0, Inf) | 2 | Yes | Numerical (2D Newton) | Custom `FoldedNormal(μ, σ)` type — see [Extension distributions](#extension-distributions) |
| 17 | Pareto | [0, Inf) | 1 | Yes | Direct formula | |
| 18 | Chi | [0, Inf) | 1 | Yes | Direct formula | |
| 19 | Rayleigh | [0, Inf) | 1 | Yes | Direct formula | Special case of Chi; 1 free parameter |
| 20 | Half-truncated Normal | [a, ∞) | 2 | Yes | Numerical (standardize+2D Newton) | Use `truncated(Normal(); lower=a)` template; feasibility σ̄² < (μ̄-a)² |
| 21 | Half-truncated Laplace | [a, ∞) | 2 | Yes | Numerical (standardize+2D Newton); boundary handled directly | Boundary σ̄² = (μ̄-a)² attained as Exponential |
| 22 | Half-truncated Logistic | [a, ∞) | 2 | Yes | Numerical (standardize+2D Newton) | Same exp-tail bound as Normal |
| 23 | Beta | [0, 1] | 2 | Yes | Direct formula | |
| 24 | Uniform | [0, 1] | 0 | Yes | Direct formula | Adjusts support |
| 25 | Symmetric Triangular | [0, 1] | 0 | Yes | Direct formula | Type `SymTriangularDist`; adjusts support |
| 26 | Triangular | ℝ | 3 | Yes (mean+var+mode) | Direct formula | `TriangularDist(a, b, c)` — supply mode in addition to mean/var |
| 27 | Doubly-truncated Normal | [a, b] | 2 | Yes | Numerical (standardize to [-0.5, 0.5] + 2D Newton) | Langevin-envelope feasibility |
| 28 | Doubly-truncated Laplace | [a, b] | 2 | Yes | Numerical (standardize to [-0.5, 0.5] + 2D Newton) | Same envelope as Normal |
| 29 | Doubly-truncated Logistic | [a, b] | 2 | Yes | Numerical (standardize to [-0.5, 0.5] + 2D Newton) | Same envelope as Normal |

## Discrete distributions

| # | Distribution | Support | Free params | Implemented | Method | Notes |
|---|---|---|---|---|---|---|
| 30 | Binomial | {0, ..., n} | 1 | Yes | Direct formula | |
| 31 | Discrete Uniform | {0, ..., n} | 0 | Yes | Direct formula | Adjusts support |
| 32 | Negative Binomial | {0, 1, 2, ...} | 2 | Yes | Direct formula | |
| 33 | Geometric | {0, 1, 2, ...} | 1 | Yes | Direct formula | Special case of Neg. Binomial; 1 free parameter |
| 34 | Poisson | {0, 1, 2, ...} | 1 | Yes | Direct formula | 1 free parameter |
| 35 | Discrete Triangular | {a, ..., b} | 3 | Yes (mean+var+mode) | Approx. (continuous solve + integer search) | Custom `DiscreteTriangular(a, b, c)` type |
| 36 | Discrete Sym. Triangular | {μ-n, ..., μ+n} | 2 | Yes | Direct formula | Custom `DiscreteSymmetricTriangular(μ, n)` type |
| 37 | Truncated Poisson | {a, ..., b} | 1 | Yes | Numerical (1D root-find) | Use `truncated(Poisson(0); lower=a, upper=b)` as template |

## Extension distributions

The following types are defined in `src/extensions/` because they are *not*
provided by `Distributions.jl`. They implement the full Distributions.jl
interface (`pdf`, `logpdf`, `cdf`, `quantile`, `mean`, `var`, `rand`,
`support`, `params`) and behave as first-class distributions.

### `FoldedNormal(μ, σ)`

The distribution of `|X|` where `X ~ Normal(μ, σ)`. Support `[0, ∞)`.
Setting `μ = 0` recovers the half-normal.

```julia
d = FoldedNormal(2.0, 1.0)
mean(d), var(d)             # 2.0170, 0.9318
make_dist(FoldedNormal, mean=2.5, var=1.2)
```

### `DiscreteSymmetricTriangular(μ, n)`

Integer-valued symmetric triangular on `{μ-n, …, μ+n}` with PMF
`P(μ + k) = (n + 1 − |k|) / (n + 1)²`. Closed-form `mean = μ`,
`var = n(n+2)/6`.

```julia
make_dist(DiscreteSymmetricTriangular, mean=4, var=2.5)
```

### `DiscreteTriangular(a, b, c)`

Integer-valued asymmetric triangular on `{a, …, b}` with mode `c`. PMF is
two linear ramps meeting at `c`. The factory needs `mode` in addition to
`mean` and `var`:

```julia
make_dist(DiscreteTriangular, mean=5.0, var=2.0, mode=5)
```

## Truncated Poisson

`Distributions.jl` already provides `truncated(Poisson(λ); lower=a, upper=b)`;
the factory follows the same template-instance pattern as the truncated
Normal/Laplace/Logistic constructors:

```julia
template = truncated(Poisson(0.0); lower=0, upper=10)
make_dist(template, mean=3.0, var=v)            # var must be consistent
DistributionsFactories.dist_from_mean(template, 3.0)  # mean-only path
```

### Half-truncated Student-$t$

`Truncated{<:TDist}` on `[a, ∞)` or `(-∞, b]` has a feasibility predicate
based on the Pareto-tail bound `σ̄² < ν/(ν-2)·gap²` (requires `ν > 2`). The
constructor is **not yet implemented** — only feasibility checking. See
the GitHub issue tracker for the open construction work.

## Experimental / work in progress

- **Doubly-truncated Gamma on [a, b]** — generic 2D Newton solver in user
  coordinates (the standardize-and-solve refactor only applies to
  location-scale families). One happy-path test, no roundtrip sweep.
- **Doubly-truncated Student-$t$** — only feasibility (loose Bhatia–Davis
  fallback). Tight Pareto-tail dome and constructor are open.
- **Solver convergence near the canonical boundary** for
  `Truncated{<:Normal/Laplace/Logistic}` — fails when `|μ̄_std|` is close
  to `±0.5` and variance is small. Documented in GitHub issue #1.
