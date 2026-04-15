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
| 16 | Folded Normal | [0, Inf) | 1 | Yes | Numerical (2D Newton) | Returns parent Normal |
| 17 | Pareto | [0, Inf) | 1 | Yes | Direct formula | |
| 18 | Chi | [0, Inf) | 1 | Yes | Direct formula | |
| 19 | Rayleigh | [0, Inf) | 1 | Yes | Direct formula | Special case of Chi; 1 free parameter |
| 20 | Half-truncated Normal | [0, Inf) | 2 | No | | |
| 21 | Half-truncated Laplace | [0, Inf) | 2 | No | | |
| 22 | Half-truncated Logistic | [0, Inf) | 2 | No | | |
| 23 | Beta | [0, 1] | 2 | Yes | Direct formula | |
| 24 | Uniform | [0, 1] | 0 | Yes | Direct formula | Adjusts support |
| 25 | Symmetric Triangular | [0, 1] | 0 | Yes | Direct formula | Type `SymTriangularDist`; adjusts support |
| 26 | Triangular | [0, 1] | 1 | No | | |
| 27 | Doubly-truncated Normal | [0, 1] | 2 | No | | |
| 28 | Doubly-truncated Laplace | [0, 1] | 2 | No | | |
| 29 | Doubly-truncated Logistic | [0, 1] | 2 | No | | |

## Discrete distributions

| # | Distribution | Support | Free params | Implemented | Method | Notes |
|---|---|---|---|---|---|---|
| 30 | Binomial | {0, ..., n} | 1 | Yes | Direct formula | |
| 31 | Discrete Uniform | {0, ..., n} | 0 | Yes | Direct formula | Adjusts support |
| 32 | Negative Binomial | {0, 1, 2, ...} | 2 | Yes | Direct formula | |
| 33 | Geometric | {0, 1, 2, ...} | 1 | Yes | Direct formula | Special case of Neg. Binomial; 1 free parameter |
| 34 | Poisson | {0, 1, 2, ...} | 1 | Yes | Direct formula | 1 free parameter |
| 35 | Discrete Triangular | {0, ..., n} | 1 | No | | |
| 36 | Discrete Sym. Triangular | {0, ..., n} | 0 | No | | |
| 37 | Truncated Poisson | {0, ..., n} | 1 | No | | |

## Experimental / work in progress

The following have initial code but are not fully validated:

- **Truncated Normal on [a,b]** -- 2D Newton solver with quadrature. Passes basic tests.
- **Truncated Gamma on [a,b]** -- same solver. Passes basic tests.
- **Truncated Laplace on [a,b]** -- same solver. Not yet tested end-to-end.
- **Truncated Logistic on [a,b]** -- same solver. Not yet tested end-to-end.
