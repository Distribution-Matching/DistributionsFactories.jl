# API Reference

## Public API

```@docs
make_dist
dist_exists
available_distributions
@dist
DistSpec
fixed_params
free_params
```

## Internal functions

These are not exported but can be accessed via `DistributionsFactories.func_name(...)`.

```@docs
DistributionsFactories.dist_from_mean_var
DistributionsFactories.dist_from_mean
DistributionsFactories.dist_from_var
DistributionsFactories.dist_from_quantile
DistributionsFactories.dist_from_quantiles
DistributionsFactories.dist_from_mean_quantile
DistributionsFactories.exists_dist_from_mean_var
DistributionsFactories.dist_from_mean_var_on_support
```
