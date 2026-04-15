"""
    @dist expr

Create a distribution specification. Returns either a `Type`, a distribution instance,
or a `PartialDist` depending on the form:

- `@dist Gamma` → `Gamma` (the type itself)
- `@dist TDist(7)` → `TDist(7)` (a full instance)
- `@dist Gamma(3.0, _)` → `PartialDist{Gamma}((α=3.0, θ=missing))` (partial spec)

The result can be passed to any `dist_from_*` function:

```julia
d = dist_from_mean_var(@dist(Gamma(3.0, _)), 5.0, 3.0)
d = dist_from_mean(@dist(Gamma(3.0, _)), 5.0)
d = dist_from_var(@dist(Logistic(2.0, _)), 22.3)
d = dist_from_mean_var(@dist(TDist(7)), 5.0, 2.0)
d = dist_from_mean_var(@dist(Gamma), 5.0, 3.0)
```
"""
macro dist(expr)
    if expr isa Symbol
        # Bare type: @dist Gamma → Gamma
        return esc(expr)

    elseif expr isa Expr && expr.head == :call
        D = expr.args[1]
        args = expr.args[2:end]

        has_placeholder = any(a -> a == :_ || a == :(_), args)

        if has_placeholder
            # Partial: @dist Gamma(3.0, _) → PartialDist{Gamma}(...)
            param_exprs = [a == :_ || a == :(_) ? :missing : a for a in args]
            params_tuple = Expr(:tuple, param_exprs...)
            return esc(:(PartialDist{$D}(NamedTuple{fieldnames($D)}($params_tuple))))
        else
            # Full instance: @dist TDist(7) → TDist(7)
            return esc(expr)
        end
    else
        error("@dist: expected a distribution type or constructor call, got $expr")
    end
end
