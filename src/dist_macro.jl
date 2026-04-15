"""
    @dist expr kwargs...

Construct a distribution from a concise specification. Supports three forms:

1. **Type only** — all parameters solved from moments:
   ```julia
   @dist Gamma mean=5.0 var=3.0
   ```

2. **Partial specification** — `_` marks parameters to be solved:
   ```julia
   @dist Gamma(3.0, _) mean=5.0
   @dist Beta(2.0, _) mean=0.4
   ```

3. **Full instance** — all parameters given, LocationScale wrap if moments specified:
   ```julia
   @dist TDist(7) mean=5.0 var=2.0
   ```

Additional keyword arguments: `support` for arbitrary domain placement.

# Examples
```julia
@dist Gamma mean=5.0 var=3.0                    # standard
@dist Gamma(3.0, _) mean=5.0                    # fix α=3, solve θ
@dist TDist(7) mean=5.0 var=2.0                 # instance + LocationScale
@dist Beta mean=3.5 var=0.5 support=2..7        # with support
```
"""
macro dist(expr, kwargs...)
    # Parse keyword arguments
    kw_dict = Dict{Symbol, Any}()
    for kw in kwargs
        if kw isa Expr && kw.head == :(=)
            kw_dict[kw.args[1]] = kw.args[2]
        else
            error("@dist: expected keyword arguments (e.g. mean=5.0), got $kw")
        end
    end

    has_mean = haskey(kw_dict, :mean)
    has_var = haskey(kw_dict, :var)
    has_std = haskey(kw_dict, :std)
    has_cv = haskey(kw_dict, :cv)
    has_scv = haskey(kw_dict, :scv)
    has_second_moment = haskey(kw_dict, :second_moment)
    has_support = haskey(kw_dict, :support)

    # Determine variance from whichever dispersion keyword is given
    var_expr = if has_var
        kw_dict[:var]
    elseif has_std
        :($(kw_dict[:std])^2)
    elseif has_cv && has_mean
        :(($(kw_dict[:cv]) * $(kw_dict[:mean]))^2)
    elseif has_scv && has_mean
        :($(kw_dict[:scv]) * $(kw_dict[:mean])^2)
    elseif has_second_moment && has_mean
        :($(kw_dict[:second_moment]) - $(kw_dict[:mean])^2)
    else
        nothing
    end

    mean_expr = has_mean ? kw_dict[:mean] : nothing
    support_expr = has_support ? kw_dict[:support] : nothing

    if expr isa Symbol
        # Form 1: bare type — @dist Gamma mean=5.0 var=3.0
        D = expr
        return _emit_type_call(D, mean_expr, var_expr, support_expr)

    elseif expr isa Expr && expr.head == :call
        D = expr.args[1]
        args = expr.args[2:end]

        has_placeholder = any(a -> a == :_ || a == :(_), args)

        if has_placeholder
            # Form 2: partial — @dist Gamma(3.0, _) mean=5.0
            return _emit_partial_call(D, args, mean_expr, var_expr, support_expr)
        else
            # Form 3: full instance — @dist TDist(7) mean=5.0 var=2.0
            return _emit_instance_call(D, args, mean_expr, var_expr, support_expr)
        end
    else
        error("@dist: expected a distribution type or constructor call, got $expr")
    end
end

function _emit_type_call(D, mean_expr, var_expr, support_expr)
    if mean_expr !== nothing && var_expr !== nothing
        if support_expr !== nothing
            return esc(:(dist_from_mean_var_on_support($D, $mean_expr, $var_expr; support=$support_expr)))
        else
            return esc(:(dist_from_mean_var($D, $mean_expr, $var_expr)))
        end
    elseif mean_expr !== nothing
        return esc(:(dist_from_mean($D, $mean_expr)))
    elseif var_expr !== nothing
        return esc(:(dist_from_var($D, $var_expr)))
    else
        error("@dist: must specify at least one moment (mean, var, std, cv, scv, or second_moment)")
    end
end

function _emit_partial_call(D, args, mean_expr, var_expr, support_expr)
    # Build NamedTuple expression: map positional args to fieldnames at runtime
    param_exprs = []
    for (i, a) in enumerate(args)
        if a == :_ || a == :(_)
            push!(param_exprs, :missing)
        else
            push!(param_exprs, a)
        end
    end

    # Generate: PartialDist{D}(NamedTuple{fieldnames(D)}((args...)))
    params_tuple = Expr(:tuple, param_exprs...)
    partial_expr = :(PartialDist{$D}(NamedTuple{fieldnames($D)}($params_tuple)))

    if mean_expr !== nothing && var_expr !== nothing
        return esc(:(dist_from_mean_var($partial_expr, $mean_expr, $var_expr)))
    elseif mean_expr !== nothing
        return esc(:(dist_from_mean($partial_expr, $mean_expr)))
    else
        error("@dist with partial specification: must provide at least mean=...")
    end
end

function _emit_instance_call(D, args, mean_expr, var_expr, support_expr)
    instance_expr = Expr(:call, D, args...)

    if mean_expr !== nothing && var_expr !== nothing
        if support_expr !== nothing
            return esc(:(dist_from_mean_var_on_support($instance_expr, $mean_expr, $var_expr; support=$support_expr)))
        else
            return esc(:(dist_from_mean_var($instance_expr, $mean_expr, $var_expr)))
        end
    else
        # No moments — just construct the instance
        return esc(instance_expr)
    end
end
