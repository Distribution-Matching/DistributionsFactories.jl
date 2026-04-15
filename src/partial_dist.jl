"""
    DistSpec{D}

A partial distribution specification where some parameters are fixed and others
are `missing` (to be solved from moment constraints).

# Fields
- `params::NamedTuple` — parameter values, with `missing` for parameters to be fitted

# Examples
```julia
DistSpec{Gamma}((α=3.0, θ=missing))   # fix α, fit θ
DistSpec{Logistic}((μ=2.0, θ=missing)) # fix μ, fit θ
```
"""
struct DistSpec{D<:Distribution, NT<:NamedTuple}
    params::NT
end

DistSpec{D}(nt::NT) where {D<:Distribution, NT<:NamedTuple} = DistSpec{D, NT}(nt)

"""
    fixed_params(p::DistSpec)

Return a `NamedTuple` of the parameters that are fixed (not `missing`).

```julia
fixed_params(@dist Gamma(3.0, _))   # (α = 3.0,)
```
"""
function fixed_params(p::DistSpec)
    keys_fixed = [k for k in keys(p.params) if !ismissing(p.params[k])]
    vals_fixed = [p.params[k] for k in keys_fixed]
    return NamedTuple{Tuple(keys_fixed)}(Tuple(vals_fixed))
end

"""
    free_params(p::DistSpec)

Return a `Tuple` of the parameter names that are `missing` (to be solved).

```julia
free_params(@dist Gamma(3.0, _))   # (:θ,)
```
"""
function free_params(p::DistSpec)
    return Tuple(k for k in keys(p.params) if ismissing(p.params[k]))
end

function _param_names(::Type{D}) where {D<:Distribution}
    return fieldnames(D)
end

"""
    dist_from_mean(p::DistSpec, μ̄)

Construct a distribution from a `DistSpec` with one `missing` parameter,
solving it from the mean `μ̄`.
"""
function dist_from_mean(p::DistSpec{D}, μ̄::Number) where {D}
    fixed, free = _split_params(p)
    length(free) == 1 || throw(ArgumentError(
        "dist_from_mean with DistSpec requires exactly 1 missing parameter, got $(length(free)): $free"))

    free_name = first(free)

    sol = _solve_for_param(p, free_name,
        x -> mean(D(_fill_params(p, free_name, x)...)) - μ̄,
        _initial_guess(D, free_name, μ̄),
        "mean=$μ̄")

    return D(_fill_params(p, free_name, sol)...)
end

"""
    dist_from_mean_var(p::DistSpec, μ̄, σ̄²)

Construct a distribution from a `DistSpec`. If all parameters are `missing`,
delegates to the standard type-based method. If some are fixed, solves the
remaining from moment constraints.
"""
function dist_from_mean_var(p::DistSpec{D}, μ̄::Number, σ̄²::Number) where {D}
    fixed, free = _split_params(p)

    if length(free) == 0
        # All params given — construct directly, wrap in LocationScale if moments don't match
        d = D(values(p.params)...)
        if isapprox(mean(d), μ̄, rtol=1e-8) && isapprox(var(d), σ̄², rtol=1e-8)
            return d
        else
            # Wrap in LocationScale to achieve target moments
            σ_scale = √(σ̄² / var(d))
            μ_shift = μ̄ - σ_scale * mean(d)
            return μ_shift + σ_scale * d
        end
    elseif length(free) == length(_param_names(D))
        # All missing — delegate to standard type dispatch
        return dist_from_mean_var(D, μ̄, σ̄²)
    elseif length(free) == 1
        # One free parameter, two constraints (overdetermined).
        # Try solving from variance first (scale params affect variance),
        # then verify mean. If that fails, solve from mean and verify variance.
        free_name = first(free)
        guess = _initial_guess(D, free_name, √σ̄²)

        # Attempt 1: solve from variance, verify mean
        local d
        try
            sol = _solve_for_param(p, free_name,
                x -> var(D(_fill_params(p, free_name, x)...)) - σ̄²,
                guess, "var=$σ̄²")
            d = D(_fill_params(p, free_name, sol)...)
            if isapprox(mean(d), μ̄, rtol=1e-4)
                return d
            end
        catch
        end

        # Attempt 2: solve from mean, verify variance
        try
            sol = _solve_for_param(p, free_name,
                x -> mean(D(_fill_params(p, free_name, x)...)) - μ̄,
                guess, "mean=$μ̄")
            d = D(_fill_params(p, free_name, sol)...)
            if isapprox(var(d), σ̄², rtol=1e-4)
                return d
            end
        catch
        end

        throw(DomainError("$D with fixed $(Dict(fixed)): cannot satisfy both μ̄=$μ̄ and σ̄²=$σ̄² with 1 free parameter"))
    elseif length(free) == 2
        # Two free parameters — solve from mean and variance via 2D Newton
        free_names = collect(free)
        h = 1e-7
        x = [_initial_guess(D, free_names[1], μ̄), _initial_guess(D, free_names[2], √σ̄²)]

        function residual(x1, x2)
            ps = _fill_params(p, free_names, [x1, x2])
            d = D(ps...)
            return [mean(d) - μ̄, var(d) - σ̄²]
        end

        for _ in 1:200
            F = residual(x[1], x[2])
            maximum(abs.(F)) < 1e-10 && break

            J = zeros(2, 2)
            for j in 1:2
                xp = copy(x); xp[j] += h
                Fp = residual(xp[1], xp[2])
                J[:, j] = (Fp - F) / h
            end

            dx = J \ (-F)
            step = 1.0
            for _ in 1:20
                x_new = x + step * dx
                try
                    F_new = residual(x_new[1], x_new[2])
                    if maximum(abs.(F_new)) < maximum(abs.(F))
                        break
                    end
                catch
                end
                step *= 0.5
            end
            x .+= step * dx
        end

        full_params = _fill_params(p, free_names, x)
        return D(full_params...)
    else
        throw(ArgumentError("Cannot solve $D with $(length(free)) free parameters from 2 moment constraints"))
    end
end

"""
    dist_from_var(p::DistSpec, σ̄²)

Construct a distribution from a `DistSpec` with one `missing` parameter,
solving it from the variance `σ̄²`.
"""
function dist_from_var(p::DistSpec{D}, σ̄²::Number) where {D}
    fixed, free = _split_params(p)
    length(free) == 1 || throw(ArgumentError(
        "dist_from_var with DistSpec requires exactly 1 missing parameter, got $(length(free)): $free"))

    free_name = first(free)

    sol = _solve_for_param(p, free_name,
        x -> var(D(_fill_params(p, free_name, x)...)) - σ̄²,
        _initial_guess(D, free_name, √σ̄²),
        "var=$σ̄²")

    return D(_fill_params(p, free_name, sol)...)
end

# Convenience wrappers for DistSpec — delegate to dist_from_mean_var or dist_from_var

dist_from_std(p::DistSpec, σ̄::Number) = dist_from_var(p, σ̄^2)
dist_from_mean_std(p::DistSpec, μ̄::Number, σ̄::Number) = dist_from_mean_var(p, μ̄, σ̄^2)
dist_from_mean_cv(p::DistSpec, μ̄::Number, cv::Number) = dist_from_mean_var(p, μ̄, (cv * μ̄)^2)
dist_from_mean_scv(p::DistSpec, μ̄::Number, scv::Number) = dist_from_mean_var(p, μ̄, scv * μ̄^2)
dist_from_mean_second_moment(p::DistSpec, μ̄::Number, m2::Number) = dist_from_mean_var(p, μ̄, m2 - μ̄^2)

# --- Internal helpers ---

function _split_params(p::DistSpec{D}) where {D}
    names = _param_names(D)
    fixed = [(n => p.params[n]) for n in keys(p.params) if !ismissing(p.params[n])]
    free = [n for n in keys(p.params) if ismissing(p.params[n])]
    return fixed, free
end

function _fill_params(p::DistSpec{D}, free_name::Symbol, value) where {D}
    names = keys(p.params)
    return Tuple(n == free_name ? value : p.params[n] for n in names)
end

function _fill_params(p::DistSpec{D}, free_names::Vector{Symbol}, values) where {D}
    names = keys(p.params)
    free_map = Dict(zip(free_names, values))
    return Tuple(haskey(free_map, n) ? free_map[n] : p.params[n] for n in names)
end

function _initial_guess(::Type{D}, param_name::Symbol, hint::Number) where {D}
    # Generic initial guess — use hint magnitude, stay positive
    return max(abs(hint), 0.1)
end

function _solve_for_param(p::DistSpec{D}, free_name::Symbol, objective::Function, guess::Number, target_desc::String) where {D}
    try
        return find_zero(objective, guess)
    catch e
        if e isa Roots.ConvergenceFailed
            fixed = Dict(k => p.params[k] for k in keys(p.params) if !ismissing(p.params[k]))
            throw(DomainError("$D with fixed $fixed: cannot solve free parameter :$free_name to match target $target_desc. The free parameter may not affect this moment."))
        end
        rethrow()
    end
end



"""
    @dist expr

Create a distribution specification. Returns either a `Type`, a distribution instance,
or a `DistSpec` depending on the form:

- `@dist Gamma` → `Gamma` (the type itself)
- `@dist TDist(7)` → `TDist(7)` (a full instance)
- `@dist Gamma(3.0, _)` → `DistSpec{Gamma}((α=3.0, θ=missing))` (partial spec)

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
            # Partial: @dist Gamma(3.0, _) → DistSpec{Gamma}(...)
            param_exprs = [a == :_ || a == :(_) ? :missing : a for a in args]
            params_tuple = Expr(:tuple, param_exprs...)
            return esc(:(DistSpec{$D}(NamedTuple{fieldnames($D)}($params_tuple))))
        else
            # Full instance: @dist TDist(7) → TDist(7)
            return esc(expr)
        end
    else
        error("@dist: expected a distribution type or constructor call, got $expr")
    end
end
