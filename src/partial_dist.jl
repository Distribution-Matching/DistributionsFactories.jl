"""
    PartialDist{D}

A partial distribution specification where some parameters are fixed and others
are `missing` (to be solved from moment constraints).

# Fields
- `params::NamedTuple` — parameter values, with `missing` for parameters to be fitted

# Examples
```julia
PartialDist{Gamma}((α=3.0, θ=missing))   # fix α, fit θ
PartialDist{Logistic}((μ=2.0, θ=missing)) # fix μ, fit θ
```
"""
struct PartialDist{D<:Distribution, NT<:NamedTuple}
    params::NT
end

PartialDist{D}(nt::NT) where {D<:Distribution, NT<:NamedTuple} = PartialDist{D, NT}(nt)

function _param_names(::Type{D}) where {D<:Distribution}
    return fieldnames(D)
end

"""
    dist_from_mean(p::PartialDist, μ̄)

Construct a distribution from a `PartialDist` with one `missing` parameter,
solving it from the mean `μ̄`.
"""
function dist_from_mean(p::PartialDist{D}, μ̄::Number) where {D}
    fixed, free = _split_params(p)
    length(free) == 1 || throw(ArgumentError(
        "dist_from_mean with PartialDist requires exactly 1 missing parameter, got $(length(free)): $free"))

    free_name = first(free)

    # Solve: find the value of the free parameter such that mean(D(params...)) == μ̄
    sol = find_zero(
        x -> begin
            full_params = _fill_params(p, free_name, x)
            mean(D(full_params...)) - μ̄
        end,
        _initial_guess(D, free_name, μ̄)
    )

    full_params = _fill_params(p, free_name, sol)
    return D(full_params...)
end

"""
    dist_from_mean_var(p::PartialDist, μ̄, σ̄²)

Construct a distribution from a `PartialDist`. If all parameters are `missing`,
delegates to the standard type-based method. If some are fixed, solves the
remaining from moment constraints.
"""
function dist_from_mean_var(p::PartialDist{D}, μ̄::Number, σ̄²::Number) where {D}
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

        # Attempt 1: solve from variance
        local d
        try
            sol = find_zero(
                x -> begin
                    fp = _fill_params(p, free_name, x)
                    var(D(fp...)) - σ̄²
                end, guess)
            d = D(_fill_params(p, free_name, sol)...)
            if isapprox(mean(d), μ̄, rtol=1e-4)
                return d
            end
        catch
        end

        # Attempt 2: solve from mean
        sol = find_zero(
            x -> begin
                fp = _fill_params(p, free_name, x)
                mean(D(fp...)) - μ̄
            end, guess)
        d = D(_fill_params(p, free_name, sol)...)
        if !isapprox(var(d), σ̄², rtol=1e-4)
            throw(DomainError("$D with fixed $(Dict(fixed)): cannot satisfy both μ̄=$μ̄ and σ̄²=$σ̄² with 1 free parameter"))
        end
        return d
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
    dist_from_var(p::PartialDist, σ̄²)

Construct a distribution from a `PartialDist` with one `missing` parameter,
solving it from the variance `σ̄²`.
"""
function dist_from_var(p::PartialDist{D}, σ̄²::Number) where {D}
    fixed, free = _split_params(p)
    length(free) == 1 || throw(ArgumentError(
        "dist_from_var with PartialDist requires exactly 1 missing parameter, got $(length(free)): $free"))

    free_name = first(free)
    sol = find_zero(
        x -> begin
            full_params = _fill_params(p, free_name, x)
            var(D(full_params...)) - σ̄²
        end,
        _initial_guess(D, free_name, √σ̄²)
    )

    full_params = _fill_params(p, free_name, sol)
    return D(full_params...)
end

# Convenience wrappers for PartialDist — delegate to dist_from_mean_var or dist_from_var

dist_from_std(p::PartialDist, σ̄::Number) = dist_from_var(p, σ̄^2)
dist_from_mean_std(p::PartialDist, μ̄::Number, σ̄::Number) = dist_from_mean_var(p, μ̄, σ̄^2)
dist_from_mean_cv(p::PartialDist, μ̄::Number, cv::Number) = dist_from_mean_var(p, μ̄, (cv * μ̄)^2)
dist_from_mean_scv(p::PartialDist, μ̄::Number, scv::Number) = dist_from_mean_var(p, μ̄, scv * μ̄^2)
dist_from_mean_second_moment(p::PartialDist, μ̄::Number, m2::Number) = dist_from_mean_var(p, μ̄, m2 - μ̄^2)

# --- Internal helpers ---

function _split_params(p::PartialDist{D}) where {D}
    names = _param_names(D)
    fixed = [(n => p.params[n]) for n in keys(p.params) if !ismissing(p.params[n])]
    free = [n for n in keys(p.params) if ismissing(p.params[n])]
    return fixed, free
end

function _fill_params(p::PartialDist{D}, free_name::Symbol, value) where {D}
    names = keys(p.params)
    return Tuple(n == free_name ? value : p.params[n] for n in names)
end

function _fill_params(p::PartialDist{D}, free_names::Vector{Symbol}, values) where {D}
    names = keys(p.params)
    free_map = Dict(zip(free_names, values))
    return Tuple(haskey(free_map, n) ? free_map[n] : p.params[n] for n in names)
end

function _initial_guess(::Type{D}, param_name::Symbol, hint::Number) where {D}
    # Generic initial guess — use hint magnitude, stay positive
    return max(abs(hint), 0.1)
end
