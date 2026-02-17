using Plots
using SpecialFunctions
using Distributions
using DistributionsFactories

# ============================================================
# Explore Weibull dist_from_mean_var numerical stability
# across a wide range of CV (coefficient of variation) values.
#
# CV depends only on the shape parameter k (not on λ),
# so we fix λ = 1 and sweep k.
# ============================================================

# --- helpers ---
function weibull_cv(k)
    # CV = √(Γ(1+2/k)/Γ(1+1/k)² - 1)
    return sqrt(gamma(1 + 2/k) / gamma(1 + 1/k)^2 - 1)
end

function try_roundtrip(μ, v, method)
    try
        d = dist_from_mean_var(Weibull, μ, v; method=method)
        return params(d)  # (k, λ)
    catch e
        return nothing
    end
end

# --- parameter sweep ---
k_values = vcat(0.05:0.05:1.0, 1.1:0.1:3.0, 3.5:0.5:10.0)
λ_true = 1.0

# storage
cv_vals      = Float64[]
k_true_vals  = Float64[]
k_sb_vals    = Union{Float64, Nothing}[]
k_og_vals    = Union{Float64, Nothing}[]
err_sb_vals  = Union{Float64, Nothing}[]
err_og_vals  = Union{Float64, Nothing}[]

println("=" ^ 105)
println(rpad("k_true", 8), rpad("CV", 10), rpad("k_sb", 14), rpad("k_og", 14),
        rpad("relerr_sb", 14), rpad("relerr_og", 14), "status")
println("-" ^ 105)

for k in k_values
    d_true = Weibull(k, λ_true)
    μ = mean(d_true)
    v = var(d_true)
    cv = sqrt(v) / μ

    push!(cv_vals, cv)
    push!(k_true_vals, k)

    # simple bracketing
    res_sb = try_roundtrip(μ, v, :simple_bracketing)
    if res_sb !== nothing
        k_sb, λ_sb = res_sb
        err_sb = abs(k_sb - k) / k
        push!(k_sb_vals, k_sb)
        push!(err_sb_vals, err_sb)
    else
        push!(k_sb_vals, nothing)
        push!(err_sb_vals, nothing)
    end

    # oscar garcia
    res_og = try_roundtrip(μ, v, :oscar_garcia)
    if res_og !== nothing
        k_og, λ_og = res_og
        err_og = abs(k_og - k) / k
        push!(k_og_vals, k_og)
        push!(err_og_vals, err_og)
    else
        push!(k_og_vals, nothing)
        push!(err_og_vals, nothing)
    end

    # status
    sb_str = res_sb !== nothing ? "ok" : "FAIL"
    og_str = res_og !== nothing ? "ok" : "FAIL"
    status = "sb=$sb_str  og=$og_str"

    k_sb_s = res_sb !== nothing ? string(round(res_sb[1], sigdigits=8)) : "---"
    k_og_s = res_og !== nothing ? string(round(res_og[1], sigdigits=8)) : "---"
    err_sb_s = err_sb_vals[end] !== nothing ? string(round(err_sb_vals[end], sigdigits=4)) : "---"
    err_og_s = err_og_vals[end] !== nothing ? string(round(err_og_vals[end], sigdigits=4)) : "---"

    println(rpad(round(k, sigdigits=4), 8),
            rpad(round(cv, sigdigits=6), 10),
            rpad(k_sb_s, 14),
            rpad(k_og_s, 14),
            rpad(err_sb_s, 14),
            rpad(err_og_s, 14),
            status)
end

println("=" ^ 105)

# --- plots ---
# filter out failures for plotting
idx_sb = findall(x -> x !== nothing, err_sb_vals)
idx_og = findall(x -> x !== nothing, err_og_vals)

p1 = plot(title="Relative error in recovered k vs CV",
          xlabel="CV", ylabel="Relative error in k",
          yscale=:log10, legend=:topleft, size=(900, 500))

if !isempty(idx_sb)
    plot!(p1, cv_vals[idx_sb], Float64[err_sb_vals[i] for i in idx_sb] .+ 1e-16,
          lw=2, marker=:circle, ms=3, label="simple_bracketing")
end
if !isempty(idx_og)
    plot!(p1, cv_vals[idx_og], Float64[err_og_vals[i] for i in idx_og] .+ 1e-16,
          lw=2, marker=:diamond, ms=3, label="oscar_garcia")
end

# mark failures
idx_sb_fail = findall(x -> x === nothing, err_sb_vals)
idx_og_fail = findall(x -> x === nothing, err_og_vals)

if !isempty(idx_sb_fail)
    vline!(p1, cv_vals[idx_sb_fail], ls=:dash, lc=:red, alpha=0.4, label="sb failure")
end
if !isempty(idx_og_fail)
    vline!(p1, cv_vals[idx_og_fail], ls=:dot, lc=:orange, alpha=0.4, label="og failure")
end

p2 = plot(title="Relative error in recovered k vs k_true",
          xlabel="k_true", ylabel="Relative error in k",
          yscale=:log10, legend=:topleft, size=(900, 500))

if !isempty(idx_sb)
    plot!(p2, k_true_vals[idx_sb], Float64[err_sb_vals[i] for i in idx_sb] .+ 1e-16,
          lw=2, marker=:circle, ms=3, label="simple_bracketing")
end
if !isempty(idx_og)
    plot!(p2, k_true_vals[idx_og], Float64[err_og_vals[i] for i in idx_og] .+ 1e-16,
          lw=2, marker=:diamond, ms=3, label="oscar_garcia")
end

p_all = plot(p1, p2, layout=(2, 1), size=(900, 900))
savefig(p_all, "temp/weibull_stability.pdf")
println("\nPlot saved to temp/weibull_stability.pdf")
