using Plots
using SpecialFunctions

# The equation to solve: C = γ / B(1/γ, 1/γ)
# B(1/γ,1/γ) = Γ(1/γ)²/Γ(2/γ)
# So: C = γ * Γ(2/γ) / Γ(1/γ)²

f(γ) = γ * gamma(2/γ) / gamma(1/γ)^2
logf(γ) = log(γ) + loggamma(2/γ) - 2*loggamma(1/γ)

# --- Plot 1: log-scale C vs γ  ---
γ_range = 0.05:0.01:20.0
C_values = f.(γ_range)

p1 = plot(γ_range, C_values, lw=2, label="C(γ)", yscale=:log10,
    xlabel="γ", ylabel="C  (log scale)",
    title="C = γ / B(1/γ, 1/γ)")
hline!([0.5], ls=:dash, lc=:red, label="C = 1/2 (limit as γ→∞)")
hline!([1.0], ls=:dash, lc=:green, label="C = 1  (γ = 1)")

# --- Plot 2: log(C) vs γ (the log version) ---
logC_values = logf.(γ_range)

p2 = plot(γ_range, logC_values, lw=2, label="log C(γ)",
    xlabel="γ", ylabel="log(C)",
    title="log(C) = log(γ) + logΓ(2/γ) - 2logΓ(1/γ)")
hline!([log(0.5)], ls=:dash, lc=:red, label="log(1/2) ≈ -0.693")
hline!([0.0], ls=:dash, lc=:green, label="log(1) = 0")

# --- Plot 3: zoomed into the "action" region γ ∈ [0.3, 10] ---
γ_zoom = 0.3:0.005:10.0
C_zoom = f.(γ_zoom)

p3 = plot(γ_zoom, C_zoom, lw=2, label="C(γ)",
    xlabel="γ", ylabel="C",
    title="C vs γ  (zoomed, γ ∈ [0.3, 10])")
hline!([0.5], ls=:dash, lc=:red, label="C = 1/2")
hline!([1.0], ls=:dash, lc=:green, label="C = 1")

# Mark some reference points
scatter!([1.0], [f(1.0)], ms=5, mc=:green, label="γ=1 → C=1")
scatter!([0.5], [f(0.5)], ms=5, mc=:orange, label="γ=0.5 → C=3")
scatter!([2.0], [f(2.0)], ms=5, mc=:blue, label="γ=2 → C≈0.637")

# --- Plot 4: inverse view — γ vs C (what we solve for) ---
# Use the same data but flip axes, and use log scale for C
p4 = plot(C_zoom, γ_zoom, lw=2, label="γ(C)", xscale=:log10,
    xlabel="C  (log scale)", ylabel="γ",
    title="Inverse: solving for γ given C")
vline!([0.5], ls=:dash, lc=:red, label="C = 1/2 (γ → ∞)")
vline!([1.0], ls=:dash, lc=:green, label="C = 1 (γ = 1)")

p_all = plot(p1, p2, p3, p4, layout=(2,2), size=(1100,850),
    plot_title="Exploring C = γ / B(1/γ, 1/γ)", plot_titlefontsize=11)
savefig(p_all, "temp/beta_eq_exploration.pdf")

# --- Key values table ---
println("=== Key values of f(γ) = γ / B(1/γ, 1/γ) ===")
println("  γ → ∞ :  C → 1/2")
for γ in [0.1, 0.2, 0.25, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0]
    println("  γ = $(lpad(γ,5))  =>  C = $(round(f(γ), sigdigits=6)),  log(C) = $(round(logf(γ), sigdigits=6))")
end

println("\n=== Properties ===")
println("  f is strictly decreasing on (0, ∞)")
println("  Range of C: (1/2, ∞)")
println("  f(1) = 1 exactly (since B(1,1) = 1)")
println("  As γ→0⁺: C→∞")
println("  As γ→∞ : C→1/2")
