using CairoMakie

xs = range(1, 2, length=500)
ys = range(0, 100, length=500)

tol = 1e-10  # mask denominator values close to 0 (adjust if needed)

Zmasked = [begin
    denom = y*(2 - x) - x^2*(x - 1)

    # make holes at singularities (and near-singularities)
    if abs(denom) < tol
        NaN
    else
        z = 2*x^2 / denom
        (isfinite(z) && z >= 0) ? min(z, 10) : NaN   # clamp to [0,10], otherwise hole
    end
end for x in xs, y in ys]

fig = Figure(size=(900, 700))
fig = Figure(size=(900, 700))
ax = Axis3(fig[1, 1],
    xlabel="μ", ylabel="var", zlabel="v₁",
    azimuth = 1.25,
    elevation = 0.35,
    title = "Feasible Region of v₁ defined by μ and σ²"
)

surface!(ax, xs, ys, Zmasked; colormap=:viridis, nan_color=:transparent)
save("v1.png", fig) 