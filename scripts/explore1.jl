using Plots
using Distributions

mean_points = Float64[]
var_points = Float64[]

for d₁ ∈ 1:0.1:50
    for d₂ ∈ 10:10:5000
        dist = FDist(d₁, d₂)
        push!(mean_points, mean(dist))
        push!(var_points, var(dist))
    end
end

scatter(mean_points, var_points, msw=0, ms= 1, legend=false, xlim=(0,2), ylim = (0,10))