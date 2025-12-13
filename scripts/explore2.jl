using Plots
using Distributions

d_1 = Gamma(0.75,1/0.5)
d_2 = FDist(12, 6)
f_1(x) = pdf(d_1, x)
f_2(x) = pdf(d_2, x)
x_range = 0:0.001:5

@show mean(d_1), var(d_1)
@show mean(d_2), var(d_2)

plot(x_range, f_1.(x_range), lw = 2, label = "Gamma density", ylim = (0,2))
plot!(x_range, f_2.(x_range),lw = 2,  label = "F density",
        xlabel = "x", ylabel = "Density")

savefig("temp/gamma_vs_F.pdf")
