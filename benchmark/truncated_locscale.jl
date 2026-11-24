#=
Per-solve timing of truncated_locscale_params / truncated_locscale_moments.

    julia --project=benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
    julia --project=benchmark benchmark/truncated_locscale.jl
=#
using DistributionsFactories, BenchmarkTools, Printf

function report(label, b)
    @printf("%-62s %8.2f µs  %3d allocs\n", label, median(b).time / 1e3, b.allocs)
end

cases = [
    ("Normal, [0, 1], interior",               Normal,      0.0, 1.0, 0.4,  0.03),
    ("Normal, [0, 1], 0.999 of max variance",  Normal,      0.0, 1.0, 0.2,
        0.999 * DistributionsFactories._truncexp_max_var(0.0, 1.0, 0.2)),
    ("Normal, [0, ∞)",                         Normal,      0.0, Inf, 1.5,  1.0),
    ("Normal, (-∞, 2]",                        Normal,     -Inf, 2.0, 0.5,  0.8),
    ("TDist(4.5), [0, 1]",                     TDist(4.5),  0.0, 1.0, 0.4,  0.03),
    ("TDist(4.5), [0, ∞)",                     TDist(4.5),  0.0, Inf, 1.5,  1.0),
    ("TDist(3), (-∞, 2]",                      TDist(3.0), -Inf, 2.0, 0.5,  0.8),
    ("TDist(30), [-1, 2]",                     TDist(30.0), -1.0, 2.0, 0.3, 0.4),
]
println("truncated_locscale_params (one solve, after compilation)")
for (label, fam, a, b, m, v) in cases
    truncated_locscale_params(fam, a, b, m, v)
    report(label, @benchmark truncated_locscale_params($fam, $a, $b, $m, $v))
end
println("\ntruncated_locscale_moments (one evaluation)")
for (label, fam, a, b) in (("Normal, [0, 1]", Normal, 0.0, 1.0),
                          ("Normal, [0, ∞)", Normal, 0.0, Inf),
                          ("TDist(4.5), [0, 1]", TDist(4.5), 0.0, 1.0),
                          ("TDist(4.5), [0, ∞)", TDist(4.5), 0.0, Inf))
    report(label, @benchmark truncated_locscale_moments($fam, 0.3, 0.7, $a, $b))
end
