#!/usr/bin/env julia
# Build a JSON oracle of (dist, kwargs) -> (mean, var, ppf, params) used by
# the Python (distsfactory) and R (distsfactory) downstream packages to verify
# that their make_dist / dist_exists implementations agree numerically with
# this Julia reference.
#
# Usage (from the DistributionsFactories.jl repo root):
#   julia --project=. scripts/build_cross_oracle.jl
#
# By default, writes the oracle JSON to the checked-in test-data locations of
# the sibling repos when they exist:
#   ../distsfactory-python/tests/data/cross_oracle.json
#   ../distsfactory-r/tests/testthat/data/cross_oracle.json
#
# Both consumers read the same schema. R filters cases at load time to the
# families it implements.
#
# Optional first argument: a single output path to write to instead of the
# sibling defaults.
#
# The DistributionsFactories.jl project needs JSON3 in its Project.toml:
#   julia --project=. -e 'using Pkg; Pkg.add("JSON3")'
#
# To add cases, append to the CASES vector below. Each entry is a Case with
# the canonical distribution name (shared across Python and R), the Julia
# type, and the make_dist kwargs as a NamedTuple. Discrete distributions emit
# moments only; continuous distributions also emit ppf at quantiles
# (0.1, 0.5, 0.9).

using DistributionsFactories
using Distributions
using JSON3

struct Case
    py_name::String          # canonical Python distribution name
    julia_type::Type         # Julia type to pass to make_dist
    kwargs::NamedTuple       # kwargs to pass to make_dist (must round-trip to Python identically)
end

# A feasibility check oracle entry: same shape as a constructor case, but
# records `dist_exists(...)`'s bool answer instead of moments. Used by the
# Python pytest to verify both packages agree on the feasibility region.
struct FeasibilityCase
    py_name::String
    julia_type::Any            # Type or Distribution instance
    kwargs::NamedTuple
    julia_support::Any         # nothing, an Interval, or a UnitRange
    py_support::Any            # nothing, [lo, hi] tuple, or {start, stop} dict
end

# ---------------------------------------------------------------------------
# Battery of test cases
# ---------------------------------------------------------------------------

const CASES = Case[
    # Real-line continuous: mean+var
    Case("normal",         Normal,             (mean=0.0,  var=1.0)),
    Case("normal",         Normal,             (mean=5.0,  var=3.0)),
    Case("normal",         Normal,             (mean=-2.0, var=4.0)),
    Case("laplace",        Laplace,            (mean=0.0,  var=2.0)),
    Case("laplace",        Laplace,            (mean=5.0,  var=8.0)),
    Case("logistic",       Logistic,           (mean=0.0,  var=π^2/3)),
    Case("logistic",       Logistic,           (mean=3.0,  var=10.0)),
    Case("gumbel",         Gumbel,             (mean=5.0,  var=3.0)),
    Case("gumbel",         Gumbel,             (mean=0.0,  var=π^2/6)),
    Case("uniform",        Uniform,            (mean=5.0,  var=3.0)),
    Case("uniform",        Uniform,            (mean=0.0,  var=1/12)),
    Case("sym_triangular", SymTriangularDist,  (mean=2.0,  var=6.0)),

    # TDist
    Case("tdist",          TDist,              (mean=0.0,  var=2.0)),
    Case("tdist",          TDist,              (mean=0.0,  var=5.0)),

    # Positive continuous: mean+var
    Case("gamma",          Gamma,              (mean=5.0,  var=3.0)),
    Case("gamma",          Gamma,              (mean=10.0, var=50.0)),
    Case("gamma",          Gamma,              (mean=0.5,  var=0.1)),
    Case("exponential",    Exponential,        (mean=3.0,  var=9.0)),
    Case("erlang",         Erlang,             (mean=4.0,  var=4.0)),  # k=4, θ=1
    Case("lognormal",      LogNormal,          (mean=5.0,  var=3.0)),
    Case("lognormal",      LogNormal,          (mean=2.0,  var=10.0)),
    Case("weibull",        Weibull,            (mean=5.0,  var=3.0)),
    Case("weibull",        Weibull,            (mean=2.0,  var=10.0)),
    Case("frechet",        Frechet,            (mean=5.0,  var=3.0)),
    Case("rayleigh",       Rayleigh,           (mean=5.0,  var=5.0^2 * (4-π)/π)),
    Case("inverse_gamma",  InverseGamma,       (mean=5.0,  var=3.0)),
    Case("pareto",         Pareto,             (mean=5.0,  var=3.0)),
    Case("folded_normal",  FoldedNormal,       (mean=2.5,  var=1.2)),
    Case("fdist",          FDist,              (mean=1.5,  var=10.0)),
    Case("chi",            Chi,                # Build a known feasible point from nu=4
                                                (mean=mean(Chi(4)), var=var(Chi(4)))),
    Case("chisq",          Chisq,              (mean=4.0,  var=8.0)),

    # Unit-interval continuous
    Case("beta",           Beta,               (mean=0.5,  var=0.05)),
    Case("beta",           Beta,               (mean=0.3,  var=0.02)),
    Case("beta",           Beta,               (mean=0.8,  var=0.04)),

    # Discrete
    Case("binomial",       Binomial,           (mean=5.0,  var=2.5)),
    Case("poisson",        Poisson,            (mean=5.0,  var=5.0)),
    Case("negative_binomial", NegativeBinomial,(mean=5.0,  var=8.0)),
    Case("geometric",      Geometric,          (mean=2.0,  var=2.0*(1+2.0))),
    Case("discrete_uniform", DiscreteUniform,  (mean=5.0,  var=10.0)),

    # Two-quantile constructions
    Case("normal",         Normal,             (q1=-1.0,   q3=1.0)),
    Case("logistic",       Logistic,           (q1=-1.0,   q3=1.0)),
    Case("laplace",        Laplace,            (q1=-1.0,   q3=1.0)),
    Case("lognormal",      LogNormal,          (quantiles=[(0.1, 1.0), (0.9, 10.0)],)),
    Case("weibull",        Weibull,            (quantiles=[(0.1, 1.0), (0.9, 5.0)],)),
    Case("pareto",         Pareto,             (quantiles=[(0.1, 1.0), (0.9, 10.0)],)),
    Case("gamma",          Gamma,              (q1=2.0,    q3=8.0)),
    Case("beta",           Beta,               (q1=0.3,    q3=0.7)),

    # Mean + quantile
    Case("gamma",          Gamma,              (mean=5.0,  median=4.5)),
    Case("beta",           Beta,               (mean=0.4,  median=0.38)),
    Case("lognormal",      LogNormal,          (mean=3.0,  median=2.5)),

    # Mode-based
    Case("gamma",          Gamma,              (mean=5.0,  mode=3.0)),
    Case("gamma",          Gamma,              (mode=3.0,  iqr=4.0)),
    Case("beta",           Beta,               (mean=0.4,  mode=0.35)),
    Case("rayleigh",       Rayleigh,           (mode=2.0,)),
    Case("normal",         Normal,             (mode=3.0,  var=4.0)),

    # 1-parameter from mean
    Case("exponential",    Exponential,        (mean=3.0,)),
    Case("poisson",        Poisson,            (mean=5.0,)),
    Case("rayleigh",       Rayleigh,           (mean=5.0,)),
    Case("geometric",      Geometric,          (mean=3.0,)),

    # Single-quantile
    Case("exponential",    Exponential,        (median=2.0,)),
    Case("rayleigh",       Rayleigh,           (median=2.0,)),

    # 3-parameter triangular (mean + var + mode)
    Case("triangular",     TriangularDist,        (mean=5.0, var=2.0, mode=4.0)),
    Case("triangular",     TriangularDist,        (mean=0.0, var=1.0, mode=0.0)),
    Case("discrete_sym_triangular", DiscreteSymmetricTriangular, (mean=5.0, var=4.0)),
    Case("discrete_triangular",     DiscreteTriangular,         (mean=5.0, var=2.0, mode=5.0)),
]

# Feasibility-region oracle: each case records the bool that `dist_exists`
# returns on each side. Mix of in-region, on-the-boundary, and out-of-region
# inputs so the Python predicate is forced to match Julia case-by-case.
const FEAS_CASES = FeasibilityCase[
    # Beta: var < mu*(1-mu)
    FeasibilityCase("beta",        Beta,        (mean=0.5,  var=0.1),   nothing, nothing),
    FeasibilityCase("beta",        Beta,        (mean=0.5,  var=0.3),   nothing, nothing),
    FeasibilityCase("beta",        Beta,        (mean=1.5,  var=0.1),   nothing, nothing),

    # Exponential: var must equal mean^2
    FeasibilityCase("exponential", Exponential, (mean=2.5,  var=6.25),  nothing, nothing),
    FeasibilityCase("exponential", Exponential, (mean=2.5,  var=1.5),   nothing, nothing),

    # Poisson: var must equal mean
    FeasibilityCase("poisson",     Poisson,     (mean=5.0,  var=5.0),   nothing, nothing),
    FeasibilityCase("poisson",     Poisson,     (mean=5.0,  var=3.0),   nothing, nothing),

    # NegativeBinomial: var > mean
    FeasibilityCase("negative_binomial", NegativeBinomial, (mean=5.0, var=8.0), nothing, nothing),
    FeasibilityCase("negative_binomial", NegativeBinomial, (mean=5.0, var=3.0), nothing, nothing),

    # TDist: mu = 0, var > 1
    FeasibilityCase("tdist",       TDist,       (mean=0.0,  var=2.0),   nothing, nothing),
    FeasibilityCase("tdist",       TDist,       (mean=1.0,  var=2.0),   nothing, nothing),
    FeasibilityCase("tdist",       TDist,       (mean=0.0,  var=0.5),   nothing, nothing),

    # FDist: 1 < mu < 2, plus a lower-bound on var
    FeasibilityCase("fdist",       FDist,       (mean=1.5,  var=10.0),  nothing, nothing),
    FeasibilityCase("fdist",       FDist,       (mean=0.9,  var=5.0),   nothing, nothing),
    FeasibilityCase("fdist",       FDist,       (mean=2.5,  var=5.0),   nothing, nothing),

    # Cauchy: no finite mean/var ever
    FeasibilityCase("cauchy",      Cauchy,      (mean=0.0,  var=1.0),   nothing, nothing),

    # Truncated Normal on a bounded interval — Langevin dome
    # On [-1, 1] at mu=0: dome boundary ~= 0.333
    FeasibilityCase("normal", Normal, (mean=0.0, var=0.10), -1.0..1.0, (-1.0, 1.0)),
    FeasibilityCase("normal", Normal, (mean=0.0, var=0.30), -1.0..1.0, (-1.0, 1.0)),
    FeasibilityCase("normal", Normal, (mean=0.0, var=0.50), -1.0..1.0, (-1.0, 1.0)),

    # Half-truncated Normal on [0, inf): bound is (mu - lo)^2
    FeasibilityCase("normal", Normal, (mean=2.0, var=3.0), 0.0..Inf,   (0.0, Inf)),
    FeasibilityCase("normal", Normal, (mean=2.0, var=5.0), 0.0..Inf,   (0.0, Inf)),

    # Half-truncated Laplace: boundary is attained
    FeasibilityCase("laplace", Laplace, (mean=2.0, var=4.0), 0.0..Inf, (0.0, Inf)),  # exactly (mu-lo)^2
    FeasibilityCase("laplace", Laplace, (mean=2.0, var=5.0), 0.0..Inf, (0.0, Inf)),
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Translate a Julia kwarg name/value to its JSON-canonical form for Python.
# Tuples in `quantiles` need explicit listification.
function _kwargs_to_json(kw::NamedTuple)
    out = Dict{String,Any}()
    for (k, v) in pairs(kw)
        key = String(k)
        if v isa AbstractVector && eltype(v) <: Tuple
            out[key] = [collect(t) for t in v]
        elseif v isa Tuple
            out[key] = collect(v)
        else
            out[key] = v
        end
    end
    return out
end

# For a Julia distribution, return basic moments + a small ppf table.
# Discrete distributions skip ppf (scipy and Julia disagree on ties; not a
# bug, just a parameterization detail).
_safe(x) = isfinite(x) ? x : (isnan(x) ? "nan" : (x > 0 ? "inf" : "-inf"))

function _dist_summary(d)
    info = Dict{String,Any}(
        "mean" => _safe(mean(d)),
        "var"  => _safe(var(d)),
        "discrete" => d isa DiscreteUnivariateDistribution,
    )
    if !(d isa DiscreteUnivariateDistribution)
        info["ppf"] = Dict(
            "0.1" => _safe(quantile(d, 0.1)),
            "0.5" => _safe(quantile(d, 0.5)),
            "0.9" => _safe(quantile(d, 0.9)),
        )
    end
    return info
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function build(output_path)
    records = []
    for case in CASES
        d = make_dist(case.julia_type; case.kwargs...)
        push!(records, Dict(
            "py_name" => case.py_name,
            "julia_type" => string(case.julia_type),
            "kwargs" => _kwargs_to_json(case.kwargs),
            "expected" => _dist_summary(d),
        ))
    end

    feas_records = []
    for case in FEAS_CASES
        ok = if case.julia_support === nothing
            dist_exists(case.julia_type; case.kwargs...)
        else
            dist_exists(case.julia_type; support=case.julia_support, case.kwargs...)
        end
        py_support = if case.py_support === nothing
            nothing
        else
            [_safe(case.py_support[1]), _safe(case.py_support[2])]
        end
        push!(feas_records, Dict(
            "py_name" => case.py_name,
            "kwargs"  => _kwargs_to_json(case.kwargs),
            "support" => py_support,
            "expected" => ok,
        ))
    end

    mkpath(dirname(output_path))
    open(output_path, "w") do io
        JSON3.pretty(io, Dict("constructors" => records, "feasibility" => feas_records))
    end

    @info "Wrote $(length(records)) constructor cases and $(length(feas_records)) feasibility cases to $output_path"
end

# Run
let
    here = @__DIR__
    if length(ARGS) >= 1
        build(ARGS[1])
    else
        # Default: write to whichever sibling consumer repos exist.
        sibling_targets = [
            normpath(joinpath(here, "..", "..", "distsfactory-python", "tests", "data", "cross_oracle.json")),
            normpath(joinpath(here, "..", "..", "distsfactory-r", "tests", "testthat", "data", "cross_oracle.json")),
        ]
        any_written = false
        for target in sibling_targets
            parent = dirname(target)
            if isdir(dirname(parent))  # tests/data/ — accept if tests/ exists
                build(target)
                any_written = true
            else
                @info "Skipping (sibling repo not found): $target"
            end
        end
        if !any_written
            @warn "No sibling consumer repos found. Re-run with an explicit output path."
        end
    end
end
