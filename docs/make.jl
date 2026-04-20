using Documenter
using DistributionsFactories

makedocs(
    sitename = "DistributionsFactories.jl",
    modules = [DistributionsFactories],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
        "Distribution Table" => "distributions.md",
    ],
    warnonly = [:missing_docs],
)

deploydocs(
    repo = "github.com/Distribution-Matching/DistributionsFactories.jl.git",
    devbranch = "main",
)
