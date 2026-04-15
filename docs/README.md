# Building documentation locally

```bash
# From the repository root:
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The output is in `docs/build/`. Open `docs/build/index.html` in a browser.

# Deploying to GitHub Pages

Documentation is automatically built and deployed by GitHub Actions on every push to `main`.

To enable, go to the GitHub repo settings:
**Settings > Pages > Source: Deploy from a branch > Branch: `gh-pages`**

The docs will be available at:
`https://ron-ash.github.io/DistributionsFactories.jl/`
