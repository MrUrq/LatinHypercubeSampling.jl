using Documenter, LatinHypercubeSampling

makedocs(
        modules = [LatinHypercubeSampling],
        format = :html,
        sitename = "LatinHypercubeSampling.jl",
        strict = true,
        assets = ["assets/favicon.ico"],
        clean = true,
        checkdocs = :none,
        pages = Any[
                "Home" => "index.md",
                "Manual" => Any[
                        "man/lhcoptim.md",
                        "man/sublhcoptim.md",
                        "man/refinelhc.md",
                ]
        ]
)

deploydocs(
    repo = "github.com/MrUrq/LatinHypercubeSampling.jl.git",
    target = "build",
    julia = "0.6",
    deps = nothing,
    make = nothing
)
