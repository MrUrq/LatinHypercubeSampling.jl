using Documenter, LatinHypercubeSampling, Random

makedocs(;
        modules = [LatinHypercubeSampling],
        format = Documenter.HTML(   assets = ["assets/favicon.ico"],
                                prettyurls = get(ENV, "CI", nothing) == "true"),
        sitename = "LatinHypercubeSampling.jl",
        strict = true,
        clean = true,
        checkdocs = :none,
        pages = [
                "Home" => "index.md",
                "Manual" => [
                        "man/lhcoptim.md",
                        "man/sublhcoptim.md",
                        "man/categorical.md",
                ]
        ]
)

deploydocs(
    repo = "github.com/MrUrq/LatinHypercubeSampling.jl.git"
)