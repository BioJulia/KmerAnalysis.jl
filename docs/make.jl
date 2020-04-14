using Documenter, MerCounting, Pkg

makedocs(
    modules = [MerCounting],
    format = Documenter.HTML(),
    sitename = "MerCounting.jl",
    authors = replace(join(Pkg.TOML.parsefile("Project.toml")["authors"], ", "), r" <.*?>" => "" ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
        ],
        "API" => [
            "MerCount" => "api/MerCount.md"
        ]
    ],
    
)

deploydocs(
    repo = "github.com/BioJulia/GenomeGraphs.jl.git",
    push_preview = true,
    deps = nothing,
    make = nothing
)