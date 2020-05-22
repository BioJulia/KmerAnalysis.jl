using Documenter, KmerAnalysis, Pkg

makedocs(
    modules = [KmerAnalysis],
    format = Documenter.HTML(),
    sitename = "KmerAnalysis.jl",
    authors = replace(join(Pkg.TOML.parsefile("Project.toml")["authors"], ", "), r" <.*?>" => "" ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Installation" => "man/installation.md",
            "Basic types and counting" => "man/basic_counting.md",
            "Counters for read datasets" => "man/counting_reads.md",
            "Analysing Kmer count data" => "man/analyses.md"
        ],
        "API" => [
            "MerCount" => "api/MerCount.md",
            "KmerFrequencySpectra" => "api/KmerFrequencySpectra.md"
        ]
    ],
    
)

deploydocs(
    repo = "github.com/BioJulia/KmerAnalysis.jl.git",
    push_preview = true,
    deps = nothing,
    make = nothing
)