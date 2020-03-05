module MerCounting

export
    CountMode,
    IndexedCounts

using ReadDatastores


@enum CountMode Canonical NonCanonical

include("IndexedCounts.jl")

end # module
