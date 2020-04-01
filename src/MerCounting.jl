module MerCounting

export
    CountMode,
    IndexedCounts,
    
    MerCount,
    mer,
    freq,
    merge,
    collapse_into_counts,
    collapse_into_counts!,
    
    # Single thread in RAM kmer counter.
    count_mers_srd

using
    BioSequences,
    ReadDatastores
    


@enum CountMode Canonical NonCanonical

include("MerCount.jl")
include("counter/MinimizerTable")
include("MerCountHist.jl")
include("IndexedCounts.jl")
include("counters/srd.jl")

end # module
