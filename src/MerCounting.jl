module MerCounting

export
    CountMode,
    Canonical,
    NonCanonical,
    IndexedCounts,
    
    MerCount,
    mer,
    freq,
    merge,
    merge_into!,
    collapse_into_counts,
    collapse_into_counts!,
    
    # Single thread in RAM kmer counter.
    count_mers_srd

using
    BioSequences,
    ReadDatastores,
    Distributed
    
abstract type CountMode end

struct Canonical    <: CountMode end
struct NonCanonical <: CountMode end

@inline (::Canonical)(x) = canonical(x)
@inline (::NonCanonical)(x) = fwmer(x)

const CANONICAL = Canonical()
const NONCANONICAL = NonCanonical()

include("MerCount.jl")
include("counters/Counters.jl")
include("counter/MinimizerTable.jl")
include("MerCountHist.jl")
include("IndexedCounts.jl")

end # module
