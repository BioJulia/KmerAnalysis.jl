module KmerAnalysis

export
    CountMode,
    Canonical,
    CANONICAL,
    NonCanonical,
    NONCANONICAL,
    IndexedCounts,
    
    MerCount,
    mer,
    freq,
    merge,
    merge_into!,
    collapse_into_counts,
    collapse_into_counts!,
    collect_mers,
    
    Counters,
    
    # Kmer frequency spectra
    spectra,
    spectra!
    
    

using BioSequences, ReadDatastores, Distributed
    
abstract type CountMode end

struct Canonical    <: CountMode end
struct NonCanonical <: CountMode end

@inline (::Canonical)(x) = canonical(x)
@inline (::NonCanonical)(x) = fwmer(x)

const CANONICAL = Canonical()
const NONCANONICAL = NonCanonical()

include("MerCount.jl")
include("counters/Counters.jl")
include("KmerFrequencySpectra.jl")
include("IndexedCounts.jl")

end # module
