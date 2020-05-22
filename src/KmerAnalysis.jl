module KmerAnalysis

export
    CountMode,
    Canonical,
    CANONICAL,
    NonCanonical,
    NONCANONICAL,
    
    # MerCount
    MerCount,
    mer,
    freq,
    merge,
    merge_into!,
    collapse_into_counts,
    collapse_into_counts!,
    collect_mers,
    collect_mers!,
    
    # Counters
    serial_mem,
    count!,
    dist_mem,
    
    # K-mer frequency spectra
    KmerFrequencySpectra,
    spectra,
    spectra!,
    find_next_peak,
    find_spectra_peaks,
    
    # Indexed k-mer counts
    IndexedCounts,
    add_count!,
    project_count

using BioSequences, ReadDatastores, Distributed
    
abstract type CountMode end

struct Canonical    <: CountMode end
struct NonCanonical <: CountMode end

@inline (::Canonical)(x) = canonical(x)
@inline (::NonCanonical)(x) = fwmer(x)

const CANONICAL = Canonical()
const NONCANONICAL = NonCanonical()

abstract type AbstractKmerCounter{M<:AbstractMer} <: Function end

include("MerCount.jl")
include("counters/serial_mem.jl")
include("counters/dist_mem.jl")
include("KmerFrequencySpectra.jl")
include("IndexedCounts.jl")

end # module
