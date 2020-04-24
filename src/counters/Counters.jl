module Counters

export
    serial_mem,
    dist_mem

using
    BioSequences,
    ReadDatastores,
    Distributed

import ..MerCounting: collapse_into_counts!, CountMode, MerCount, unsafe_merge_into!, collect_mers

include("serial_mem.jl")
include("dist_mem.jl")
include("dist_batch/dist_batch.jl")


end # module