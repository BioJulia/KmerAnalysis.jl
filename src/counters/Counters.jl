module Counters

export serial_mem
    

using
    BioSequences,
    ReadDatastores,
    Distributed

import ..MerCounting: collapse_into_counts, CountMode, MerCount, unsafe_merge_into!

include("serial_mem.jl")
include("dist_mem.jl")


end # module