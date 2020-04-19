module Counters

export serial_mem
    

using
    BioSequences,
    ReadDatastores

import ..MerCounting: collapse_into_counts, CountMode

include("serial_mem.jl")
include("dist_mem.jl")


end # module