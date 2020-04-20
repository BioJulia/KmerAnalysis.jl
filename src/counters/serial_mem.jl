
function collect_mers(::Type{M}, count_mode::CountMode, input::DatastoreBuffer{<:PairedReads}, range::AbstractRange = 1:length(input)) where {M<:AbstractMer}
    max_read_size = max_read_length(ReadDatastores.datastore(sbuf))
    v = Vector{M}(undef, length(input) * (max_read_size - ksize(M) + 1))
    wi = firstindex(v)
    read_sequence = eltype(input)()
    @inbounds for i in range
        for mer in each(M, load_sequence!(sbuf, i, read_sequence))
            v[wi] = count_mode(mer)
            wi = wi + 1
        end
    end
    return resize!(v, wi - 1)
end

"""
    serial_mem(::Type{M}, input::DatastoreBuffer{<:ReadDatastore}, count_mode::CountMode, range::AbstractRange = 1:length(input)) where {M<:AbstractMer}

MerCounting's simplest kmer counting method.

Build a sorted list (vector) of kmer counts (MerCount), serially and entirely in memory.

!!! warning
    This function is a serial and in memory `MerCount` list builder that can build a
    kmer count from a ReadsDatastore on its own (if you have memory and time),
    but it is also intended to be composed into other multi-process or multi-threaded
    kmer counting strategies.

    This method estimates roughly how many kmers will be generated by the reads
    specified by `range` in the dataset. It then pre-allocates an array to contain
    them. It then collects the kmers, sorts, them, and then collapses them into a
    list of counts sorted by the kmer.
    
    So if you want to count kmers and have the resources to throw at it, this is
    the simplest method, and possibly even the quickest given that simplicity.
"""
function serial_mem(::Type{M}, input::DatastoreBuffer{<:ReadDatastore}, count_mode::CountMode, range::AbstractRange = 1:length(input)) where {M<:AbstractMer}
    @info "Collecting kmers from all reads in $(name(ReadDatastores.datastore(input)))"
    all_mers = collect_mers(M, count_mode, input, range)
    return collapse_into_counts(all_mers)
end

macro serial_mem(K::Int, count_mode::Symbol, filename::String, range::Expr)
    # Figure out type of read datastore
    dstype = ReadDatastores.deduce_datastore_type(filename)
    # Process the range argument
    @assert range.head === :call
    @assert range.args[1] === :(:)
    interval = range.args[2]:range.args[3]
    quote
        open($dstype, $filename) do ds
            serial_mem(DNAMer{$K}, $count_mode, ds, $interval)
        end
    end
end