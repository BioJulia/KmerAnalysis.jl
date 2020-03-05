struct IndexedCounts{M<:AbstractMer}
    mer_index::Vector{M}            # Ordered list of kmer that contain counts
    count_names::Vector{Symbol}     # Names of the count vectors
    counts::Vector{Vector{UInt16}}  # Count vectors, each contains an entry per mer in the index
    name::Symbol                    # The name of this MerCounts
    counting_mode::CountMode
end

# Construction of empty MerCounts with a name
@inline function IndexedCounts{M}(name::Symbol, mode::CountMode = Canonical) where {M<:AbstractMer}
    return IndexedCounts{M}(Vector{M}(), Vector{Symbol}(), Vector{Vector{UInt16}}(), name, mode)
end

# Construction from a sequence distance graph
function IndexedCounts{M}(name::Symbol, sequences, mode::CountMode = Canonical) where {M<:AbstractMer}
    @info string("Creating an empty kmer count datastore called ", name)
    mc = IndexedCounts{M}(name, mode)
    index!(mc, sequences)
    return mc
end

function _determine_kindex_size(::Type{M}, sequences) where {M<:AbstractMer}
    # TODO: This can probably be done in parallel very simply with the distributed
    # for loop macro.
    t = 0
    for seq in sequences
        if length(node) >= BioSequences.ksize(M)
            t = t + length(node) + 1 - BioSequences.ksize(M)
        end
    end
    return t
end

function _count_and_collapse!(mers::Vector{M}, counts::Vector{UInt16}) where {M<:AbstractMer}
    wi = 1
    ri = 1
    stop = lastindex(mers) + 1
    @inbounds while ri < stop
        mers[wi] = mers[ri]
        ci = one(UInt16)
        while (ri += 1) < stop && mers[ri] == mers[wi]
            ci = ci + one(UInt16)
        end
        push!(counts, ci)
        wi = wi + 1
    end
    resize!(mers, ri)
    return nothing
end

function index!(mc::IndexedCounts{M}, sequences) where {M<:AbstractMer}
    @info "Indexing kmer counts of reference sequences"
    @info "Populating index with reference sequence mers"
    # Add all Kmers from reference sequences.
    t = _determine_kindex_size(M, sequences)
    mcindex = mc.mer_index
    resize!(mcindex, t)
    index_i = 1
    if mc.counting_mode === Canonical
        for seq in sequences
            if length(seq) >= BioSequences.ksize(M)
                for mer in each(M, seq)
                    mcindex[index_i] = canonical(mer)
                    index_i = index_i + 1
                end
            end
        end
    else
        for seq in sequences
            if length(seq) >= BioSequences.ksize(M)
                for mer in each(M, seq)
                    mcindex[index_i] = fwmer(mer)
                    index_i = index_i + 1
                end
            end
        end
    end
    @info "Sorting the index"
    # Sort the index
    sort!(mcindex)
    @info "Counting and collapsing mers in index"
    # Create a first count;
    # Collapse the kmer index, saving the coverage to the first count.
    first_count = Vector{UInt16}()
    sizehint!(first_count, length(mcindex))
    push!(mc.counts, first_count)
    push!(mc.count_names, :sdg)
    _count_and_collapse!(mcindex, first_count)
    return mc
end

function Base.summary(io::IO, mc::IndexedCounts{M}) where {M}
    print(io, BioSequences.ksize(M), "-mer Indexed Counts Datastore '",  mc.name, "': ", length(mc.counts), " stored ", BioSequences.ksize(M), "-mer counts")
end

function Base.show(io::IO, mc::IndexedCounts{M}) where {M}
    summary(io, mc)
end