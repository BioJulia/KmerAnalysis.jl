struct IndexedCounts{M<:AbstractMer,C<:CountMode}
    mer_index::Vector{M}            # Ordered list of kmer that contain counts
    count_names::Vector{Symbol}     # Names of the count vectors
    counts::Vector{Vector{UInt16}}  # Count vectors, each contains an entry per mer in the index
    name::Symbol                    # The name of this MerCounts
end

# Construction of empty MerCounts with a name
@inline function IndexedCounts{M,C}(name::Symbol) where {M<:AbstractMer,C<:CountMode}
    return IndexedCounts{M,C}(Vector{M}(), Vector{Symbol}(), Vector{Vector{UInt16}}(), name)
end

@inline function IndexedCounts(::Type{M}, name::Symbol, mode::C = CANONICAL) where {M<:AbstractMer,C<:CountMode}
    return IndexedCounts{M,C}(name)
end

# Construction from a sequence distance graph
function IndexedCounts(name::Symbol, sequences, mode::C = CANONICAL) where {M<:AbstractMer,C<:CountMode}
    @info string("Creating an empty kmer count datastore called ", name)
    mc = IndexedCounts{M,C}(name, mode)
    index!(mc, sequences)
    return mc
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

function index!(mc::IndexedCounts{M,C}, ref_sequences::Vector{<:LongSequence}) where {M<:AbstractMer,C<:CountMode}
    @info "Indexing kmer counts of reference sequences"
    @info "Populating index with reference sequence mers"
    # Add all Kmers from reference sequences.
    # t = expected number of kmers.
    t = sum(length(x) for x in ref_sequences) - BioSequences.ksize(M) + 1
    mcindex = mc.mer_index
    resize!(mcindex, t)
    vec_i = 1
    f = C()
    for seq in ref_sequences
        for mer in each(M, seq)
            mcindex[vec_i] = f(mer)
            vec_i = vec_i + 1
        end
    end
    resize!(mcindex, vec_i)
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