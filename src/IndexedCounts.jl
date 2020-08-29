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
    resize!(mers, wi - 1)
    return nothing
end

function index!(mc::IndexedCounts{M,C}, ref_sequences) where {M<:AbstractMer,C<:CountMode}
    @assert eltype(ref_sequences) <: LongSequence
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

"""
    add_count!(mc::IndexedCounts{M}, name::Symbol, reads::ReadDatastore) where {M}

Count the k-mers in a ReadDatastore and add them to the IndexedCounts{M}.

Currently IndexedCounts use a serial and in-memory counting method, as 

!!! note
    Only the k-mers in the IndexedCounts' index shall be counted from the reads. 
"""
function add_count!(mc::IndexedCounts{M,C}, name::Symbol, reads::ReadDatastore) where {M,C}
    mode = C()
    # TODO: Check if name already exists.
    push!(mc.count_names, name)
    c = zeros(UInt16, length(mc.mer_index))
    push!(mc.counts, c)
    # Populate lookup map.
    kmer_map = Dict{M, Int}()
    sizehint!(kmer_map, length(mc.mer_index))
    for (i, j) in enumerate(mc.mer_index)
        kmer_map[j] = i
    end
    for read in reads
        for mer in each(M, read)
            findk = get(kmer_map, mode(mer), 0)
            @inbounds if findk > 0
                old = c[findk]
                inc = old + one(UInt16)
                c[findk] = ifelse(old > inc, typemax(UInt16), inc)
            end
        end
    end
end

function project_count(mc::IndexedCounts{M,C}, name::Symbol, sequence::BioSequence) where {M,C}
    countidx = findfirst(isequal(name), mc.count_names)
    if isnothing(countidx)
        throw(ArgumentError(string("No count dataset called", string(name), "exists")))
    end
    counts = mc.counts[countidx]
    querymers = collect_mers(M, C(), sequence)
    projection = zeros(UInt16, length(querymers))
    indexsize = length(mc.mer_index)
    index = mc.mer_index
    for (i, mer) in enumerate(querymers)
        kidx = min(searchsortedfirst(index, mer), indexsize)
        if index[kidx] === mer
            projection[i] = counts[kidx]
        end
    end
    return projection
end