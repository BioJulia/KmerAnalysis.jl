"""
A simple mer count struct.

MerCount is a simple struct that binds a mer value to a count of the number
of times it has been observed.

This type, (sorted) vectors of them, and some additional utility methods, can
form the basic building blocks of the higher-level mer counting functionality.
This struct can also be an eltype of kmer hashes or more involved specialized
types that store counts of Kmers.

!!! note
    The count is stored as an UInt8 because often once the count is more than
    255 we hardly care anymore.
"""
struct MerCount{M<:AbstractMer}
    mer::M
    count::UInt8
    function MerCount{M}(mer::M, count::Integer) where {M<:AbstractMer}
        return new(mer, convert(UInt8, min(typemax(UInt8), count)))
    end
end

"Shorthand for `MerCount{DNAMer{K}}`"
const DNAMerCount{K} = MerCount{DNAMer{K}}

"Shorthand for `MerCount{RNAMer{K}}`"
const RNAMerCount{K} = MerCount{RNAMer{K}}

"Get the mer from a `MerCount`."
@inline mer(x::MerCount{<:AbstractMer}) = x.mer

"Get the count from a `MerCount`."
@inline freq(x::MerCount{<:AbstractMer}) = x.count

"Get the count from a `MerCount`, and convert it to type R."
@inline freq(::Type{R}, x::MerCount{<:AbstractMer}) where {R<:Real} = convert(R, freq(x))

"""
    Base.isless(x::MerCount{M}, y::MerCount{M}) where {M<:AbstractMer}

Check if a MerCount x is less than MerCount y.

!!! note
    The comparison only considers the mer, not the count.
"""
Base.isless(x::MerCount{M}, y::MerCount{M}) where {M<:AbstractMer} = mer(x) < mer(y)

"""
    Base.:(>)(x::MerCount{M}, y::MerCount{M}) where {M<:AbstractMer}

Check if a MerCount x is greater than MerCount y.

!!! note
    The comparison only considers the mer, not the count.
"""
Base.:(>)(x::MerCount{M}, y::MerCount{M}) where {M<:AbstractMer} = mer(x) > mer(y)

function merge(x::MerCount{M}, y::MerCount{M}) where {M<:AbstractMer}
    return MerCount{M}(x.mer, UInt16(x.count) + UInt16(y.count))
end

function Base.show(io::IO, mfreq::MerCount{<:AbstractMer})
    print(io, mer(mfreq), " occurs ", freq(mfreq), " times")
end

function collect_mers(count_mode::CountMode, input::BioSequences.EveryMerIterator)
    return collect(count_mode(x) for x in input)
end

function collect_mers(::Type{M}, count_mode::CountMode, input::BioSequence) where {M<:AbstractMer}
    return collect_mers(count_mode, each(M, input))
end

function collect_mers(::Type{M}, count_mode::CountMode, input::DatastoreBuffer{<:PairedReads}, range::AbstractRange = 1:length(input)) where {M<:AbstractMer}
    max_read_size = max_read_length(ReadDatastores.datastore(input))
    v = Vector{M}(undef, length(input) * (max_read_size - BioSequences.ksize(M) + 1))
    wi = firstindex(v)
    read_sequence = eltype(input)()
    @inbounds for i in range
        for mer in each(M, load_sequence!(input, i, read_sequence))
            v[wi] = count_mode(mer)
            wi = wi + 1
        end
    end
    return resize!(v, wi - 1)
end

collect_mers(v::Vector{M}) where {M<:AbstractMer} = v

"""
    unsafe_collapse_into_counts!(result::Vector{MerCount{M}}, mers::Vector{M}) where {M<:AbstractMer}

!!! warning
    This method is internal and marked as unsafe because it assumes that the
    `mers` input vector is already sorted.
"""
function unsafe_collapse_into_counts!(result::Vector{MerCount{M}}, mers::Vector{M}) where {M<:AbstractMer}
    wi = 1
    ri = 1
    stop = lastindex(mers) + 1
    empty!(result)
    @inbounds while ri < stop
        ci = one(UInt16)
        while (ri += 1) < stop && mers[ri] == mers[wi]
            ci = ci + one(UInt16)
        end
        push!(result, MerCount{M}(mers[wi], ci)) # TODO: See about removing this push!.
        wi = ri
    end
    return result
end

"""
    collapse_into_counts!(result::Vector{MerCount{M}}, mers::Vector{M}) where {M<:AbstractMer}

Build a vector of sorted `MerCount`s from a Vector of a mer type.

This is a basic kernel function used for any higher level and more complex
kmer counting procedures.

This is like `collapse_into_counts`, except it's first argument is a `result`
vector that is cleared and filled with the result.

!!! note
    The input vector `mers` will be sorted by this method.
"""
function collapse_into_counts!(result::Vector{MerCount{M}}, mers::Vector{M}) where {M<:AbstractMer}
    sort!(mers)
    return unsafe_collapse_into_counts!(result, mers)
end

"""
    collapse_into_counts(mers::Vector{M}) where {M<:AbstractMer}

Build a vector of sorted `MerCount`s from a Vector of a mer type.

This is a basic kernel function used for any higher level and more complex
kmer counting procedures.

!!! note
    The input vector `mers` will be sorted by this method.
"""
function collapse_into_counts!(mers::Vector{M}) where {M<:AbstractMer}
    return collapse_into_counts!(Vector{MerCount{M}}(), mers)
end

# The dispatch here could be a little better. - Distinguising between Mer iterators and
# non mer iterators perhaps.
function Vector{MerCount{M}}(mode::CountMode, input::BioSequences.EveryMerIterator{M}) where {M<:AbstractMer}
    mers = collect_mers(mode, input)
    return collapse_into_counts!(mers)
end

function Vector{MerCount{M}}(mode::CountMode, input::Any, args...) where {M<:AbstractMer}
    mers = collect_mers(M, mode, input, args...)
    return collapse_into_counts!(mers)
end

"""
    unsafe_merge_into!(a::Vector{MerCount{M}}, b::Vector{MerCount{M}}) where {M<:AbstractMer}

Merge the `MerCount`s from vector `b` into the vector `a`.

!!! warning
    This method is marked as unsafe as it assumes both of the input vectors `a`
    and `b` are already sorted.
"""
function unsafe_merge_into!(a::Vector{MerCount{M}}, b::Vector{MerCount{M}}) where {M<:AbstractMer}
    a_i = firstindex(a)
    a_end = lastindex(a) + 1
    b_i = b_i2 = firstindex(b)
    b_end = lastindex(b) + 1
    
    # Merge, accumulating counts on `a`.
    @inbounds while b_i < b_end
        # Move cursor a_i forward so long as mers in a are less than the mer
        # in b[b_i]
        while a_i < a_end && a[a_i] < b[b_i]
            a_i = a_i + 1
        end
        # Check if the mers at the a_i and b_i cursors are the same. If so,
        # the count in b needs to be merged into the count in a.
        if a_i < a_end && mer(a[a_i]) == mer(b[b_i])
            # Combine entries
            a[a_i] = merge(a[a_i], b[b_i])
            a_i = a_i + 1
            b_i = b_i + 1
        end
        while b_i < b_end && (a_i == a_end || b[b_i] < a[a_i])
            b[b_i2] = b[b_i]
            b_i2 = b_i2 + 1
            b_i = b_i + 1
        end
    end
    # Shrink `b` to the size of the remaining contents.
    resize!(b, b_i2 - 1)
    
    # Expand `a` to allow the insertion of unique values in `b`.
    oldsize = length(a)
    resize!(a, oldsize + length(b))
    r_a = oldsize
    
    # Merge-sort from the bottom into `a`.
    wr_a = lastindex(a)
    rend_a = firstindex(a)
    r_b = lastindex(b)
    r_end_b = firstindex(b)
    @inbounds while wr_a >= rend_a
        if r_b >= r_end_b && (r_a < rend_a || b[r_b] > a[r_a])
            a[wr_a] = b[r_b]
            r_b = r_b - 1
        else
            a[wr_a] = a[r_a]
            r_a = r_a - 1
        end
        wr_a = wr_a - 1
    end
    empty!(b)
    return a
end

"""
    merge_into!(a::Vector{MerCount{M}}, b::Vector{MerCount{M}}) where {M<:AbstractMer}

Merge the `MerCount`s from vector `b` into the vector `a`.

!!! note
    This will sort the input vectors `a` and `b`.
"""
function merge_into!(a::Vector{MerCount{M}}, b::Vector{MerCount{M}}) where {M<:AbstractMer}
    sort!(a)
    sort!(b)
    return unsafe_merge_into!(a, b)
end

function unsafe_collapse!(freqs::Vector{MerCount{M}}) where {M<:AbstractMer}
    wi = 1
    ri = 1
    pi = 1
    stop = lastindex(freqs) + 1
    @inbounds while ri < stop
        ci = one(UInt16)
        while (ri += 1) < stop && mer(freqs[ri]) == mer(freqs[pi])
            ci = ci + one(UInt16)
        end
        freqs[wi] = MerCount{M}(mer(freqs[pi]), ci)
        pi = ri
        wi = wi + 1
    end
    resize!(freqs, wi - 1)
    return freqs
end

collapse!(freqs::Vector{MerCount{M}}) where {M<:AbstractMer} = collapse_sorted!(sort!(freqs))

