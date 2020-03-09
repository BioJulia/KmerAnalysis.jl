module MerCounting

export
    CountMode,
    IndexedCounts

using
    BioSequences,
    ReadDatastores
    


@enum CountMode Canonical NonCanonical

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

"""
    unsafe_collapse_into_counts!(result::Vector{MerCount{M}}, mers::Vector{M}) where {M<:AbstractMer}

!!! warning
    This method is marked as unsafe because it assumes that the `mers` input
    vector is already sorted.
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
        push!(result, MerCount{M}(mers[wi], ci))
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
"""
function collapse_into_counts(mers::Vector{M}) where {M<:AbstractMer}
    return collapse_into_counts!(Vector{MerCount{M}}(), mers)
end


include("IndexedCounts.jl")

end # module
