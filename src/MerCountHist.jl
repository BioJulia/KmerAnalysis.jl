"""
A type for storing a frequency histogram of MerCounts, also sometimes referred
to as a Kmer spectra. 
"""
struct MerCountHist
    data::Vector{UInt64}
    min::UInt8
end

function Base.summary(io::IO, hist::MerCountHist)
    print(io, "Count histogram of motifs appearing more than ", hist.min, "times")
end

"""
    hist!(freqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}

Build a histogram of kmer frequencies, excluding any kmer counts that don't meet
`min_count`.
"""
function hist(freqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}
    hist = zeros(UInt64, 256)
    for x in freqs
        f = freq(x)
        if f ≥ min_count
            old = hist[f]
            hist[f] = old + 1
        end
    end
    return MerFreqHist(hist, convert(UInt8, min_count))
end

"""
    hist!(freqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}

A version of `hist` that builds a MerCountHist, AND mutates the input `freqs` by
filtering it to remove MerCounts that did not meet `min_count` and so were not
included in the constructed MerCountHist.
"""
function hist!(freqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}
    hist = zeros(UInt64, 256)
    wi = firstindex(freqs)
    used = 0
    for x in freqs
        f = freq(x)
        if f ≥ min_count
            old = hist[f]
            hist[f] = old + 1
            freqs[wi] = x
            wi = wi + 1
            used = used + 1
        end
        resize!(freqs, used)
    end
    return MerCountHist(hist, convert(UInt8, min_count))
end