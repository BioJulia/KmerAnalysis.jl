"""
A type for storing a frequency histogram of MerCounts, also sometimes referred
to as a Kmer spectra. 
"""
struct KmerFrequencySpectra{N}
    data::Array{UInt64, N}
    min::UInt8
end

function KmerFrequencySpectra{1}(min::Integer = 0)
    return KmerFrequencySpectra{1}(zeros(UInt64, 257), convert(UInt8, min))
end

"""
    KmerFrequencySpectra{1}(freqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}

Build a 1 dimensional k-mer frequency spectra, from a vector of kmer counts,
excluding any k-mer counts that don't meet `min_count`.
"""
function KmerFrequencySpectra{1}(freqs::Vector{MerCount{M}}, min_count::Integer) where {M<:AbstractMer}
    spec = KmerFrequencySpectra{1}(min_count)
    sdat = spec.data
    for x in freqs
        f = freq(x)
        if f ≥ min_count
            i = f + one(UInt8)
            old = sdat[i]
            sdat[i] = old + one(UInt8)
        end
    end
    return spec
end

function KmerFrequencySpectra{1}(freqs::Vector{MerCount{M}}) where {M<:AbstractMer}
    spec = KmerFrequencySpectra{1}(0)
    sdat = spec.data
    for x in freqs
        i = freq(x) + one(UInt8)
        old = sdat[i]
        sdat[i] = old + one(UInt8)
    end
    return spec
end

function KmerFrequencySpectra{2}(min::Integer = 0)
    return KmerFrequencySpectra{2}(zeros(UInt64, 257, 257), convert(UInt8, min))
end

"""
    KmerFrequencySpectra{2}(xfreqs::Vector{MerCount{M}}, yfreqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}

Build a 2 dimensional k-mer frequency spectra, from two vectors of kmer counts,
excluding any k-mer counts that don't meet the `min_count`.

!!! tip
    Typically, `xfreqs` should be counts from a big dataset, typically something like
    raw reads, and `yfreqs` should be counts from a smaller or more collapsed sequence
    dataset like a reference genome or genome assembly graph.
    
    If you are comparring two similar sized datasets like two read datasets, then
    which one is passed as `xfreqs` and which is passed as `yfreqs` does not really matter.
"""
function KmerFrequencySpectra{2}(xfreqs::Vector{MerCount{M}}, yfreqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}
    spec = KmerFrequencySpectra{2}(min_count)
    sdat = spec.data
    
    i = firstindex(xfreqs)
    j = firstindex(yfreqs)
    i_end = lastindex(xfreqs) + 1
    j_end = lastindex(yfreqs) + 1
    
    while i < i_end && j < j_end
        fi = xfreqs[i]
        fj = yfreqs[j]
        if mer(fi) === mer(fj)
            sdat[freq(fi) + one(UInt8), freq(fj) + one(UInt8)] += 1
            i += 1
            j += 1
        elseif mer(fi) < mer(fj)
            sdat[freq(fi) + one(UInt8), one(UInt8)] += 1
            i += 1
        elseif mer(fj) < mer(fi)
            sdat[one(UInt8), freq(fj) + one(UInt8)] += 1
            j += 1
        end
    end
    
    while i < i_end
        fi = xfreqs[i]
        sdat[freq(fi) + one(UInt8), one(UInt8)] += 1
        i += 1
    end
    
    while j < j_end
        fj = yfreqs[j]
        sdat[one(UInt8), freq(fj) + one(UInt8)] += 1
        j += 1
    end
    
    return spec
end

function Base.summary(io::IO, hist::KmerFrequencySpectra)
    print(io, "Count histogram of motifs appearing more than ", hist.min, " times")
end

Base.show(io::IO, hist::KmerFrequencySpectra) = summary(io, hist)

Base.firstindex(x::KmerFrequencySpectra{1}) = 0
Base.lastindex(x::KmerFrequencySpectra{1}) = 256
Base.eachindex(x::KmerFrequencySpectra{1}) = firstindex(x):lastindex(x)

@inline function boundscheck(x::KmerFrequencySpectra{1}, i)
    if i < firstindex(x.data) || i > lastindex(x.data)
        throw(BoundsError(x, i))
    end
    return true
end

@inline function Base.getindex(x::KmerFrequencySpectra{1}, i)
    i′ = i + 1
    @boundscheck checkbounds(x, i′)
    return @inbounds x.data[i′]
end

@inline function Base.setindex(x::KmerFrequencySpectra{1}, i, val)
    i′ = i + 1
    @boundscheck checkbounds(x, i′)
    @inbounds x.data[i′] = val
end

"""
    spectra(freqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}

Build a 1 dimensional kmer frequency spectra, from a vector of kmer counts,
excluding any kmer counts that don't meet `min_count`.
"""
function spectra(freqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}
    return KmerFrequencySpectra{1}(freqs, min_count)
end

function spectra(c::C, input, min_count::Integer = 0) where {C<:AbstractKmerCounter}
    return spectra(c(input), min_count)
end

"""
    spectra(xfreqs::Vector{MerCount{M}}, yfreqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}

Build
"""
function spectra(xfreqs::Vector{MerCount{M}}, yfreqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}
    return KmerFrequencySpectra{2}(xfreqs, yfreqs, min_count)
end

"""
    spectra(x::Pair{I,C}, min_count::Integer = 0) where {I,C<:AbstractKmerCounter}

Build a 1 dimensional kmer frequency spectra, using a pair of (input data => kmer counter)
as the input argument `x`. The spectra function uses the kmer counter to get the counts
from the input data, before computing the spectra.

Any kmer counts that don't meet `min_count` are excluded.
"""
function spectra(x::Pair{I,C}, min_count::Integer = 0) where {I,C<:AbstractKmerCounter}
    input = first(x)
    counter = last(x)
    counts = counter(input)
    return spectra(counts, min_count)
end

function spectra(x::Pair{X,A}, y::Pair{Y,B}, min_count::Integer = 0) where {X,Y,A<:AbstractKmerCounter,B<:AbstractKmerCounter}
    xinput = first(x)
    xcounter = last(x)
    xcounts = xcounter(xinput)
    yinput = first(y)
    ycounter = last(y)
    ycounts = ycounter(yinput)
    
    return spectra(xcounts, ycounts, min_count)
end

mutable struct Peak
    born::Int
    died::Int
    left::Int
    right::Int
end

function build_persistent_homology(seq)
    peaks = Peak[]
    idx2peak = [-1 for s in seq]
    indicies = sortperm(seq, rev = true)
    for idx in indicies
        lftdone = idx > 1 && idx2peak[idx - 1] != -1
        rgtdone = idx < lastindex(seq) && idx2peak[idx + 1] != -1
        il = lftdone ? idx2peak[idx - 1] : -1 
        ir = rgtdone ? idx2peak[idx + 1] : -1
        
        # New peak born.
        if !lftdone && !rgtdone
            push!(peaks, Peak(idx, -1, idx, idx))
            idx2peak[idx] = length(peaks)
        end
        
        # Directly merge to next peak left.
        if lftdone && !rgtdone
            peaks[il].right += 1
            idx2peak[idx] = il
        end
        
        # Directly merge to next peak right.
        if !lftdone && rgtdone
            peaks[ir].left -= 1
            idx2peak[idx] = ir
        end
        
        # Merge left and right peaks
        if lftdone && rgtdone
            # Left was born earlier: merge right to left.
            if seq[peaks[il].born] > seq[peaks[ir].born]
                peaks[ir].died = idx
                peaks[il].right = peaks[ir].right
                idx2peak[peaks[il].right] = idx2peak[idx] = il
            else
                peaks[il].died = idx
                peaks[ir].left = peaks[il].left
                idx2peak[peaks[ir].left] = idx2peak[idx] = ir
            end
        end
    end
    return peaks
end

@inline persistence(p::Peak, seq) = p.died == -1 ? Inf : Float64(seq[p.born]) - Float64(seq[p.died])

"""
A simple struct describing the location and significance (or persistence)
of a peak detected in kmer frequency spectra data.

It has the following fields:

1. peak: The k-mer frequency at which the peak is at its highest.
2. left: The k-mer frequency constituting the left bound of the peak.
3. right: The k-mer frequency constituting the right bound of the peak.
4. persistence: The persistence of the peak in persistent homology analysis -
   essentially a significance value for the peak. Higher value = more significant.
"""
struct SpectraPeak
    peak::Int
    left::Int
    right::Int
    persistence::Float64
end

"""
    find_spectra_peaks(s::KmerFrequencySpectra{1})

Automatically detect all the peaks present in the signal of a 1D k-mer
frequency spectra, using a persistent homology method.

Returns a vector of [`SpectraPeak`](@ref).
"""
function find_spectra_peaks(s::KmerFrequencySpectra{1})
    ph = build_persistent_homology(s.data)
    peaks = [SpectraPeak(h.born - 1, h.left - 1, h.right - 1, persistence(h, s.data)) for h in ph]
    return peaks
end