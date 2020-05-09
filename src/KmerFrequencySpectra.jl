"""
A type for storing a frequency histogram of MerCounts, also sometimes referred
to as a Kmer spectra. 
"""
struct KmerFrequencySpectra{N}
    data::Array{UInt64, N}
    min::UInt8
end

function Base.summary(io::IO, hist::KmerFrequencySpectra)
    print(io, "Count histogram of motifs appearing more than ", hist.min, " times")
end

Base.show(io::IO, hist::KmerFrequencySpectra) = summary(io, hist)

Base.firstindex(x::KmerFrequencySpectra) = firstindex(x.data)
Base.lastindex(x::KmerFrequencySpectra) = lastindex(x.data)
Base.eachindex(x::KmerFrequencySpectra) = Base.OneTo(lastindex(x))

@inline function boundscheck(x::KmerFrequencySpectra, i)
    if i < 1 || i > lastindex(x.data)
        throw(BoundsError(x, i))
    end
    return true
end

@inline function Base.getindex(x::KmerFrequencySpectra, i)
    @boundscheck checkbounds(x, i)
    return @inbounds x.data[i]
end

"""
    spectra(freqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}

Build a kmer frequency spectra, excluding any kmer counts that don't meet
`min_count`.
"""
function spectra(freqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}
    hist = zeros(UInt64, 256)
    for x in freqs
        f = freq(x)
        if f ≥ min_count
            old = hist[f]
            hist[f] = old + 1
        end
    end
    return KmerFrequencySpectra{1}(hist, convert(UInt8, min_count))
end

function spectra(xfreqs::Vector{MerCount{M}}, yfreqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}
    hist = zeros(UInt64, 256, 256)
    
    i = j = 1
    nx = length(xfreqs)
    ny = length(yfreqs)
    fi = xfreqs[i]
    fj = yfreqs[j]
    
    while i <= nx && j <= ny
        if mer(fi) === mer(fj)
            hist[freq(fi) + 1, freq(fj) + 1] += 1
            fi = xfreqs[i += 1]
            fj = yfreqs[j += 1]
        elseif mer(fi) < mer(fj)
            hist[freq(fi) + 1, 1] += 1
            fi = xfreqs[i += 1]
        elseif mer(fj) < mer(fi)
            hist[1, freq(fj) + 1] += 1
            fj = yfreqs[j += 1]
        end
    end
    
    while i <= nx
        hist[freq(fi) + 1, 1] += 1
        fi = xfreqs[i += 1]
    end
    
    while j <= ny
        hist[1, freq(fj) + 1] += 1
        fj = yfreqs[j += 1]
    end
    
    return KmerFrequencySpectra{2}(hist, convert(UInt8, min_count))
end

"""
    spectra!(freqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}

A version of `spectra` that builds a kmer frequency spectra, AND mutates the
input `freqs` by filtering it to remove MerCounts that did not meet `min_count`
and so were not included in the constructed MerCountHist.
"""
function spectra!(freqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}
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
    return KmerFrequencySpectra{1}(hist, convert(UInt8, min_count))
end

function AbstractPlotting.convert_arguments(::AbstractPlotting.PointBased, spec::KmerFrequencySpectra{1})
    return ([Point2f0(i, j) for (i, j) in enumerate(spec.data)],)
end

AbstractPlotting.plottype(::KmerFrequencySpectra{1}) = AbstractPlotting.BarPlot
