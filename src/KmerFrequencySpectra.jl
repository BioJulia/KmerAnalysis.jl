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
    
    i = firstindex(xfreqs)
    j = firstindex(yfreqs)
    i_end = lastindex(xfreqs) + 1
    j_end = lastindex(yfreqs) + 1
    
    while i < i_end && j < j_end
        fi = xfreqs[i]
        fj = yfreqs[j]
        if mer(fi) === mer(fj)
            hist[freq(fi) + 1, freq(fj) + 1] += 1
            i += 1
            j += 1
        elseif mer(fi) < mer(fj)
            hist[freq(fi) + 1, 1] += 1
            i += 1
        elseif mer(fj) < mer(fi)
            hist[1, freq(fj) + 1] += 1
            j += 1
        end
    end
    
    while i < i_end
        fi = xfreqs[i]
        hist[freq(fi) + 1, 1] += 1
        i += 1
    end
    
    while j < j_end
        fj = yfreqs[j]
        hist[1, freq(fj) + 1] += 1
        j += 1
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

function _find_plottable_subset(spec::KmerFrequencySpectra{1})
    return view(spec.data, firstindex(spec.data):findprev(x -> x > 0, spec.data, lastindex(spec.data)))
end

function _find_plottable_subset(spec::KmerFrequencySpectra{2})
    mat = spec.data
    furthest_row = 0
    furthest_col = 0
    for i in 1:256
        for j in 1:256
            v = mat[i, j]
            furthest_row = ifelse((j > furthest_row) & !iszero(v), j, furthest_row)
            furthest_col = ifelse((i > furthest_col) & !iszero(v), i, furthest_col)
        end
    end
    return view(mat, 1:furthest_col, 1:furthest_row)
end

function AbstractPlotting.convert_arguments(::AbstractPlotting.PointBased, spec::KmerFrequencySpectra{1})
    tv = _find_plottable_subset(spec)
    return ([Point2f0(i, j) for (i, j) in enumerate(tv)],)
end

function AbstractPlotting.convert_arguments(p::AbstractPlotting.SurfaceLike, spec::KmerFrequencySpectra{2})
    return AbstractPlotting.convert_arguments(p, _find_plottable_subset(spec))
end

AbstractPlotting.plottype(::KmerFrequencySpectra{1}) = AbstractPlotting.BarPlot
AbstractPlotting.plottype(::KmerFrequencySpectra{2}) = AbstractPlotting.Heatmap