# SerialMemCounter
# ================

"""
    SerialMemCounter{M<:AbstractMer,C<:CountMode}

KmerAnalysis' simplest kmer counter.

Build a sorted list (vector) of kmer counts (MerCount), serially and entirely in memory.

Build a SerialMemCounter using the [`serial_mem`](@ref) method. This returns
a SerialMemCounter.

A SerialMemCounter can be treated as a function or functor and passed to other
functions as an argument, and it can be called on other arguments.

It is a general purpose kmer counting method:

It basically collects all the kmers from iterating over the input, before sorting
the collected kmers and collapsing them into counts. To reduce allocations and 
GC from many repeated calls of a SerialMemCounter, they keep a few buffers that
are used often for kmer collection and count collapsing. 

Whilst general purpose for a variety of argument types, it is most suitable for:

1. Counting the kmers in a small number of long-ish `BioSequences` (e.g. a
   reference genome or low(ish)-complexity genome graph.).

2. Counting the kmers in a ReadDatastore if the dataset size is small enough
   such that all the kmers collected from every dataset can fit in a single
   machines memory. If you have a decent MacBook Pro, for example, datasets like
   _E. coli_ paired end sequencing reads will be no problem.

3. Other larger ReadDatastores IF you have the RAM to throw at it, e.g. you're on
   a HPC machine configured for large memory jobs.
"""
struct SerialMemCounter{M<:AbstractMer,C<:CountMode} <: AbstractKmerCounter{M}
    mer_buffer::Vector{M}
    count_buffer::Vector{MerCount{M}}
end

serial_mem(::Type{M}, count_mode::C) where {M<:AbstractMer,C<:CountMode} = SerialMemCounter{M,C}(Vector{M}(), Vector{MerCount{M}}())

@inline function count!(smc::SerialMemCounter{M,C}, input, args...) where {M<:AbstractMer,C<:CountMode}
    _serial_mem!(smc.count_buffer, smc.mer_buffer, C(), input, args...)
end

@inline function (smc::SerialMemCounter{M,C})(input, args...) where {M<:AbstractMer,C}
    count!(smc, input, args...)
    return copy(smc.count_buffer)
end

@inline function _serial_mem!(count_buffer::Vector{MerCount{M}}, mer_buffer::Vector{M}, count_mode::CountMode, input, args...) where {M<:AbstractMer}
    collect_mers!(mer_buffer, count_mode, input, args...)
    collapse_into_counts!(count_buffer, mer_buffer)
end

@inline function serial_mem(::Type{M}, input::PairedReads{A}, count_mode::CountMode) where {M,A}
    counts = Vector{MerCount{M}}()
    mers = Vector{M}()
    _serial_mem!(counts, mers, count_mode, input)
    return counts
end