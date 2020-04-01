struct MinimizerTable{A<:NucleicAcidAlphabet{2},K,M}
    nbins::UInt32 # Number of parts onto which the kmer space will be subdivided.
    minimizer_table::Vector{UInt32} # The number of times each minimizer has been seen.
    bin_map::Vector{UInt32} # Indexing a minimizer returns the bin it belongs to.
end

function MinimizerTable{A,K,M}(nbins::Integer) where {A<:NucleicAcidAlphabet{2},K,M}
    @assert K > M
    return MinimizerTable{A,K,M}(
        convert(UInt32, nbins),
        zeros(UInt32, (1 << (2 * M)) + 1),
        Vector{UInt32}()
    )
end

function MinimizerTable(::Type{K}, ::Type{M}, nbins::Integer) where {K<:AbstractMer,M<:AbstractMer}
    @assert typeof(Alphabet(K)) === typeof(Alphabet(M))
    return MinimizerTable{typeof(Alphabet(K)),BioSequences.ksize(K),BioSequences.ksize(M)}(nbins)
end

@inline function Base.setindex(mt::MinimizerTable, i, val::Integer)
    return mt.minimizer_table[i] = val
end

function minimizer(::Type{DNAMer{J}}, mer::DNAMer{K}) where {J,K}
    mask = (one(UInt64) << (2 * J)) - 1
    mer_data = BioSequences.encoded_data(mer)
    min_bits = mer_data & mask
    #min_i = 0
    for i in (K - J + 1):-1:1
        new_bits = mer_data & mask
        #min_i = ifelse(new_bits < min_bits, i, min_i)
        min_bits = ifelse(new_bits < min_bits, new_bits, min_bits)
        mer_data = mer_data >> 2
    end
    return DNAMer{J}(min_bits)
end

function count_minimizers!(mt::MinimizerTable, seq::LongSequence{DNAAlphabet{2}})
    @inbounds for mer in each(DNAMer{31}, seq)
        mmer = encoded_data(minimizer(DNAMer{7}, mer))
        mt[mmer] = mt[mmer] + 1
    end
end

# Unspecified type for sequences to allow compiler to adapt to any iterable that
# yields a LongSequence{DNAAlphabet{2}} each iteration.
function initialize_minimizer_table!(mt::MinimizerTable, sequences)
    for seq in sequences
        count_minimizers!(mt, seq)
    end
    return mt
end