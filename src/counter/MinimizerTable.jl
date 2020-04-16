struct MinimizerTable{K,M}
    nbins::UInt32 # Number of parts onto which the kmer space will be subdivided.
    minimizer_table::Vector{UInt32} # The number of times each minimizer has been seen.
    bin_map::Vector{UInt32} # Indexing a minimizer returns the bin it belongs to.
end

function MinimizerTable{K,M}(nbins::Integer) where {K,M}
    @assert K > M
    return MinimizerTable{K,M}(
        convert(UInt32, nbins),
        zeros(UInt32, (1 << (2 * M)) + 1),
        Vector{UInt32}()
    )
end

Base.copy(mt::M) where {M<:MinimizerTable} = M(mt.nbins, copy(mt.minimizer_table), copy(mt.bin_map))

@inline function merge_tables(vec::Vector{M}) where {M<:MinimizerTable}
    mt = copy(vec[1])
    @inbounds for i in 2:lastindex(vec)
        tbl = vec[i]
        mt.minimizer_table .+= tbl.minimizer_table
    end
    return mt
end

function minimizer(::Type{DNAMer{J}}, mer::DNAMer{K}) where {J,K}
    mask = (one(UInt64) << (2 * J)) - 1
    mer_data = BioSequences.encoded_data(mer)
    min_bits = mer_data & mask
    for i in (K - J + 1):-1:1
        new_bits = mer_data & mask
        min_bits = ifelse(new_bits < min_bits, new_bits, min_bits)
        mer_data = mer_data >> 2
    end
    return DNAMer{J}(min_bits)
end

@inline function count_minimizers!(mt::MinimizerTable{K,M}, count_mode::CountMode, seq::LongSequence{A}) where {K,M,A<:DNAAlphabet}
    tbl = mt.minimizer_table
    @inbounds for elem in each(DNAMer{K}, seq)
        mer = count_mode(elem)
        mmer = BioSequences.encoded_data(minimizer(DNAMer{M}, mer))
        tbl[mmer] = tbl[mmer] + 1
    end
end

function initialize_minimizer_table!(mt::MinimizerTable, count_mode::CountMode, datastore::DatastoreBuffer{<:ReadDatastore}, part::AbstractRange)
    seq = eltype(datastore)()
    for i in part
        load_sequence!(datastore, i, seq)
        count_minimizers!(mt, count_mode, seq)
    end
    return mt
end

function calculate_bins!(mt::MinimizerTable)
    minicounts = mt.minimizer_table
    binmap = mt.bin_map
    resize!(binmap, length(minicounts))
    sorted_minimers = sortperm(mt.minimizer_table, rev = true)
    binsums = zeros(Int, mt.nbins)
    for minimer in sorted_minimers
        selected_bin = argmin(binsums)
        minimer_count = minicounts[minimer]
        binsums[selected_bin] += minimer_count
        binmap[minimer] = selected_bin
    end
    return mt
end

function minimizer_table_maker(::Type{M}, ::Type{D}, filename::String, count_mode::CountMode, nbins::Int) where {M<:MinimizerTable,D<:ReadDatastore}
    mt = M(nbins)
    ds = buffer(open(D, filename))
    part = myid():nprocs():length(ds)
    initialize_minimizer_table!(mt, count_mode, ds, part)
    return mt
end

function calculate_minimizer_distribution(::Type{M}, ::Type{D}, filename::String, count_mode::CountMode, nbins::Int) where {M<:MinimizerTable,D<:ReadDatastore}
    local_tables = Vector{M}(undef, nprocs())
    @sync for p in procs()
        @async local_tables[p] = remotecall_fetch(minimizer_table_maker, p, M, D, filename, count_mode, nbins)
    end
    global_table = merge_tables(local_tables)
    return calculate_bins!(global_table)
end
