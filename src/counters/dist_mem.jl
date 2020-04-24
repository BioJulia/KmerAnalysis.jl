
function _do_serial_mem(::Type{M}, ::Type{D}, filename::String, count_mode::CountMode, i) where {M<:AbstractMer,D<:ReadDatastore}
    ds = open(D, filename)
    part = i:nworkers():length(ds)
    counts = serial_mem(M, buffer(ds), count_mode, part)
    close(ds)
    return counts
end

"""
The multi-process parallel version of the [`serial_mem`](@ref) kmer counter.

Like [`serial_mem`](@ref) this counter does all counting and processing in memory,
so it is not designed to limit memory use, but if using a cluster of machines the
memory used per worker machine would be less. Of course the data is copied back
over to the master process and the results combined and so the master node has
to have enough RAM available to hold all the result.
"""
function dist_mem end

function dist_mem(::Type{M}, ::Type{D}, filename::String, count_mode::CountMode) where {M<:AbstractMer,D<:ReadDatastore}
    nw = nworkers()
    @info "Splitting all reads in $filename, across $nw worker processes and counting kmers"
    if nw === 1
        @warn "There are not any extra worker processes available. This method will not provide any speed benefit over serial_mem."
    end
    local_counts = Vector{Vector{MerCount{M}}}(undef, nw)
    @sync begin
        @inbounds for (i, w) in enumerate(workers())
            @async local_counts[i] = remotecall_fetch(_do_serial_mem, w, M, D, filename, count_mode, i)
        end
        @info "Waiting for workers to finish"
    end
    @info "Merging all counts from $nw worker processes"
    final_count = local_counts[1]
    @inbounds for i in 2:lastindex(local_counts)
        unsafe_merge_into!(final_count, local_counts[i])
    end
    return final_count
end

function dist_mem(::Type{M}, input::ReadDatastore, count_mode::CountMode) where {M<:AbstractMer}
    return dist_mem(M, typeof(input), input.filename, count_mode)
end

function dist_mem(::Type{M}, input::DatastoreBuffer{<:ReadDatastore}, count_mode::CountMode) where {M<:AbstractMer}
    return dist_mem(M, ReadDatastores.datastore(input), count_mode)
end