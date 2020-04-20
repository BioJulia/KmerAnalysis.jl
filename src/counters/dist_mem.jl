#function do_serial_mem(::Type{M}, ::Type{D}, filename::String, count_mode::CountMode,)
#    ds = buffer(open(D, filename))
#    part = myid():nprocs():length(ds)
#    serial_mem(M, count_mode, ds, part)
#end

function dist_mem(::Type{M}, input::DatastoreBuffer{<:ReadDatastore}, count_mode::CountMode) where {M<:AbstractMer}
    @info "Splitting all reads in $(name(ReadDatastores.datastore(input))), across $(nprocs()) processes and counting kmers"
    local_counts = Vector{Vector{MerCount{M}}}(undef, nprocs())
    @sync for p in procs()
        p:nprocs():length(input)
        @async local_counts[p] = remotecall_fetch(serial_mem, p, M, input, count_mode, p:nprocs():length(input))
    end
    final_count = local_counts[1]
    @inbounds for i in 2:lastindex(local_counts)
        unsafe_merge_into!(final_count, local_counts[i])
    end
    return final_count
end