#function do_serial_mem(::Type{M}, ::Type{D}, filename::String, count_mode::CountMode,)
#    ds = buffer(open(D, filename))
#    part = myid():nprocs():length(ds)
#    serial_mem(M, count_mode, ds, part)
#end

function dist_mem(::Type{M}, input::DatastoreBuffer{<:ReadDatastore}, count_mode::CountMode) where {M<:AbstractMer}
    @info "Splitting all reads in $(datastore(name(input))), across $(nprocs()) processes and counting kmers"
    local_counts = Vector{Vector{MerCount{M}}}(undef, nprocs())
    @sync for p in procs()
        p:nprocs():length(input)
        @async local_counts[p] = remotecall_fetch(serial_mem, M, input, count_mode, p:nprocs:length(input))
    end
    return local_counts
    #unsafe_merge_into!
    #all_mers = collect_mers(M, count_mode, input, range)
    #return collapse_into_counts(all_mers)
end