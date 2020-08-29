
function kci!(freqs::Vector{<:Integer}, kci_peak::Float64)
    sort!(freqs)
    return length(freqs) > 10 ? freqs[div(length(freqs), 2)] / kci_peak : -1.0
end

kci(freqs::Vector{<:Integer}, kci_peak::Float64) = kci!(copy(freqs), kci_peak)

function kci(indexed_counts::IndexedCounts{M,C}, name::Symbol, sequence::BioSequence, kci_peak::Float64) where {M,C}
    freqs = project_count(indexed_counts, name, sequence)
    filter!(!iszero, freqs)
    return kci!(freqs, kci_peak)
end