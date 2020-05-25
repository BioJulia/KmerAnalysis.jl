module TestKmerAnalysis

using KmerAnalysis, BioSequences, Test

@testset "Count Mode" begin
    iterres = BioSequences.MerIterResult(1, mer"TCGGAAGACTAAGGAGCCTATCAGATATGTC", mer"GACATATCTGATAGGCTCCTTAGTCTTCCGA")
    @test CANONICAL(iterres) === mer"GACATATCTGATAGGCTCCTTAGTCTTCCGA"
    @test NONCANONICAL(iterres) === mer"TCGGAAGACTAAGGAGCCTATCAGATATGTC"
end

include("MerCount.jl")
include("Counters.jl")

end # module