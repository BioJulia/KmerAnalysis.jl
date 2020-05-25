@testset "Counters" begin
    @testset "SerialMem" begin
        seq = LongSequence{DNAAlphabet{2}}("ACTGAAGAGC")
        @test serial_mem(DNAMer{4}, NONCANONICAL)(seq) == [
            MerCount{DNAMer{4}}(mer"AAGA",1), MerCount{DNAMer{4}}(mer"ACTG",1),
            MerCount{DNAMer{4}}(mer"AGAG",1), MerCount{DNAMer{4}}(mer"CTGA",1),
            MerCount{DNAMer{4}}(mer"GAAG",1), MerCount{DNAMer{4}}(mer"GAGC",1),
            MerCount{DNAMer{4}}(mer"TGAA",1)
        ]
    end
    
end