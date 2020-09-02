@testset "systems.jl" begin
    chainA = homopolymer_chain()
    @test islinearchain(chainA.architecture) == true
    @test ischarged(chainA.charged) == false

    chainAB = diblock_chain()
    @test islinearchain(chainAB.architecture) == true
    @test ischarged(chainAB.charged) == false
    @test chainAB.blocks[1].E2 == chainAB.blocks[2].E2
    @test chainAB.blocks[1].E1 != chainAB.blocks[2].E1
end