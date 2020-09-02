@testset "systems.jl: Components" begin
    chainA = homopolymer_chain()
    @test islinearchain(chainA.architecture) == true
    @test ischarged(chainA.charged) == false

    chainAB = diblock_chain()
    @test islinearchain(chainAB.architecture) == true
    @test ischarged(chainAB.charged) == false
    @test chainAB.blocks[1].E2 == chainAB.blocks[2].E2
    @test chainAB.blocks[1].E1 != chainAB.blocks[2].E1

    sol = solvent()
    @test sol.label == :S
end

@testset "systems.jl: Systems" begin
    ab = AB_system()
    @test systemtype(ab) == NeatPolymer()
    @test multicomponent(ab) == false
    @test ncomponents(ab) == 1
    @test species(ab) == [:A, :B]
    @test nspecies(ab) == 2

    aba = AB_A_system()
    @test systemtype(aba) == PolymerBlend()
    @test multicomponent(aba) == true
    @test ncomponents(aba) == 2
    @test species(aba) == [:A, :B]
    @test nspecies(aba) == 2

    abs = AB_S_system()
    @test systemtype(abs) == PolymerSolution()
    @test multicomponent(abs) == true
    @test ncomponents(abs) == 2
    @test species(abs) == [:A, :B, :S]
    @test nspecies(abs) == 3
end