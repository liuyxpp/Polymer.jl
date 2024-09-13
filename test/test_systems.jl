@testset "systems.jl: Components" begin
    chainA = homopolymer_chain()
    @test length(chainA.blocks) == 1

    chainAB = diblock_chain()
    @test length(chainAB.blocks) == 2
    @test chainAB.blocks[1].E2 == chainAB.blocks[2].E2
    @test chainAB.blocks[1].E1 != chainAB.blocks[2].E1

    chainABC = linearABC()
    @test chainABC isa BlockCopolymer

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

    abc = ABC_system()
    @test abc isa PolymerSystem

    a_b = A_B_system()
    @test a_b isa PolymerSystem

    aba = AB_A_system()
    @test systemtype(aba) == PolymerBlend()
    @test multicomponent(aba) == true
    @test ncomponents(aba) == 2
    @test species(aba) == [:A, :B]
    @test nspecies(aba) == 2

    a_ab = A_AB_system()
    @test a_ab isa PolymerSystem

    ab_a_b = AB_A_B_system()
    @test ab_a_b isa PolymerSystem

    a_b_ab = A_B_AB_system()
    @test a_b_ab isa PolymerSystem

    abs = AB_S_system()
    @test systemtype(abs) == PolymerSolution()
    @test multicomponent(abs) == true
    @test ncomponents(abs) == 2
    @test species(abs) == [:A, :B, :S]
    @test nspecies(abs) == 3

    a_b_s = A_B_S_system()
    @test a_b_s isa PolymerSystem

    a_b_s1_s2 = A_B_S1_S2_system()
    @test a_b_s1_s2 isa PolymerSystem
end