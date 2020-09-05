@testset "types.jl: Confinement" begin
    @test isconfined(BulkConfinement()) == false
    @test isconfined(BoxConfinement()) == true
    @test isconfined(SlabConfinement()) == true
    @test isconfined(DiskConfinement()) == true
    @test isconfined(SphereConfinement()) == true
    @test isconfined(CylinderConfinement()) == true
end

@testset "types.jl: Charge" begin
    @test ischarged(Neutral()) == false
    @test ischarged(SmearedCharge()) == true
    @test ischarged(DiscreteCharge()) == true
end

@testset "types.jl: Polymer" begin
    @test islinearchain(LinearArchitecture()) == true
    @test islinearchain(StarArchitecture()) == false
    @test islinearchain(CombArchitecture()) == false
    @test islinearchain(RingArchitecture()) == false
end

@testset "types.jl: Specie" begin
    segment = KuhnSegment(:A)
    @test segment.b == 1.0
    @test segment.M == 1.0
end

@testset "types.jl: Block" begin
    sA = KuhnSegment(:A)
    e1A = FreeEnd(:A1)
    e2A = FreeEnd(:A2)
    A = PolymerBlock(:A, sA, 1.0, e1A, e2A)
    @test A.f == 1.0
end

@testset "types.jl: Component" begin
    sA = KuhnSegment(:A)
    sB = KuhnSegment(:B)
    e1A = FreeEnd(:A1)
    e2A = BranchPoint(:A2)
    e1B = FreeEnd(:B1)
    e2B = BranchPoint(:B2)
    A = PolymerBlock(:A, sA, 0.4, e1A, e2A)
    B = PolymerBlock(:B, sB, 0.6, e1B, e2B)
    c = BlockCopolymer(:AB, [A,B])
    @test islinearchain(c.architecture) == true
    @test ischarged(c.charged) == false
    pc = Component(c)
    @test pc.α == 1.0
    @test pc.ϕ == 1.0

    smol = SmallMolecule(:S)
    @test smol.b == 1.0
    @test smol.M == 1.0
    sc = Component(smol)
    @test sc.α == 1.0
    @test sc.ϕ == 1.0
end

@testset "types.jl: Polymer System AB+S" begin
    sA = KuhnSegment(:A)
    sB = KuhnSegment(:B)
    e1A = FreeEnd(:A1)
    e2A = BranchPoint(:A2)
    e1B = FreeEnd(:B1)
    e2B = BranchPoint(:B2)
    A = PolymerBlock(:A, sA, 0.4, e1A, e2A)
    B = PolymerBlock(:B, sB, 0.6, e1B, e2B)
    chain = BlockCopolymer(:AB, [A,B])
    smol = SmallMolecule(:S)
    polymer = Component(chain, 1.0, 0.1)
    solvent = Component(smol, 0.01, 0.9)
    χN_map = Dict(
        Set([:A,:B])=>20.0,
        Set([:A,:S])=>80.0,
        Set([:B,:S])=>120.0)
    ABS = PolymerSystem([polymer, solvent]; χN_map=χN_map)
    @test isconfined(ABS.confinement) == false
    @test ABS.C == 1.0

    @test multicomponent(ABS) == true
    @test ncomponents(ABS) == 2
    @test species(chain) == [:A, :B]
    @test nspecies(chain) == 2
    @test species(smol) == [:S]
    @test nspecies(smol) == 1
    @test species(polymer) == species(chain)
    @test nspecies(polymer) == nspecies(chain)
    @test species(solvent) == species(smol)
    @test nspecies(solvent) == nspecies(smol)
    @test species(ABS) == [:A, :B, :S]
    @test nspecies(ABS) == 3

    @test systemtype(ABS) == PolymerSolution()
end

@testset "types.jl: Polymer System AB+A" begin
    sA = KuhnSegment(:A)
    sB = KuhnSegment(:B)
    eA = FreeEnd(:A1)
    eAB = BranchPoint(:AB)
    eB = FreeEnd(:B1)
    A = PolymerBlock(:A, sA, 0.4, eA, eAB)
    B = PolymerBlock(:B, sB, 0.6, eB, eAB)
    chainAB = BlockCopolymer(:AB, [A,B])
    hA = PolymerBlock(:hA, sA, 1.0, FreeEnd(:hA1), FreeEnd(:hA2))
    chainA = BlockCopolymer(:hA, [hA])
    polymerAB = Component(chainAB, 1.0, 0.5)
    polymerA = Component(chainA, 0.5, 0.5)
    χN_map = Dict(Set([:A,:B])=>20.0)
    AB_A = PolymerSystem([polymerAB, polymerA]; χN_map=χN_map)

    @test multicomponent(AB_A) == true
    @test ncomponents(AB_A) == 2
    @test species(chainAB) == [:A, :B]
    @test nspecies(chainAB) == 2
    @test species(chainA) == [:A]
    @test nspecies(chainA) == 1
    @test species(polymerAB) == species(chainAB)
    @test nspecies(polymerAB) == nspecies(chainAB)
    @test species(polymerA) == species(chainA)
    @test nspecies(polymerA) == nspecies(chainA)
    @test species(AB_A) == [:A, :B]
    @test nspecies(AB_A) == 2

    @test systemtype(AB_A) == PolymerBlend()
end

