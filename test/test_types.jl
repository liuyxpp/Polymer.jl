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

    smol = SmallMolecule(:S)
    @test smol.b == 1.0
    @test smol.M == 1.0
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
    c = PolymerComponent(:AB, [A,B])
    @test islinearchain(c.architecture) == true
    @test ischarged(c.charged) == false
    @test c.α == 1.0
    @test c.ϕ == 1.0

    sc = SmallMoleculeComponent(:S)
    @test sc.α == 0.01
    @test sc.ϕ == 0.0
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
    chain = PolymerComponent(:AB, [A,B])
    solvent = SmallMoleculeComponent(:S)
    ABS = PolymerSystem([chain, solvent])
    @test ABS.dim == D1()
    @test isconfined(ABS.confinement) == false
    @test ABS.C == 1.0

    @test multicomponent(ABS) == true
    @test ncomponents(ABS) == 2
    @test species(chain) == [:A, :B]
    @test nspecies(chain) == 2
    @test species(solvent) == [:S]
    @test nspecies(solvent) == 1
    @test species(ABS) == [:A, :B, :S]
    @test nspecies(ABS) == 3

    @test systemtype(ABS) == PolymerSolution()
end

@testset "types.jl: Polymer System AB+A" begin
    sA = KuhnSegment(:A)
    sB = KuhnSegment(:B)
    e1A = FreeEnd(:A1)
    e2A = BranchPoint(:A2)
    e1B = FreeEnd(:B1)
    e2B = BranchPoint(:B2)
    A = PolymerBlock(:A, sA, 0.4, e1A, e2A)
    B = PolymerBlock(:B, sB, 0.6, e1B, e2B)
    chainAB = PolymerComponent(:AB, [A,B]; ϕ=0.1)
    e1hA = FreeEnd(:hA1)
    e2hA = FreeEnd(:hA2)
    hA = PolymerBlock(:hA, sA, 1.0, e1hA, e2hA)
    chainA = PolymerComponent(:hA, [hA]; α=0.5, ϕ=0.9)
    AB_A = PolymerSystem([chainAB, chainA])

    @test multicomponent(AB_A) == true
    @test ncomponents(AB_A) == 2
    @test species(chainAB) == [:A, :B]
    @test nspecies(chainAB) == 2
    @test species(chainA) == [:A]
    @test nspecies(chainA) == 1
    @test species(AB_A) == [:A, :B]
    @test nspecies(AB_A) == 2

    @test systemtype(AB_A) == PolymerBlend()
end

