@testset "make.jl: make" begin
    configfile = joinpath(@__DIR__, "ABS.yml")
    config = load_config(configfile, top="system")
    abs = make(config)
    @test systemtype(abs) == PolymerSolution()
end

@testset "make.jl: constructors" begin
    abs = AB_S_system()
    config = to_config(abs)

    sp1 = KuhnSegment(config.species[1])
    @test sp1.label == :A
    sp2 = KuhnSegment(config.species[2])
    @test sp2.label == :B
    sp3 = SmallMolecule(config.species[3])
    @test sp3.label == :S
    sps = [sp1, sp2, sp3]
    blockA = PolymerBlock(config.components[1].blocks[1], sps)
    @test blockA.label == :A
    blockB = PolymerBlock(config.components[1].blocks[2], sps)
    @test blockB.label == :B
    bcp = BlockCopolymer(config.components[1], sps)
    @test bcp.label == :AB
    smol = SmallMolecule(config.components[2], sps)
    @test smol.label == :S
    c1 = Component(config.components[1], sps)
    @test c1.molecule isa BlockCopolymer
    c2 = Component(config.components[2], sps)
    @test c2.molecule isa SmallMolecule
    abs1 = PolymerSystem(config)
    @test systemtype(abs) == PolymerSolution()
end